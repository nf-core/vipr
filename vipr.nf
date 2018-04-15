#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
 *
 * AW:
 * - provide config -params-file
 */

/*
println "Full params debug output:"
params.each{ k, v -> println " ${k}:${v}" }
println ""
params.readunits.each{
    println " ${it.key}: ${it.value}"
}
*/


input_ref_fasta = file(params.ref_fasta)
if (!input_ref_fasta.exists()) 
  exit 1, "Missing input reference fasta file: ${input_ref_fasta}"

cont_fasta = file(params.cont_fasta)
if (!cont_fasta.exists()) 
  exit 1, "Missing contamination fasta file: ${cont_fasta}"
cont_fasta_ch = Channel.from(file(cont_fasta), file(cont_fasta + ".amb"), file(cont_fasta + ".ann"),
             file(cont_fasta + ".bwt"), file(cont_fasta + ".pac"),
             file(cont_fasta + ".sa")).toList()
//cont_fasta_ch.subscribe { println "cont_fasta_ch $it" }

if (params.publishdir == null)
  exit 1, "publishdir missing from params"

/*
fq_pairs = Channel.fromFilePairs('fastq_input/*_R{1,2}*.fastq.gz', flat:true)
fq_pairs.into { fq_pairs_tadpole; fq_pairs_debug; fq_pairs_polish_assembly; fq_pairs_final_mapping }
fq_pairs_debug.subscribe { println "FastQ pairs: $it" }
 */

/*
fqs_from_readunit = { sample_key, ru_key:
                     params.samples[sample_key].readunits[ru_key].fq1, params.samples[sample_key].readunits[ru_key].fq2}

Channel
    .from(sample_keys)
    .map { sk -> tuple(sk, tuple(rk, params.readunits[rk]['fq1'], params.readunits[rk]['fq2'])) }
    .groupTuple() 
    .set { fastq_ch }
 */



sample_keys = params.samples.keySet()
println "List of samples: " +  sample_keys.join(", ")
/*
readunits_keys = params.readunits.keySet()
println "List of readunits: " +  readunits_keys.join(", ")
*/
/* FIXME the readunits channel below breaks with empty fq1 or fq2
  (null fq2 shows up as input.2).  The following looks like a possible
  solution: https://github.com/nextflow-io/nextflow/issues/236 For now
 we assume and check for paired end

params.readunits.each{
    k,v -> assert v.containsKey('fq1') && v.containsKey('fq2')}
readunits = Channel
    .from ( readunits_keys )
    .map {
    [it, file(params.readunits[it]['fq1']), file(params.readunits[it]['fq2']) ] }
    //.subscribe onNext: { println "DEBUG " + it }, onComplete: { println 'Done' }
*/


/* old params/samples
Channel
    .from(readunits_keys)
    .map { rk -> tuple(params.samples.find{ it.value.contains(rk) }?.key,
                       tuple(rk, params.readunits[rk]['fq1'], params.readunits[rk]['fq2'])) }
    .groupTuple() 
    .set { fastq_ch }

fastq_ch.subscribe { println "$it" }
 */

//params.samples.keySet().each{ println it + " " + params.samples[it].readunits.keySet() }


/*
readunit_keys = []
params.samples.keySet().each{ readunit_keys += params.samples[it].readunits.keySet() }
 */


def GetReadPair = { sk, rk ->
    tuple(file(params.samples[sk].readunits[rk]['fq1']),
          file(params.samples[sk].readunits[rk]['fq2']))
}

def GetReadUnitKeys = { sk ->
    params.samples[sk].readunits.keySet()
}

/* channel with sample name as key and all read pairs following 
 * see https://groups.google.com/forum/#!topic/nextflow/CF7Joh5xrkU
*/
Channel
    .from(sample_keys)
    .map { sk -> tuple(sk, GetReadUnitKeys(sk).collect{GetReadPair(sk, it)}.flatten()) }
    .set { fastq_ch }

//fastq_ch.subscribe { println "$it" }


process trim_and_combine {
    tag { "Preprocessing of " + reads.size()/2 + "  read pairs for " + sample_key }
    cpus 4//runs in minutes on this type of data
    memory '6 GB'
    time = '2h'
    beforeScript 'source /mnt/projects/rpd/rc/init.2017-04'// FIXME can we define this globally?
    module 'skewer/0.2.2'
    //publishDir '/data/chunks' FIXME where

    input:
        set sample_key, file(reads) from fastq_ch
    output:
        set sample_key, file('R1-trimmed.fastq.gz'), file('R2-trimmed.fastq.gz') into trim_and_combine_ch
    script:
        """
        # loop over readunits in pairs per sample
        pairno=0
        echo ${reads.join(" ")} | xargs -n2 | while read fq1 fq2; do
            let pairno=pairno+1
            # note: don't make reads smaller than assembler kmer length
            skewer --quiet -t ${task.cpus} -m pe -q 3 -n -l 31 -z -o pair\${pairno}-skewer-out \$fq1 \$fq2;
            cat *-trimmed-pair1.fastq.gz >> R1-trimmed.fastq.gz; 
            cat *-trimmed-pair2.fastq.gz >> R2-trimmed.fastq.gz; 
            rm *-trimmed-pair[12].fastq.gz;
        done
        """
}

process decont {
    tag { "Decontaminating " + sample_key }
    cpus 4
    memory '8 GB'
    time = '6h'
    beforeScript 'source /mnt/projects/rpd/rc/init.2017-04'// FIXME can we define this globally?
    module 'decont/0.5'
    publishDir "${params.publishdir}/${sample_key}/", mode: 'copy'
    
    input:
        set sample_key, file(fq1), file(fq2) from trim_and_combine_ch
        set file(cont_fasta), file(cont_amb), file(cont_ann), file(cont_bwt), \
            file(cont_pac), file(cont_sa) from cont_fasta_ch
    output:
        set sample_key, file("decont_1.fastq.gz"), file("decont_2.fastq.gz") into \
          decont_ch, fq_pairs_polish_assembly_ch, fq_pairs_final_mapping
    script:
        """
        # bbduk should be faster but uses plenty of memory (>32GB for human)
        decont.py -i ${fq1} ${fq2} -t ${task.cpus} -c 0.5 -r ${cont_fasta} -o decont
        """
}


process tadpole {
    tag { "Tadpole assembly of " + sample_id }
    cpus 1//runs in minutes on this type of data
    memory '10 GB'
    time = '1h'
    beforeScript 'source /mnt/projects/rpd/rc/init.2017-04'// FIXME can we define this globally?
    module 'bbmap/37.17'
    //publishDir '/data/chunks' FIXME where
    
    input:
        set sample_id, file(fq1), file(fq2) from decont_ch
    
    output:
        set sample_id, file("${sample_id}_contigs.fa") into contigs_ch
    
    script:
        """
        tadpole.sh -Xmx10g threads=${task.cpus} in=${fq1} in2=${fq2} out=${sample_id}_contigs.fa,
        """
}


process gap_fill_assembly {
    tag { "Orienting and gap filling contigs for " + sample_id }
    cpus 1
    memory '1 GB'
    time = '1h'
    beforeScript 'source /mnt/projects/rpd/rc/init.2017-04'// FIXME can we define this globally?
    module 'simple-contig-joiner/91ab4d9'
    //publishDir '/data/chunks' FIXME where

    input:
        set sample_id, file(contigs_fa) from contigs_ch
        file(input_ref_fasta)
    
    output:
        set sample_id, file("${sample_id}-gap-filled-assembly.fa"), file("${sample_id}-gap-filled-assembly.gaps.bed") into gap_filled_assembly_ch
    script:
        """
        simple_contig_joiner.py -c ${contigs_fa} -r ${input_ref_fasta} \
          -s "${sample_id}-gap-filled-assembly" -o ${sample_id}-gap-filled-assembly.fa \
          -b "${sample_id}-gap-filled-assembly.gaps.bed"
        """
}


process polish_assembly {
    tag { "Polishing assembly for " + sample_id }
    cpus 4
    memory '4 GB'
    time = '2h'
    beforeScript 'source /mnt/projects/rpd/rc/init.2017-04'// FIXME can we define this globally?
    module 'polish-viral-ref/a7b09b9:seqtk/4feb6e8'
    //publishDir '/data/chunks' FIXME where

    input:
        set sample_id, file(assembly_fa), file(assembly_gaps_bed), file(fq1), file(fq2) \
            from gap_filled_assembly_ch.join(fq_pairs_polish_assembly_ch)
    output:
        set sample_id, file("polished_assembly.fa") into polished_assembly_ch
    script:
        """
        # downsample to 1M reads to increase runtime
        seqtk sample -s 666 ${fq1} 1000000 | gzip > R1_ds.R1.fastq.gz;
        seqtk sample -s 666 ${fq2} 1000000 | gzip > R2_ds.R2.fastq.gz;
        polish_viral_ref.sh -t ${task.cpus} -1 R1_ds.R1.fastq.gz -2 R2_ds.R2.fastq.gz \
            -r ${assembly_fa} -o polished_assembly.fa
        """
}


process final_mapping {
    tag { "Mapping to polished assembly for " + sample_id }
    cpus 4
    memory '4 GB'
    time = '2h'
    beforeScript 'source /mnt/projects/rpd/rc/init.2017-04'// FIXME can we define this globally?
    module "lofreq/2.1.3.1:bwa/0.7.15:samtools/1.3"
    //publishDir '/data/chunks' FIXME where

    input:
        set sample_id, file(ref_fa), file(fq1), file(fq2) \
            from polished_assembly_ch.join(fq_pairs_final_mapping)
    output:
        set sample_id, file(ref_fa), file("${sample_id}.cons.bam") into final_mapping_for_vcf_ch, final_mapping_for_cov_ch
    script:
        """
        bwa index ${ref_fa};
        samtools faidx ${ref_fa};
        bwa mem -t ${task.cpus} ${ref_fa} ${fq1} ${fq2} | \
            lofreq viterbi -f ${ref_fa} - | \
            lofreq alnqual -u - ${ref_fa} | \
            lofreq indelqual --dindel -f ${ref_fa} - | \
            samtools sort -o ${sample_id}.cons.bam -T ${sample_id}.final.tmp -
        """        
}


process var_calling {
    tag { "Final variant calling for " + sample_id }
    cpus 4
    memory '4 GB'
    time = '6h'
    beforeScript 'source /mnt/projects/rpd/rc/init.2017-04'// FIXME can we define this globally?
    module "lofreq/2.1.3.1:samtools/1.3"
    //publishDir '/data/chunks' FIXME where

    input:
        set sample_id, file(ref_fa), file(bam) from final_mapping_for_vcf_ch
    output:
        set sample_id, file("${sample_id}.vcf.gz")
    script:
        """
        samtools faidx ${ref_fa};
        samtools index ${bam};
        lofreq call-parallel --pp-threads ${task.cpus} -f ${ref_fa} \
           --call-indels -o ${sample_id}.vcf.gz ${bam}
        """        
}


process genomecov {
    tag { "Genome coverage for " + sample_id }
    cpus 4
    memory '2 GB'
    time = '2h'
    beforeScript 'source /mnt/projects/rpd/rc/init.2017-04'// FIXME can we define this globally?
    module "bedtools/2.25.0"
    publishDir "${params.publishdir}/${sample_id}/", mode: 'copy'

    input:
        set sample_id, file(ref_fa), file(bam) from final_mapping_for_cov_ch
    output:
        file("${sample_id}.cov.gz")
    script:
        """
        bedtools genomecov -d -ibam ${bam} | gzip > ${sample_id}.cov.gz;
        """
}
    
/* Introspection
 *
 * https://www.nextflow.io/docs/latest/metadata.html
 */
workflow.onComplete {
    println """Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Work Dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """
}
workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}


