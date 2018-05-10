#!/usr/bin/env nextflow
/*
 * vim: syntax=groovy
 * -*- mode: groovy;-*-
 *
 ===============================================================================
 ViPR: Viral amplicon analysis and intrahost / low-frequency variant calling
 ===============================================================================
 # Homepage / Documentation
 https://github.com/nf-core/vipr
 # Authors
 Andreas Wilm <wilma@gis.a-star.edu.sg>
 -------------------------------------------------------------------------------
 */


// Check that Nextflow version is up to date enough
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException("Nextflow version too old: ${workflow.nextflow.version} < $params.nf_required_version")
    }
} catch (all) {
  log.error "====================================================\n" +
            "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}



/* Input validation
 */
input_ref_fasta = file(params.ref_fasta)
if (!input_ref_fasta.exists()) 
  exit 1, "Missing input reference fasta file: ${input_ref_fasta}"

cont_fasta = file(params.cont_fasta)
if (!cont_fasta.exists()) 
  exit 1, "Missing contamination fasta file: ${cont_fasta}"

if(!params.skip_kraken) {
    kraken_db = file(params.kraken_db)
    if (!kraken_db.exists()) 
    exit 1, "Missing contamination fasta file: ${kraken_db}"
}

log.info "=================================================="
log.info " nf-core/vipr : Viral amplicon analysis v${params.version}"
log.info "=================================================="


def GetReadPair = { sk, rk ->
    // FIXME if files don't exist, their path might be relative to the input yaml
    // see https://gist.github.com/ysb33r/5804364
    tuple(file(params.samples[sk].readunits[rk]['fq1']),
          file(params.samples[sk].readunits[rk]['fq2']))
}


def GetReadUnitKeys = { sk ->
    params.samples[sk].readunits.keySet()
}


/* input channel for references used in decontamination step
 */
cont_fasta_ch = Channel.from(
    file(cont_fasta), file(cont_fasta + ".amb"), file(cont_fasta + ".ann"),
    file(cont_fasta + ".bwt"), file(cont_fasta + ".pac"),
    file(cont_fasta + ".sa")).toList()
//cont_fasta_ch.subscribe { println "cont_fasta_ch $it" }


/* FIXME allow other means of defining input, e.g. CSV.
 * input for fastq files. channel has sample name as key and all read pairs following 
 * see https://groups.google.com/forum/#!topic/nextflow/CF7Joh5xrkU
 */
sample_keys = params.samples.keySet()
println "List of samples: " +  sample_keys.join(", ")
Channel
    .from(sample_keys)
    .map { sk -> tuple(sk, GetReadUnitKeys(sk).collect{GetReadPair(sk, it)}.flatten()) }
    .set { fastq_ch }
//fastq_ch.subscribe { println "$it" }



/* Trim and combine read-pairs per sample
 */
process trim_and_combine {
    tag { "Preprocessing of " + reads.size()/2 + "  read pairs for " + sample_id }
    //publishDir "${params.outdir}/${sample_id}/reads/", mode: 'copy'

    input:
        set sample_id, file(reads) from fastq_ch
    output:
        set sample_id, file("${sample_id}_R1-trimmed.fastq.gz"), file("${sample_id}_R2-trimmed.fastq.gz") into trim_and_combine_ch
    script:
        """
        # loop over readunits in pairs per sample
        pairno=0
        echo ${reads.join(" ")} | xargs -n2 | while read fq1 fq2; do
            let pairno=pairno+1
            # note: don't make reads smaller than assembler kmer length
            skewer --quiet -t ${task.cpus} -m pe -q 3 -n -l 31 -z -o pair\${pairno}-skewer-out \$fq1 \$fq2;
            cat *-trimmed-pair1.fastq.gz >> ${sample_id}_R1-trimmed.fastq.gz; 
            cat *-trimmed-pair2.fastq.gz >> ${sample_id}_R2-trimmed.fastq.gz; 
            rm *-trimmed-pair[12].fastq.gz;
        done
        fastqc -t {task.cpus} ${sample_id}_R1-trimmed.fastq.gz ${sample_id}_R2-trimmed.fastq.gz;
        """
}


/* Decontaminate reads against host reference
 */
process decont {
    tag { "Decontaminating " + sample_id }
    publishDir "${params.outdir}/${sample_id}/reads/", mode: 'copy'
    
    input:
        set sample_id, file(fq1), file(fq2) from trim_and_combine_ch
        set file(cont_fasta), file(cont_amb), file(cont_ann), file(cont_bwt), \
            file(cont_pac), file(cont_sa) from cont_fasta_ch
    output:
        set sample_id, file("${sample_id}_trimmed_decont_1.fastq.gz"), file("${sample_id}_trimmed_decont_2.fastq.gz") into \
            fastq_for_tadpole, fastq_for_polish_assembly_ch, fastq_for_mapping_ch, fastq_for_kraken_ch
        set file("${sample_id}_trimmed_decont_1_fastqc.zip"), file("${sample_id}_trimmed_decont_2_fastqc.zip"), \
            file("${sample_id}_trimmed_decont_1_fastqc.html"), file("${sample_id}_trimmed_decont_2_fastqc.html") into fastqc_ch
    script:
        // bbduk should be faster but uses plenty of memory (>32GB for human)
        """
        decont.py -i ${fq1} ${fq2} -t ${task.cpus} -c 0.5 -r ${cont_fasta} -o ${sample_id}_trimmed_decont;
        # since this is the last fastqc processing step, let's run fastqc here
        fastqc -t {task.cpus} ${sample_id}_trimmed_decont_1.fastq.gz ${sample_id}_trimmed_decont_2.fastq.gz;
        """
}


/* Metagenomics classification: QC for sample purity
 */
if(!params.skip_kraken) {
    process kraken {
        tag { "Running Kraken on " + sample_id }
        publishDir "${params.outdir}/${sample_id}/", mode: 'copy'
        
        input:
            set sample_id, file(fq1), file(fq2) from fastq_for_kraken_ch
        output:
            file("${sample_id}_kraken.report")
        script:
            """
            kraken --threads ${task.cpus} --preload --db ${kraken_db} \
              -paired ${fq1} ${fq2} > kraken.out;
            # do not gzip! otherwise kraken-report happily runs (with some warnings) and produces rubbish results
            kraken-report --db ${kraken_db} kraken.out > ${sample_id}_kraken.report
            """
    }
}


/* Assembly of reads. only few programs can assemble at the depth
 * possible for such samples (>100k). tadpole gave consistently better
 * results (fewer errors/variants) than spades (with diginorm)
 */
process tadpole {
    tag { "Tadpole assembly of " + sample_id }
    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'
    
    input:
        set sample_id, file(fq1), file(fq2) from fastq_for_tadpole
    output:
        set sample_id, file("${sample_id}_contigs.fa") into contigs_ch    
    script:
        """
        tadpole.sh -Xmx10g threads=${task.cpus} in=${fq1} in2=${fq2} out=${sample_id}_contigs.fa,
        """
}


/* Orient contigs according to reference and fill gaps with reference
 */
process gap_fill_assembly {
    tag { "Orienting and gap filling contigs for " + sample_id }
    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'
    // capture special error code, telling us that we cannot proceed for valid reasons.
    // meaning of consequently missing downstream files needs to be reflected in docs
    errorStrategy = { task.exitStatus == 3 ? 'ignore' : 'terminate' }
    
    input:
        set sample_id, file(contigs_fa) from contigs_ch
        file(input_ref_fasta)    
    output:
        set sample_id, file("${sample_id}-gap-filled-assembly.fa"), file("${sample_id}-gap-filled-assembly.gaps.bed") \
            into gap_filled_assembly_ch
    script:
        """
        set +e;
        log=${sample_id}-gap-filled-assembly.log;
        simple_contig_joiner.py -c ${contigs_fa} -r ${input_ref_fasta} \
          -s "${sample_id}-gap-filled-assembly" -o ${sample_id}-gap-filled-assembly.fa \
          -b "${sample_id}-gap-filled-assembly.gaps.bed" >& \$log;

        rc=\$?
        if [ \$rc -ne 0 ]; then
            # nothing to join, means we cannot continue. so stop here.
            grep 'Nothing to join' \$log && exit 3;
            exit \$rc;
        fi
        """
}


/* Polish assembly by repeated mapping, variant calling and variant incorporation 
 */
process polish_assembly {
    tag { "Polishing assembly for " + sample_id }
    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'
    
    input:
        set sample_id, file(assembly_fa), file(assembly_gaps_bed), file(fq1), file(fq2) \
            from gap_filled_assembly_ch.join(fastq_for_polish_assembly_ch)
    output:
        set sample_id, file("${sample_id}_polished_assembly.fa") into polished_assembly_ch
    script:
        """
        # downsample to 1M reads to increase runtime
        seqtk sample -s 666 ${fq1} 1000000 | gzip > R1_ds.R1.fastq.gz;
        seqtk sample -s 666 ${fq2} 1000000 | gzip > R2_ds.R2.fastq.gz;
        polish_viral_ref.sh -t ${task.cpus} -1 R1_ds.R1.fastq.gz -2 R2_ds.R2.fastq.gz \
            -r ${assembly_fa} -o ${sample_id}_polished_assembly.fa
        """
}


/* Mapping against polished assembly
 */
process final_mapping {
    tag { "Mapping to polished assembly for " + sample_id }
    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    input:
        set sample_id, file(ref_fa), file(fq1), file(fq2) \
            from polished_assembly_ch.join(fastq_for_mapping_ch)
    output:
        set sample_id, file(ref_fa), file("${sample_id}.bam"), file("${sample_id}.bam.bai") into \
            final_mapping_for_vcf_ch, final_mapping_for_cov_ch
        set sample_id, "${sample_id}.bam.stats" into final_mapping_bamstats_ch
    script:
        """
        bwa index ${ref_fa};
        samtools faidx ${ref_fa};
        bwa mem -t ${task.cpus} ${ref_fa} ${fq1} ${fq2} | \
            lofreq viterbi -f ${ref_fa} - | \
            lofreq alnqual -u - ${ref_fa} | \
            lofreq indelqual --dindel -f ${ref_fa} - | \
            samtools sort -o ${sample_id}.bam -T ${sample_id}.final.tmp -;
        samtools index ${sample_id}.bam;
        samtools stats ${sample_id}.bam > ${sample_id}.bam.stats
        """        
}


/* Low frequency variant calling
 */
process var_calling {
    tag { "Final variant calling for " + sample_id }
    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    input:
        set sample_id, file(ref_fa), file(bam), file(bai) from final_mapping_for_vcf_ch
    output:
        set sample_id, file(ref_fa), file("${sample_id}.vcf.gz") into vcf_ch
    script:
        """
        samtools faidx ${ref_fa};
        lofreq call-parallel --pp-threads ${task.cpus} -f ${ref_fa} \
           -d 1000000 --call-indels -o ${sample_id}.vcf.gz ${bam}
        """        
}


/* Compute coverage 
*
 * FIXME: fast process, best joined with another process
 */
process genomecov {
    tag { "Genome coverage for " + sample_id }
    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    input:
        set sample_id, file(ref_fa), file(bam), file(bai) from final_mapping_for_cov_ch
    output:
        set sample_id, file("${sample_id}.cov.gz") into cov_ch
    script:
        """
        # note: -d is one-based. -dz is zero-based but only non-zero values, so less explicit.
        bedtools genomecov -d -ibam ${bam} | gzip > ${sample_id}.cov.gz;
        """
}


/* Plot coverage and variant AF. Also polish reference for upload
 */
process vipr_tools {
    tag { "Plotting AF vs. coverage and readying fasta for  " + sample_id }
    publishDir "${params.outdir}/${sample_id}/", mode: 'copy'

    input:
        set sample_id, file(cov), file(ref_fa), file(vcf) from cov_ch.join(vcf_ch)
    output:
        set sample_id, file("${sample_id}_af-vs-cov.html"), file("${sample_id}_0cov2N.fa")
    script:
        """
        vipr_af_vs_cov_html.py --vcf ${vcf} --cov ${cov} --plot ${sample_id}_af-vs-cov.html;
        vipr_gaps_to_n.py -i ${ref_fa} -c ${cov} > ${sample_id}_0cov2N.fa;
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


