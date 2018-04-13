#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
 *
 * AW:
 * - provide config -params-file
 */

input_ref_fasta = file(params.ref_fasta)

fq_pairs = Channel.fromFilePairs('fastq_input/*_R{1,2}*.fastq.gz', flat:true)
fq_pairs.into { fq_pairs_tadpole; fq_pairs_debug; fq_pairs_polish_assembly; fq_pairs_final_mapping }

fq_pairs_debug.subscribe { println "FastQ pairs: $it" }


process tadpole {
    tag { "Tadpole assembly of " + sample_id }
    cpus 1//runs in minutes on this type of data
    memory '10 GB'
    time = '1h'
    beforeScript 'source /mnt/projects/rpd/rc/init.2017-04'// FIXME can we define this globally?
    module 'bbmap/37.17'
    //publishDir '/data/chunks' FIXME where
    
    input:
        set sample_id, file(fq1), file(fq2) from fq_pairs_tadpole
    
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
            from gap_filled_assembly_ch.join(fq_pairs_polish_assembly)
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
        set sample_id, file(ref_fa), file(assembly_gaps_bed), file(fq1), file(fq2) \
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
    //publishDir '/data/chunks' FIXME where

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

