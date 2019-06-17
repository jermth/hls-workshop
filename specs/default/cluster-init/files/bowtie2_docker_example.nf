/*
 * Microsoft Azure Genomics Workshop
 * Introduction to Nextflow
 */
params.reads = "$baseDir/NA12787/SRR622461_{1,2}.fastq.gz"
params.bowtie2_index = "$baseDir/bowtie2_index/grch38/grch38_1kgmaj"
params.number_of_reads = 1000
params.outdir = "results"

println ""
println "=========================================================="
println " Microsoft Azure Genomics Workshop "
println "=========================================================="
println "reads                          : ${params.reads}"
println "Number of reads to process:    : ${params.number_of_reads}"
println "bowtie2_index                  : ${params.bowtie2_index}"
println "Resultdir                      : ${params.outdir}"


/*
 * Run Bowtie2
 */
bowtie2_index_location = file(params.bowtie2_index)
number_of_reads = params.number_of_reads
 
// read_pairs = Channel.fromFilePairs(params.reads, flat: true)

Channel 
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
    .into { read_pairs_bowtie2_ch; read_pairs_samtools_ch; read_pairs_bcftools_ch } 

process runBowtie2 {
    
    container 'biocontainers/bowtie2:v2.3.1_cv1'

    publishDir params.outdir, mode:'copy'

    input:
    file bowtie2_index from bowtie2_index_location
    set pair_id, file(read1), file(read2) from read_pairs_bowtie2_ch
     
    output:
    file ("${pair_id}.sam") into bowtie2_sam_output
       
    """
    bowtie2 -t -x ${bowtie2_index_location} -p ${task.cpus} -U ${read1},${read2} -u ${number_of_reads} -S ${pair_id}.sam
    """
}

process samToBam {

    container 'biocontainers/samtools:v1.7.0_cv4'

    publishDir params.outdir, mode:'copy'

    input:
    set sample_id, file(read1), file(read2) from read_pairs_bcftools_ch
    file ("${sample_id}.sam") from bowtie2_sam_output
    
    output:
    file ("${sample_id}.bam") into samtools_bam_output

    """
    samtools view -bS ${sample_id}.sam | samtools sort > ${sample_id}.bam
    """

}

process bamToVCF {

    container 'fredhutch/bcftools:1.9'

    publishDir params.outdir, mode:'copy'

    input:
    set sample_id, file(read1), file(read2) from read_pairs_samtools_ch
    file ("${sample_id}.bam") from samtools_bam_output
    
    output:
    file ("${sample_id}.vcf") into bcftools_vcf_output

    """
    bcftools mpileup --no-reference ${sample_id}.bam > ${sample_id}.vcf
    """

}
