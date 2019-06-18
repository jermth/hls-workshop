/*
 * Microsoft Azure Genomics Workshop
 * Introduction to Nextflow
 */

params.averedatadir = "/avere/data/genomics_workshop_data"
params.reads = "${params.averedatadir}/NA12787/SRR622461_{1,2}.fastq.gz"
params.bowtie2_index = "${params.averedatadir}/bowtie2_index/grch38/grch38_1kgmaj"
params.number_of_reads = 100000
params.outdir = "results"
params.prefix = "exercise4"

println ""
println "=========================================================="
println " Microsoft Azure Genomics Workshop "
println "=========================================================="
println "Avere Data Dir                 : ${params.averedatadir}"
println "reads                          : ${params.reads}"
println "Number of reads to process:    : ${params.number_of_reads}"
println "bowtie2_index                  : ${params.bowtie2_index}"
println "Resultdir                      : ${params.outdir}"
println "Results prefix                 : ${params.prefix}"



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
    
    publishDir params.outdir, mode:'copy'

    input:
    file bowtie2_index from bowtie2_index_location
    set pair_id, file(read1), file(read2) from read_pairs_bowtie2_ch
     
    output:
    file ("${pair_id}.${params.prefix}.sam") into bowtie2_sam_output
       
    """
    bowtie2 -t -x ${bowtie2_index_location} -p ${task.cpus} -U ${read1},${read2} -u ${number_of_reads} -S ${pair_id}.${params.prefix}.sam
    """
}

process samToBam {
    publishDir params.outdir, mode:'copy'

    input:
    set sample_id, file(read1), file(read2) from read_pairs_bcftools_ch
    file ("${sample_id}.${params.prefix}.sam") from bowtie2_sam_output
    
    output:
    file ("${sample_id}.${params.prefix}.bam") into samtools_bam_output

    """
    samtools view -bS ${sample_id}.${params.prefix}.sam | samtools sort > ${sample_id}.${params.prefix}.bam
    """

}

process bamToVCF {
    publishDir params.outdir, mode:'copy'

    input:
    set sample_id, file(read1), file(read2) from read_pairs_samtools_ch
    file ("${sample_id}.${params.prefix}.bam") from samtools_bam_output
    
    output:
    file ("${sample_id}.${params.prefix}.vcf") into bcftools_vcf_output

    """
    bcftools mpileup --no-reference ${sample_id}.${params.prefix}.bam > ${sample_id}.${params.prefix}.vcf
    """

}
