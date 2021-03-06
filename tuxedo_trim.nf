#!/usr/bin/env nextflow

/*
 * Copyright (c) 2016-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'Tuxedo-NF'.
 *
 *   Tuxedo-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Tuxedo-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Tuxedo-NF.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Main Tuxedo-NF pipeline script
 *
 * @authors
 * Evan Floden <evanfloden@gmail.com>
 */

log.info "T U X E D O - N F  ~  version 0.2"
log.info "====================================="
log.info "annotation             : ${params.annotation}"
log.info "fq files               : ${params.reads}"
log.info "output                 : ${params.output}"
log.info "index                  : ${params.index}"
log.info "====================================="
log.info "\n"

/*
 * Input parameters validation
 */

annotation_file               = file(params.annotation)
index_file                    = file(params.index)
index_file1                   = index_file + ".1.ht2"
index_name                    = index_file.getFileName()
index_dir                     = index_file.getParent()


/*
 * validate and create a channel for annotation input files
 */

Channel
    .fromPath ( annotation_file )
    .into { annotations1; annotations2 }
/*
 * Create a channel for read files
 */

if ( params.reads ) {
    Channel
        .fromFilePairs( params.reads, size: -1)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" and use_sra is false}
        .set { read_files }
}


// GENOME INDEXING
// ===============

process premade_index {

    input:
    file(index_dir)
    val(index_name)

    output:
    set val(index_name), file("genomeindex") into genome_index

    script:
    //
    // Premade HISAT2 genome index
    //
    """
    mkdir genomeindex
    cp ${index_dir}/${index_name}.*.ht2 genomeindex/.
    """
}

process trim{
    maxForks 10
    publishDir = [path: "${params.output}/trims", mode: 'link', overwrite: 'true' ]
    tag "reads: $sample_id"

    input:
    set val(sample_id), file(reads) from read_files

    output:
    set val(sample_id), file("${sample_id}_*p_trim.fastq.gz") into fatq_maps

    script:
    //
    // TRIM FASTQ
    //
    def single = reads instanceof Path
    if( single )
        """
        mkdir -p ${params.output}/trims
        java -jar -Xms8000m -Xmx8000m /Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 8 ${sample_id}_1.fastq.gz ${sample_id}_1p_trim.fastq.gz ILLUMINACLIP:/Trimmomatic-0.32/adapters/TruSeq2-SE.fa:2:30:10  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:51 > ${params.output}/trims/${sample_id}.trim.log 2>&1
       """
   else
        """
        mkdir -p ${params.output}/trims
        java -jar -Xms8000m -Xmx8000m /Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 8 ${sample_id}_1.fastq.gz ${sample_id}_2.fastq.gz \
        ${sample_id}_1p_trim.fastq.gz ${sample_id}_1u_trim.fastq.gz ${sample_id}_2p_trim.fastq.gz ${sample_id}_2u_trim.fastq.gz ILLUMINACLIP:/Trimmomatic-0.32/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:51 \
        > ${params.output}/trims/${sample_id}.trim.log 2>&1
        """

}



process mapping {
    maxForks 10
    publishDir = [path: "${params.output}/mapped_sams", mode: 'link', overwrite: 'true' ]
    tag "reads: $sample_id"

    input:
    set val(index_name), file(index_dir) from genome_index
    //set val(sample_id), file(reads) from read_files
    set val(sample_id), file(input_files) from fatq_maps

    output:
    set val(sample_id), file("${sample_id}.sam") into hisat2_sams
    file("fastqc_${sample_id}_logs") into fastqc_ch
    file("${sample_id}.maplog")
    script:
    //
    // HISAT2 mapper
    //
    def single = input_files instanceof Path
    if( single )
        """
        mkdir fastqc_${sample_id}_logs
        fastqc -o fastqc_${sample_id}_logs -f fastq -q ${input_files}
        hisat2 -p 16 -x ${index_dir}/${index_name} -U ${input_files} -S ${sample_id}.sam 2>${sample_id}.maplog
        """
    else
        """
        mkdir fastqc_${sample_id}_logs
        fastqc -o fastqc_${sample_id}_logs -f fastq -q ${input_files}
        hisat2 -p 16 -x ${index_dir}/${index_name} -1 ${input_files[0]} -2 ${input_files[1]}  -S ${sample_id}.sam 2>${sample_id}.maplog
        """
}



process sam2bam {
    maxForks 10
    publishDir = [path: "${params.output}/mapped_bams", mode: 'link', overwrite: 'true' ]
    tag "sam2bam: $name"

    input:
    set val(name), file(sam) from hisat2_sams

    output:
    set val(name), file("${name}.bam") into hisat2_bams

    script:
    //
    // SAM to sorted BAM files
    //
    """
    samtools view -S -b ${sam} | samtools sort -o ${name}.bam -
    """
}

hisat2_bams.into { hisat2_bams1; hisat2_bams2 }

process transcript_abundance {
    maxForks 10
    publishDir = [path: "${params.output}/stringtie_abundances", mode: 'link', overwrite: 'true' ]
    tag "reads: $name"

    input:
    set val(name), file(bam) from hisat2_bams2
    file annotation_f from annotations1.first()

    output:
    file("${name}") into ballgown_data

    script:
    //
    // Estimate abundances of merged transcripts in each sample
    //
    """
    stringtie  -p 16  -G ${annotation_f} -o ${name}/${name}.gtf -A ${name}/${name}.gene_abund.tab ${bam}
    """
}

workflow.onComplete {
        println ( workflow.success ? "\nDone! Open the following results in your browser --> $params.output" : "Oops .. something went wrong" )
}

