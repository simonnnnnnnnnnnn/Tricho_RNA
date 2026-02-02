#!/bin/bash

process dge_deseq2_bowtie2 {

    container "r_tricho"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"

    input:
    path files
    path deseq2
    val experiment
    path location_mapping
    path infile

    output:
    path "*.csv", emit: csvs

    script:
    """
    Rscript $deseq2 $location_mapping bowtie2 $infile ${files.join(' ')}
    """
}