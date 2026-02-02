#!/bin/bash

process dge_deseq2_star {

    container "r_tricho"
    publishDir "${experiment}/results/star", mode: "symlink"

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
    Rscript $deseq2 $location_mapping star $infile ${files.join(' ')}
    """
}