#!/bin/bash

process dge_deseq2_salmon {

    container "r_tricho"
    publishDir "${experiment}/results/salmon", mode: "symlink"

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
    Rscript $deseq2 $location_mapping salmon $infile ${files.join(' ')}
    """
}