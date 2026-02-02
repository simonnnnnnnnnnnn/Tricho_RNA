#!/bin/bash

process dge_limma_kallisto {

    container "r_tricho"
    publishDir "${experiment}/results/kallisto", mode: "symlink"

    input:
    path files
    path limma
    val experiment
    path location_mapping
    path infile

    output:
    path "*.csv", emit: csvs
    path "voomplot.png", emit: voomplot

    script:
    """
    Rscript $limma $location_mapping kallisto $infile ${files.join(' ')}
    """
}