#!/bin/bash

process dge_limma_salmon {

    container "r_tricho"
    publishDir "${experiment}/results/salmon", mode: "symlink"

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
    Rscript $limma $location_mapping salmon $infile ${files.join(' ')}
    """
}