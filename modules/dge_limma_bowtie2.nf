#!/bin/bash

process dge_limma_bowtie2 {

    container "r_tricho"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"

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
    Rscript $limma $location_mapping bowtie2 $infile ${files.join(' ')}
    """
}