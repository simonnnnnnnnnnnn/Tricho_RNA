#!/bin/bash

process dge_edger_salmon{

    container "r_tricho"
    publishDir "${experiment}/results/salmon", mode: "symlink"

    input:
    path files
    path edger
    val experiment
    path location_mapping
    path infile

    output:
    path "*.csv", emit: csvs

    script:
    """
    Rscript $edger $location_mapping salmon $infile ${files.join(' ')}
    """
}