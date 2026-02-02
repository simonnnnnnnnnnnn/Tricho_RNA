#!/bin/bash

process dge_edger_star {

    container "r_tricho"
    publishDir "${experiment}/results/star", mode: "symlink"

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
    Rscript $edger $location_mapping star $infile ${files.join(' ')}
    """
}