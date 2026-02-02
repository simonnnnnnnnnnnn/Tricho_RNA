#!/bin/bash

process dge_edger_kallisto {

    container "r_tricho"
    publishDir "${experiment}/results/kallisto", mode: "symlink"

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
    Rscript $edger $location_mapping kallisto $infile ${files.join(' ')}
    """
}