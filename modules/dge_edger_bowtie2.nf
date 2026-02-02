#!/bin/bash

process dge_edger_bowtie2 {

    container "r_tricho"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"

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
    Rscript $edger $location_mapping bowtie2 $infile ${files.join(' ')}
    """
}