#!/bin/bash

// this process calls a python script which calculates TPM values for all transcripts
// the script uses featureCounts results as input

process calculate_tpms_star {

    container "python_vis"
    publishDir "${experiment}/results/star", mode: "symlink"

    input:
    path feature_table
    path script
    val experiment

    output:
    path "*_TPMs.csv", emit: complete_table

    script:
    """
    python $script $feature_table
    """
}