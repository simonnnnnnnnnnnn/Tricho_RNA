#!/bin/bash

process compare_tpms{

    container "python_vis"
    publishDir "${experiment}/results", mode: "symlink"

    input:
    path kallisto_tpms
    path salmon_tpms
    path bowtie2_tpms
    path star_tpms
    path script
    path info
    val experiment

    output:
    path "*.png", emit: png

    script:
    """
    python3 $script $info ${kallisto_tpms.join(' ')} ${salmon_tpms.join(' ')} ${bowtie2_tpms.join(' ')} ${star_tpms.join(' ')}
    """
}