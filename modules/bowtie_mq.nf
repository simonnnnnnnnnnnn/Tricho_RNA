#!/bin/bash

process bowtie_mq{

    container "python_vis"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"

    input:
    path statsfile
    path script
    val experiment

    output:
    path "*.png", emit: png

    script:
    """
    python3 $script star $statsfile ${statsfile.simpleName}
    """
}