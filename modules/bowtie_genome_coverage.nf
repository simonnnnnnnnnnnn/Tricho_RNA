#!/bin/bash

process bowtie_genome_coverage{

    container "python_vis"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"

    input:
    path bedgraph
    path script
    val experiment

    output:
    path "*.png", emit: png

    script:
    """
    python3 $script star $bedgraph ${bedgraph.simpleName}
    """
}