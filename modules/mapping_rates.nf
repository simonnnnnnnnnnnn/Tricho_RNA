#!/bin/bash

process mapping_rates{

    container "python_vis"
    publishDir "${experiment}/results", mode: "symlink"

    input:
    path script
    path num_samples
    path kallisto
    path salmon
    path bowtie
    path star
    val experiment

    output:
    path "*.png", emit: png

    script:
    """
    python3 $script $num_samples ${kallisto.join(' ')} ${salmon.join(' ')} ${bowtie.join(' ')} ${star.join(' ')}
    """
}