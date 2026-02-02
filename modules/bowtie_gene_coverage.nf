#!/bin/bash

process bowtie_gene_coverage {

    container "python_vis"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"

    input:
    path coverfile
    path script
    val experiment

    output:
    path "*.png", emit: png

    script:
    """
    python3 $script star $coverfile ${coverfile.simpleName}
    """
}