#!/bin/bash

// this process produces the tx2gene file necessary for differential gene expression

process make_tx2gene_from_gff{

    container "python_vis"
    publishDir "${experiment}/results", mode: "symlink"

    input:
    val experiment
    path annotation
    path python_script

    output:
    path "tx2gene.csv", emit: tx2gene

    script:
    """
    python3 $python_script $annotation
    """
}