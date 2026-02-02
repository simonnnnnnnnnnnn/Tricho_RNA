#!/bin/bash

process bedtools_bowtie2{

    container "community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4"
    publishDir "${experiment}/results/bowtie", mode: "symlink"

    input:
    val experiment
    path bamfile
    path gtf

    output:
    path "*", emit: bedgraph
    path "*", emit: gene_coverage

    script:
    """
    bedtools genomecov -ibam $bamfile -bg
    bedtools coverage -a $gtf -b $bamfile
    """
}