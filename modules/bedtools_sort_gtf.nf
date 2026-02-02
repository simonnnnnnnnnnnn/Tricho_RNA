#!/bin/bash

process bedtools_sort_gtf{

    container "community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4"
    publishDir "${experiment}/results", mode: "symlink"

    input:
    path unsorted_gtf
    path genome_txt
    val experiment

    output:
    path "${unsorted_gtf.simpleName}_sorted.gtf", emit: sorted_gtf

    script:
    """
    bedtools sort -g $genome_txt -i $unsorted_gtf > ${unsorted_gtf.simpleName}_sorted.gtf
    """
}