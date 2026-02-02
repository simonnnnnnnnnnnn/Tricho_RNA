#!/bin/bash

process picard_sort_sam{
    
    container "community.wave.seqera.io/library/picard@sha256:e269216786463d44f9d83a0d6e877b34bca2c7b4d35211b4b369fe98e39ef1a5"
    publishDir "human_chr_21/results/picard", mode: "symlink"

    input:
    path bamfile

    output:
    path "${bamfile.simpleName}_sorted.bam", emit: sorted

    script:
    """
    picard SortSam \
    I=$bamfile \
    O=${bamfile.simpleName}_sorted.bam \
    SORT_ORDER=coordinate
    """
}