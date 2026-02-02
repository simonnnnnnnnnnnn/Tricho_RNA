#!/bin/bash

process bowtie2_stats {

    container "community.wave.seqera.io/library/samtools@sha256:5e5188b3dacf0117c4d4db891f2b9f4873d93288b4cfb26c1fba43eb7d241a4c"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"

    input:
    path bamfile
    val experiment

    output:
    path "${bamfile.simpleName}_stats.txt", emit: stats

    script:
    """
    samtools stats $bamfile > ${bamfile.simpleName}_stats.txt
    """
}