#!/bin/bash

// this process processes the results of STAR

process samtools_star {

    container "community.wave.seqera.io/library/samtools@sha256:5e5188b3dacf0117c4d4db891f2b9f4873d93288b4cfb26c1fba43eb7d241a4c"
    publishDir "${experiment}/results/star", mode: "symlink"

    input:
    path bamfile
    val experiment

    output:
    path "${bamfile.simpleName}_marked.bam", emit: marked_bam

    script:
    """
    samtools sort -n -o sorted_by_query.bam $bamfile
    samtools fixmate -m sorted_by_query.bam fixmate.bam
    samtools sort -o sorted_by_coordinate.bam fixmate.bam
    samtools index sorted_by_coordinate.bam
    samtools markdup -s sorted_by_coordinate.bam ${bamfile.simpleName}_marked.bam
    """

}