#!/bin/bash

// this process mirrors "samtools_for_star", here an additional step is included:
// Bowtie2 produces .sam as its output whilw star directly produces .bam
// therefore an additional conversion step is necessary

process sam_to_bam {

    container "community.wave.seqera.io/library/samtools@sha256:5e5188b3dacf0117c4d4db891f2b9f4873d93288b4cfb26c1fba43eb7d241a4c"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"

    input:
    path samfile
    val experiment

    output:
    path "${samfile.simpleName}.bam", emit: bam

    script:
    """
    samtools view -b $samfile > ${samfile.simpleName}.bam
    """
}