#!/bin/bash

// this process converts samfiles to their binary format bam
// stringtie (used to generate gene counts) takes bam as input

process samtools {

    container "community.wave.seqera.io/library/samtools@sha256:5e5188b3dacf0117c4d4db891f2b9f4873d93288b4cfb26c1fba43eb7d241a4c"
    publishDir "human_chr_21/results/samtools", mode: "symlink"

    input:
    path samfile

    output:
    path "${samfile.simpleName}.bam", emit: bam

    script:
    """
    samtools view -b $samfile | samtools sort -o ${samfile.simpleName}.bam
    """

    // samtools view does the actual conversion, samtools sort only sorts it (important for stringtie)

}