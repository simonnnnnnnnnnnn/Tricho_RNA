#!/bin/bash

// this process does the alignment with bowtie

process bowtie2 {

    container "community.wave.seqera.io/library/bowtie2@sha256:7e95aab41e539b4696750e91020b20b09e45778688ebac01b3e3d03388aa6704"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"

    input:
    tuple path(forward), path(reverse)
    path index
    val experiment

    output:
    path "*.sam", emit: sam

    script:
    """
    echo "running bowtie2 on $forward and $reverse"
    bowtie2 -x ${index}/index -1 $forward -2 $reverse -S ${forward.simpleName}.sam 2> bowtie2_errors.log
    cat bowtie2_errors.log
    """

    // next steps: samtools (sam --> bam), stringtie (make gene counts)
    // dann STAR (anderer aligner)
}