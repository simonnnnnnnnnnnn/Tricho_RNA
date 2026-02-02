#!/bin/bash

// this process only creates the index for bowtie2

process bowtie2_index{

    container "community.wave.seqera.io/library/bowtie2@sha256:7e95aab41e539b4696750e91020b20b09e45778688ebac01b3e3d03388aa6704"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"

    input:
    path fasta
    val experiment

    output:
    path "index_bowtie/", emit: idx
    path "test.txt", emit: txt

    script:
    """
    set -euxo pipefail

    mkdir -p index_bowtie
    ls -lh $fasta
    bowtie2-build -f $fasta index_bowtie/index
    touch test.txt
    """
}