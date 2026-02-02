#!/bin/bash

// this process does not align (or pseudo-align) anything --> an index for salmon is created

process salmon_index{

    container "community.wave.seqera.io/library/salmon@sha256:b4519ea6d76868516e8c545fd709bf900638cb4b9130206730e09038b2e9e274"
    publishDir "${experiment}/results/salmon", mode: "symlink"

    input:
    path fasta
    val experiment

    output:
    path "index_salmon_human/", emit: idx

    script:
    """
    salmon index -t $fasta -i index_salmon_human k 31
    """
}