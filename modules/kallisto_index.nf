#!/bin/bash

// this process exists solely to produce and index for kallisto

process kallisto_index{

    container "quay.io/biocontainers/kallisto@sha256:f2dc85d6d55e1c3bfdc437a7738b1217d0f81ee9c1c62013b5dd3b867713d482"
    publishDir "${experiment}/results/kallisto", mode: "symlink"

    input:
    path fasta
    val experiment

    output:
    path "index_kallisto_human.idx", emit: idx

    script:
    """
    kallisto index -i index_kallisto_human.idx $fasta
    """
}