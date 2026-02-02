#!/bin/bash

// this only works with a fasta and a GTF FILE
// if only gff3 is available: convert to gtf first --> gffread can be used for this

process star_index{

    container "community.wave.seqera.io/library/star@sha256:25a18a65561065970d2bee3cc481949b58711c31b48499dd4e86df18aa69e3a8"
    publishDir "${experiment}/results/star", mode: "symlink"

    input:
    path fasta
    path gtf
    val experiment

    output:
    path "index_star/", emit: idx

    // need to see how i can decide how many threads / resources one container can have
    script:
    """
    STAR --runMode genomeGenerate --genomeDir index_star --genomeFastaFiles $fasta --sjdbGTFfile $gtf --sjdbOverhang 100
    """
}