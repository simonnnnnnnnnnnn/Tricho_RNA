#!/bin/bash

process create_ref_transcriptome{

    container "quay.io/biocontainers/gffread:0.12.7--h077b44d_6"
    publishDir "${experiment}/results", mode: "symlink"

    input:
    val experiment
    path annotation
    path genome_fasta

    output:
    path "reference_transcriptome.fasta", emit: reference_transcriptome


    script:
    """
    gffread $annotation -g $genome_fasta -w reference_transcriptome.fasta
    """
}