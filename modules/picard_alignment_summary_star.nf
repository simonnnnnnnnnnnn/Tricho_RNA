#!/bin/bash

process picard_alignment_summary_star{

    container "community.wave.seqera.io/library/picard@sha256:e269216786463d44f9d83a0d6e877b34bca2c7b4d35211b4b369fe98e39ef1a5"
    publishDir "${experiment}/results/star", mode: "symlink"

    input:
    path bamfile
    path fasta
    val experiment

    output:
    path "${bamfile.simpleName}_alignment_summary.txt", emit: summary

    script:
    """
    picard CollectAlignmentSummaryMetrics \
    R=$fasta \
    I=$bamfile \
    O=${bamfile.simpleName}_alignment_summary.txt
    """
}