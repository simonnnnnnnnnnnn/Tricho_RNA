#!/bin/bash

process picard_validate {

    container "community.wave.seqera.io/library/picard@sha256:e269216786463d44f9d83a0d6e877b34bca2c7b4d35211b4b369fe98e39ef1a5"
    publishDir "human_chr_21/results/picard", mode: "symlink"

    input:
    path bamfile
    path fasta
    path dict

    output:
    path "validation_report_${bamfile.simpleName}.txt", emit: validation_report

    script:
    """
    cp $fasta ref.fa
    cp $dict ref.dict


    picard ValidateSamFile \
    -I $bamfile \
    -O validation_report_${bamfile.simpleName}.txt
    -R $fasta
    -MODE SUMMARY
    -VALIDATE_INDEX false \
    -VALIDATION_STRINGENCY LENIENT
    """
}