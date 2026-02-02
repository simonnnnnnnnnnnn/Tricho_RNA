#!/bin/bash


// this processes sole purpose is to create dict for the validation

process picard_seq_dict {


    container "community.wave.seqera.io/library/picard@sha256:e269216786463d44f9d83a0d6e877b34bca2c7b4d35211b4b369fe98e39ef1a5"
    publishDir "human_chr_21/material", mode: "symlink"

    input:
    path fasta

    output:
    path "Homo_sapiens.GRCh38.dna.chromosome.21.dict", emit: dict

    script:
    """
    picard CreateSequenceDictionary \
    -R $fasta \
    -O Homo_sapiens.GRCh38.dna.chromosome.21.dict
    """
}