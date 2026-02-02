#!/bin/bash

// this will generate a gtf file, but thats not what im after here
// with the -A flag i get a tsv file displaying the expression levels --> thats what i want

process stringtie {

    container "community.wave.seqera.io/library/stringtie@sha256:d3e1bd02a4dae3a2a7ca616e2f34892a01f0d9d41a0e80937f995aec93b9e5c8"
    publishDir "human_chr_21/results/stringtie", mode: "symlink"

    input:
    path bamfiles
    path gtf

    output:
    path "${bamfiles.simpleName}.tsv", emit: tsv
    path "${bamfiles.simpleName}.gtf", emit: gtf // i dont really need that but nice to have

    script:
    """
    stringtie $bamfiles -G $gtf -e -A ${bamfiles.simpleName}.tsv -o ${bamfiles.simpleName}.gtf
    """
}