#!/bin/bash

// this process contains the genome as well as the gene coverage
// since both use bedtools they are grouped here

process bowtie_calc_gene_cov{

    container "community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"


    input:
    val experiment
    path bamfile
    path gtf
    path txt

    output:
    path "*.txt", emit: gene_coverage

    script:
    """
    bedtools coverage -sorted -a $gtf -b $bamfile -g $txt > ${bamfile.simpleName}_coverage.txt
    """
}