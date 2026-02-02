#!/bin/bash

// this process contains the genome as well as the gene coverage
// since both use bedtools they are grouped here

process bowtie_calc_genome_cov{

    container "community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"


    input:
    val experiment
    path bamfile
    path gtf

    output:
    path "*.bedgraph", emit: bedgraph

    script:
    """
    bedtools genomecov -ibam $bamfile -bg > ${bamfile.simpleName}_genomecov.bedgraph
    """
}