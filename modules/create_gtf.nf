#!/bin/bash

// this process converts a gff file into gtf --> only used when annotation_type == gff


process create_gtf{

    container "quay.io/biocontainers/gffread:0.12.7--h077b44d_6"
    publishDir "${experiment}/results", mode: "symlink"

    input:
    path annotation
    val experiment

    output:
    path "created_gtf.gtf", emit: gtf

    script:
    """
    gffread $annotation -T -o created_gtf.gtf
    """

}