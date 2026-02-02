#!/bin/bash

// featureCounts is part of the Subreads package --> the docker image can therefore only be found under that name

process featurecounts_bowtie2{

    container "community.wave.seqera.io/library/subread@sha256:4b5569b45ab8d6f106b69cf5f682c34fce76e36871c8553f899eb54da58deb48"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"

    input:
    path bamfile
    path gtf
    val experiment

    output:
    path "${bamfile.simpleName}_counts.txt", emit: counts

    script:
    """
    featureCounts -p --ignoreDup -a $gtf -o ${bamfile.simpleName}_counts.txt $bamfile
    """
}