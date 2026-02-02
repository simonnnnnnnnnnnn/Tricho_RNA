#!/bin/bash

process sort_gtf{

    publishDir "${experiment}/results", mode: "symlink"

    input:
    path gtf
    val experiment

    output:
    path "${gtf.simpleName}_sorted.gtf", emit: sorted_gtf

    script:
    """
    sort -k1,1 -k4,4n $gtf > ${gtf.simpleName}_sorted.gtf
    """
}