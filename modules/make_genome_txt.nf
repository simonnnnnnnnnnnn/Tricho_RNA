#!/bin/bash

process make_genome_txt{

    publishDir "${experiment}/results", mode: "symlink"

    input:
    path gff
    val experiment

    output:
    path "${gff.simpleName}.txt", emit: genome_txt

    script:
    """
    grep '##sequence-region' $gff | awk '{print \$2"\t"\$4}' > ${gff.simpleName}.txt
    """
}