#!/bin/bash

// this process should not be necessary...but here we go

process clean_abundance{

    publishDir "${experiment}/results/kallisto", mode: "symlink"

    input:
    path abundance
    val experiment

    output:
    path "${abundance.simpleName}_cleaned.tsv", emit: cleaned_tsv

    script:
    """
    sed 's/transcript://g' $abundance > ${abundance.simpleName}_cleaned.tsv
    """
}