#!/bin/bash

// here the first commentline(s) is/are removed --> otherwise edger will breal

process clean_features_bowtie2 {

    publishDir "${experiment}/results/bowtie2", mode: "symlink"
    //awk 'NR==1 { \$(NF)="Counts"; print } NR>1 { print }' OFS="\\t" ${features.simpleName}_temp2.txt > ${features.simpleName}_cleaned.txt

    input:
    path features
    val experiment

    output:
    path "${features.simpleName}_cleaned.txt"

    script:
    """
    grep -v '^#' $features > ${features.simpleName}_temp1.txt
    sed 's/gene://g' ${features.simpleName}_temp1.txt > ${features.simpleName}_cleaned.txt
    rm ${features.simpleName}_temp1.txt
    """
}