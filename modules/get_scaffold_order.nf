#!/bin/bash

process get_scaffold_order{

    container "community.wave.seqera.io/library/samtools@sha256:5e5188b3dacf0117c4d4db891f2b9f4873d93288b4cfb26c1fba43eb7d241a4c"
    publishDir "${experiment}/results", mode: "symlink"

    input:
    path bamfile
    val experiment

    output:
    path "${bamfile.simpleName}_genome.txt", emit: genome_txt
    path "${bamfile.simpleName}_scaffold_order.txt", emit: scaffold_order

    script:
    """
    samtools view -H ${bamfile} | grep '^@SQ' | awk -F'\t' '
    {
        chr=""; len="";
        for(i=1;i<=NF;i++) {
            if (\$i ~ /^SN:/) { split(\$i, sn, ":"); chr = sn[2] }
            if (\$i ~ /^LN:/) { split(\$i, ln, ":"); len = ln[2] }
        }
        if (chr != "" && len != "") {
            print chr "\\t" len >> "${bamfile.simpleName}_genome.txt"
            print chr >> "${bamfile.simpleName}_scaffold_order.txt"
        }
    }
    '
    """

}