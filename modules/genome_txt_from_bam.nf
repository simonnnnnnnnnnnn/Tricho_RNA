#!/bin/bash

process genome_txt_from_bam{

    container "community.wave.seqera.io/library/samtools@sha256:5e5188b3dacf0117c4d4db891f2b9f4873d93288b4cfb26c1fba43eb7d241a4c"
    publishDir "${experiment}/results", mode: "symlink"

    input:
    path bamfile
    val experiment

    output:
    path "${bamfile.simpleName}_genome.txt", emit: gen_txt

    script:
    """
    samtools view -H $bamfile | grep '^@SQ' | \
        awk -F'\t' '{
        for(i=1;i<=NF;i++) {
            if (\$i ~ /^SN:/) {
                split(\$i, sn, ":")
                chr = sn[2]
            }
            if (\$i ~ /^LN:/) {
                split(\$i, ln, ":")
                len = ln[2]
            }
        }
        print chr "\\t" len
    }' > ${bamfile.simpleName}_genome.txt
    """
}