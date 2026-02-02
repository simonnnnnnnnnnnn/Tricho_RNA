#!/bin/bash

process sort_gtf_by_scaffold{

    publishDir "${experiment}/results", mode: "symlink"

    input:
    path gtf
    path scaffold_order
    val experiment

    output:
    path "${gtf.simpleName}_sorted.gtf", emit: sorted_gtf

    script:
    """
    awk 'BEGIN {
        while ((getline < "${scaffold_order}") > 0) {
            order[\$1] = ++i
        }
    }
    {
        if (\$0 ~ /^#/ || \$1 == "") {
            print \$0 > "tmp_gtf"
        } else {
            print order[\$1], \$0 > "tmp_gtf"
        }
    }'
    sort -k1,1n tmp_gtf | cut -f2- > ${gtf.simpleName}_sorted.gtf
    """
}