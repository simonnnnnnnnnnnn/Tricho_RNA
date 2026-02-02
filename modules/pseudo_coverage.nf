#!/bin/bash

process pseudo_coverage{

    container "python_vis"
    publishDir "${experiment}/results", mode: "symlink"

    input:
    path kallistojson
    path salmonjson
    path script
    path info
    val experiment

    output:
    path "*.png", emit: png

    script:
    """
    python3 $script $info ${kallistojson.join(' ')} ${salmonjson.join(' ')}
    """
}