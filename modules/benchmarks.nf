#!/bin/bash


// this is the last step of the pipeline and uses the generated trace.txt to benchmark the tools against each other


process benchmarks {

    // need correct container
    container "python_vis"
    publishDir "human_chr_21/results/benchmarks", mode: "symlink"

    input:
    path analysis_complete
    path script
    path trace

    output:
    path "benchmarks/"

    script:
    """
    Rscript $script $trace
    """
}