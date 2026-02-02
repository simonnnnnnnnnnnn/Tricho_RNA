#!/bin/bash

// this process is rather unusual as it serves to mark the end of the analysis, therefore it lacks some of the usual parts like an input
// this is done so that a following process can utilize the trace.txt created by nextflow

process finish_pipeline {

    publishDir "human_chr_21/results", mode: "symlink"

    output:
    path "analysis_complete.done", emit: done

    script:
    """
    touch analysis_complete.done
    """
}