#!/bin/bash

process rename_trace{

    input:
    val run_id


    script:
    """
    mv trace.txt rna_trace_${run_id}.txt
    """
}