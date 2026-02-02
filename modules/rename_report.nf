#!/bin/bash

process rename_report {

    input:
    val run_id


    script:
    """
    mv report.html rna_report_${run_id}.html
    """
}