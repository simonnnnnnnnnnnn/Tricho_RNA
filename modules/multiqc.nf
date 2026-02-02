#!/bin/bash

// multiqc is done to create a comprehensive report on quality, aggregated from fastqc, trim_galore and alignment logs (if possible)

process multiqc{

    container "community.wave.seqera.io/library/pip_multiqc@sha256:0ebb1d9605395a7df49ad0eb366b21f46afd96a5090376b0d8941cf5294a895a"
    publishDir "${experiment}/results/multiqc", mode: "symlink"

    input:
        path "*"
        val name_multi_report
        val experiment
    
    output:
        path "${name_multi_report}.html", emit: report
        path "${name_multi_report}_data", emit: data

    // this aggregates all reports and combines it into one file --> easy browsing (on the browser :D)
    script:
    """
    multiqc . -n ${name_multi_report}.html
    """
}