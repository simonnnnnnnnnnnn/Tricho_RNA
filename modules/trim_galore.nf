#!/bin/bash

process trim_galore{

    container "community.wave.seqera.io/library/trim-galore@sha256:4f00e7b2a09f3c8d8a9ce955120e177152fb1e56f63a2a6e186088b1250d9907"
    publishDir "${experiment}/results/trim_galore", mode: "symlink"

    input:
        tuple path(forward), path(reverse) //forward = 1, reverse = 2
        val experiment

    output:
        tuple path("*val_1.fq.gz"), path("*val_2.fq.gz"), emit: trimmed_reads // trim_galore always names the output _val_1-fq.gz
        path "*_trimming_report.txt", emit: trimming_reports // here the same thing with default names
        path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1 // again 1 = forward...
        path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

    script:
    """
    trim_galore --fastqc -paired ${forward} ${reverse}
    """
}
