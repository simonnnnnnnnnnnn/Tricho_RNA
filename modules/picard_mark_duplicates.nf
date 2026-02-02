#!/bin/bash

process picard_mark_duplicates {

    container "broadinstitute/picard@sha256:e57588109b67df22681acc9f6bd5be38a68c142450ccfd2d4dfc5883421cf54f"
    publishDir "human_chr_21/results/picard", mode: "symlink"

    input:
    path bamfile

    output:
    path "${bamfile.simpleName}_marked.bam", emit: bam
    path "${bamfile.simpleName}_metrics_marked.txt", emit: txt

    script:
    """
    java -Xmx4g -jar /usr/picard/picard.jar MarkDuplicates \
    I="$bamfile"
    O="${bamfile.simpleName}_marked.bam" \
    M="${bamfile.simpleName}_metrics_marked.txt" \
    VALIDATION_STRINGENCY=SILENT \
    ASSUME_SORT_ORDER=coordinate
    """
}