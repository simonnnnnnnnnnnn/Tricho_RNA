#!/bin/bash

// this process does the alignment with STAR, the index for it are created by the process star_index

process star{

    container "community.wave.seqera.io/library/star@sha256:25a18a65561065970d2bee3cc481949b58711c31b48499dd4e86df18aa69e3a8"
    publishDir "${experiment}/results/star", mode: "symlink"

    input:
    tuple path(forward), path(reverse)
    path index
    val experiment

    output:
    path "${forward.simpleName}/${forward.simpleName}_Log.out", emit: main_log
    path "${forward.simpleName}/${forward.simpleName}_Log.progress.out", emit: progress_log
    path "${forward.simpleName}/${forward.simpleName}_Log.final.out", emit: final_log
    path "${forward.simpleName}/${forward.simpleName}_Aligned.sortedByCoord.out.bam", emit: bam
    path "${forward.simpleName}/${forward.simpleName}_ReadsPerGene.out.tab", emit: gene_counts


    // slash is super important when using --outFileNamePrefix --> otherwise it write to other directory...
    script:
    """
    mkdir -p ${forward.simpleName}
    STAR --genomeDir $index \
    --readFilesIn $forward $reverse \
    --readFilesCommand zcat \
    --outFileNamePrefix ${forward.simpleName}/${forward.simpleName}_ \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts
    """
}