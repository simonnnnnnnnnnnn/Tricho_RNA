#!/bin/bash

// this process starts an R script in a docker container, here the analysis is done with EdgeR

process diff_gene_expr_edger_kallisto {

    container "diff_gene_expr"
    publishDir "${experiment}/results/${aligner}", mode: "symlink"

    input:
    path files
    path edger_script
    path limma
    path deseq2
    val experiment
    // now all variables that go into the R script
    path location_mapping
    val aligner

    output:
    path "EdgeR_control_condition1.csv", emit: edger_control_condition1
    path "EdgeR_control_condition2.csv", emit: edger_control_condition2
    path "Limma-voom_control_condition1.csv", emit: limma_voom_control_condition1
    path "Limma-voom_control_condition2.csv", emit: limma_voom_control_condition2
    path "DESeq2_control_condition1.csv", emit: deseq2_control_condition1
    path "DESeq2_control_condition2.csv", emit: deseq2_control_condition2

    script:
    """
    Rscript $edger_script $location_mapping ${files.size()} $aligner ${files.join(' ')}
    Rscript $limma $location_mapping ${files.size()} $aligner ${files.join(' ')}
    Rscript $deseq2 $location_mapping ${files.size()} $aligner ${files.join(' ')}

    """
}