#!/bin/bash

process diff_gene_expr_star{

    container "diff_gene_expr"
    publishDir "${experiment}/results/${aligner}", mode: "symlink"

    input:
    path files
    path edger
    path limma
    path deseq2
    val experiment
    // now all variables that go into the R script
    path location_mapping
    val aligner
    path infile

    output:
    path "EdgeR_*.csv", emit: edger_dge_1
    path "EdgeR_*.csv", emit: edger_dge_2
    path "Limma-voom_*.csv", emit: limma_dge_1
    path "Limma-voom_*.csv", emit: limma_dge_2
    path "voomplot.png", emit: voomplot
    path "DESeq2_*.csv", emit: deseq2_dge_1
    path "DESeq2_*.csv", emit: deseq2_dge_2

    script:
    """
    Rscript $edger $location_mapping $aligner $infile ${files.join(' ')}
    Rscript $limma $location_mapping $aligner $infile ${files.join(' ')}
    Rscript $deseq2 $location_mapping $aligner $infile ${files.join(' ')}
    """
}