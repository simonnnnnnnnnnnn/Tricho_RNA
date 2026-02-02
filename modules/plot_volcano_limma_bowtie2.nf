#!/bin/bash

process plot_volcano_limma_bowtie2{

    container "r_tricho"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"

    input:
    //tuple path(csv), path(script), val(experiment)
    // csv is actually a list that comes in --> need to loop through that still
    path csv
    path script
    val experiment

    output:
    path "*.png", emit: png

    script:
    """
    for csv in ${csv.join(' ')}; do
        base=\$(basename \$csv .csv)
        #mkdir -p plots_\$base
        Rscript $script \$base \$csv
        ls -l *.png
    done
    """
}