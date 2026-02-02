#!/bin/bash

process plot_volcano_deseq2_star{

    container "r_tricho"
    publishDir "${experiment}/results/star", mode: "symlink"

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