#!/bin/bash

process clean_quant{

    publishDir "${experiment}/results/salmon", mode: "symlink"

    input:
    path quant
    val experiment

    output:
    path "${quant.simpleName}_cleaned.sf", emit: cleaned_sf

    script:
    """
    sed 's/transcript://g' $quant > ${quant.simpleName}_cleaned.sf
    """
}