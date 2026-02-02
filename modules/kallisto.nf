#!/bin/bash


// this process ONLY DOES THE ALIGNMENT
// the indexing that has to take place beforehand is done separately {


process kallisto{

    container "quay.io/biocontainers/kallisto@sha256:f2dc85d6d55e1c3bfdc437a7738b1217d0f81ee9c1c62013b5dd3b867713d482"
    publishDir "${experiment}/results/kallisto", mode: "symlink"

    input:
        tuple path(forward), path(reverse)
        path index //usually with .idx
        val experiment

    output:
        path "${forward.simpleName}/${forward.simpleName}_abundance.tsv", emit: tsv
        path "${forward.simpleName}/${forward.simpleName}_run_info.json", emit:json

    // das sollte die jetzt immer paarweise bearbeiten
    // potential problem: publishDir specifies outdir, but kallisto itself wants that info as well --> hope they dont cause conflicts
    script:
    """
    mkdir -p ${forward.simpleName}
    kallisto quant --index=$index --output-dir=${forward.simpleName} $forward $reverse
    mv ${forward.simpleName}/abundance.tsv ${forward.simpleName}/${forward.simpleName}_abundance.tsv
    mv ${forward.simpleName}/run_info.json ${forward.simpleName}/${forward.simpleName}_run_info.json
    """
}