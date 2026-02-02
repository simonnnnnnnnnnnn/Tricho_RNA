#!/bin/bash

// this process handles the quantification with salmon
// the index for it is created beforehand with the salmon_index module

process salmon {

    container "community.wave.seqera.io/library/salmon@sha256:b4519ea6d76868516e8c545fd709bf900638cb4b9130206730e09038b2e9e274"
    publishDir "${experiment}/results/salmon", mode: "symlink"

    input:
    tuple path(forward), path(reverse)
    path index
    val experiment

    output:
    path "${forward.simpleName}/${forward.simpleName}_quant.sf", emit: sf
    path "${forward.simpleName}/${forward.simpleName}_meta_info.json", emit: json

    script:
    """
    mkdir -p ${forward.simpleName}
    salmon quant -i $index -l A -1 $forward -2 $reverse --validateMappings -o ${forward.simpleName}
    mv ${forward.simpleName}/quant.sf ${forward.simpleName}/${forward.simpleName}_quant.sf
    mv ${forward.simpleName}/aux_info/meta_info.json ${forward.simpleName}/${forward.simpleName}_meta_info.json
    """

}