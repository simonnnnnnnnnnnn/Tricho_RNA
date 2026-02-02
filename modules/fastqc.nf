#!/bin/bash

process fastqc {

    container 'community.wave.seqera.io/library/trim-galore@sha256:4f00e7b2a09f3c8d8a9ce955120e177152fb1e56f63a2a6e186088b1250d9907'
    publishDir "${experiment}/results/fastqc", mode: "symlink" //there may be a better option

    input:
        path infile //this references the defined file in the params
        val experiment


    output:
        path "*_fastqc.zip", emit: zip //need to specify all outputs here, if anything gets put out that isnt specified the pipeline isnt happy
        path "*_fastqc.html", emit: html // using glob patterns is the easiest here


    // this is basically a bash script --> i wonder if i can call scripts in other languages from here and specify the output of those
    // so basicall running an R script just like fastqc is run here
    // this should work without braces...but they sure wont hurt
    script:
    """
    echo "running on ${infile}"
    fastqc ${infile}
    echo "check output"
    ls -lh
    """
}

