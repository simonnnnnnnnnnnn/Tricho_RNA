#!/bin/bash

/* picard is used here for multiple things
1) QC --> another round of qc, this time AFTER the alignment
2) mark duplicates --> duplicated can occur but they influence following analyses --> mark them
3) sort bamfiles (when that has not taken place already) --> Bowtie2 results are already sorted by stringtie
4) finally a short validation of the bam files --> everything should be fine, but checking it is better

*/

process picard_alignment_summary_bowtie2{

    container "community.wave.seqera.io/library/picard@sha256:e269216786463d44f9d83a0d6e877b34bca2c7b4d35211b4b369fe98e39ef1a5"
    publishDir "${experiment}/results/bowtie2", mode: "symlink"

    input:
    path bamfile
    path fasta
    val experiment

    output:
    path "${bamfile.simpleName}_alignment_summary.txt", emit: summary

    script:
    """
    picard CollectAlignmentSummaryMetrics \
    R=$fasta \
    I=$bamfile \
    O=${bamfile.simpleName}_alignment_summary.txt
    """
}