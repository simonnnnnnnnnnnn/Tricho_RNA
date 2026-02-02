# RNASeq Pipeline

# Tricho_RNA


**Objective:** Creation of an RNAseq pipeline with Nextflow for _Trichoplax adhaerens_

## Setup 

The Pipeline was executed on a virtual machine running Ubuntu 22.04, the VM itself ran on SimpleVMs Servers.

### Formal Requirements

 - In order to run a Virtual Machine on the SimpleVM servers the project (and the resources) must be approved by de.NBI
 - next the user must create a private Key, this key is needed to log into the virtual machine(s) using ssh.
 - the provate key has to be downloaded, the file must have sufficiently restricted access, I made root the owner of this file using the following command:

```b̀ash
sudo chown root:root /path/to/my/private/key
```
 - the ssh command to access the virtual machine is displayed in [simplevm](https://simplevm.denbi.de/portal/webapp/#/profile) after configuring the virtual machine

### Configuring the Virtual Machine

The configuration of the virtual machine was entirely done on SimpleVMs Portal, I chose the "Flavor" de.NBI medium.
This Flavor came with the following specs:
 - 14VCPUs
 - 32 GB RAM
 - 50 GB Root Disk space

I chose Ubuntu 22.04 de.NBI as the Image (2025-05-06)

### setup of high-mem machine

Later on i had access to a high memory VM, making analyses with real data possible. All following steps regarding docker or other tools that needed to be installed to run the pipeline were performed for both machines. This machine was NOT set up by me, i only created a new nextflow environment with conda in order to not interfere with other users on the same machine:

```bash
conda create -n nextflowenv openjdk=21 git
conda activate nextflowenv
curl -s https://get.nextflow.io | bash
sudo chmod +x nextflow
sudo mv /home/ubuntu/nextflow /bin
confirm: nextflow -v # this step confirms the successful installation
```

### Customisation of the VM

In order to use the virtual machine I installed some software and ran some other commands. All commands dealing with the configuration of the VM are displayed below (ordered chronologically)

```b̀ash
mkdir project
nano test.txt
sudo apt install git
sudo apt install docker
sudo apt install openjdk-21-jdk
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv /home/ubuntu/nextflow /bin
nextflow info
sudo apt install unzip
sudo apt install tree
sudo apt install tabix
sudo apt install graphviz
```

tree and unzip are not strictly necessary to run the pipeline, but depending on the format of the reference data the may be necessary, or just informative like tree.

The commands to install nextflow follow [Nextflows Instructions](https://training.nextflow.io/latest/envsetup/02_local/#download-nextflow) for the local environment (2025-05-07).

### Transfer of Files to the VM and from the VM

To keep a local directory (where the code is written) in sync with a directory on the VM **rsync** is used. The following example shows how this was done:

```b̀ash
rsync -avu --progress -e "ssh -i ~/path/to/my/private/key.org_ecdsa -p 12345" ~/path/tomy/local/directory/ ubuntu@123.45.67.8:path/to/vm/directory/
```

To get the results from the Vm back to my local machine i used rsync as well, albeit with some changes to the flags, these were necessary as i used symlinks to save results instead of using "copy" mode. In order to get actual files not broken ones the inlusion of the "-L" flag was necessary:

```bash
rsync -avL --progress -e "ssh -i ~/path/to/my/private/key.org_ecdsa -p 12345" ubuntu@123.45.67.8:path/to/vm/directory/ ~/path/to/where/plots/should/go
```


### Nextflow

Nextflow itself requires the following the following software to run properly as is stated in [Nextflows Instructions](https://training.nextflow.io/latest/envsetup/02_local/#download-nextflow) (can be seen in the **Customisation of the VM** as well).

 - Docker
 - git
 - Java (versions 11 to 21)

### Docker

Docker of course requires images to run in its containers, the following table shows the program ran in a container, its corresponding image and the source for the image.
All of these images can be installed using this command:

```bash
docker pull name_docker_image
```

| program | image | repodigest |
| ----- | ----- | ---- |
| fastqc & trim_galore | community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18 | community.wave.seqera.io/library/trim-galore@sha256:4f00e7b2a09f3c8d8a9ce955120e177152fb1e56f63a2a6e186088b1250d9907 |
| multiqc | community.wave.seqera.io/library/pip_multiqc:ad8f247edb55897c | community.wave.seqera.io/library/pip_multiqc@sha256:0ebb1d9605395a7df49ad0eb366b21f46afd96a5090376b0d8941cf5294a895a |
| kallisto | quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1 | quay.io/biocontainers/kallisto@sha256:f2dc85d6d55e1c3bfdc437a7738b1217d0f81ee9c1c62013b5dd3b867713d482 |
| salmon | community.wave.seqera.io/library/salmon:1.10.3--fcd0755dd8abb423 | community.wave.seqera.io/library/salmon@sha256:b4519ea6d76868516e8c545fd709bf900638cb4b9130206730e09038b2e9e274 |
| bowtie2 | community.wave.seqera.io/library/bowtie2:2.5.4--d51920539234bea7 | community.wave.seqera.io/library/bowtie2@sha256:7e95aab41e539b4696750e91020b20b09e45778688ebac01b3e3d03388aa6704 |
| samtools | community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c | community.wave.seqera.io/library/samtools@sha256:5e5188b3dacf0117c4d4db891f2b9f4873d93288b4cfb26c1fba43eb7d241a4c |
| star | community.wave.seqera.io/library/star:2.7.11b--822039d47adf19a7 | community.wave.seqera.io/library/star@sha256:25a18a65561065970d2bee3cc481949b58711c31b48499dd4e86df18aa69e3a8 |
| featureCounts | community.wave.seqera.io/library/subread:2.1.1--0ac4d7e46cd0c5d7 | community.wave.seqera.io/library/subread@sha256:4b5569b45ab8d6f106b69cf5f682c34fce76e36871c8553f899eb54da58deb48 |
| gffread | quay.io/biocontainers/gffread:0.12.7--h077b44d_6 | quay.io/biocontainers/gffread@sha256:f603c5f4c8fff454ab282e917068d91060794e1daa798cbde0c125ade72f189d |
| picard | community.wave.seqera.io/library/picard | community.wave.seqera.io/library/picard@sha256:e269216786463d44f9d83a0d6e877b34bca2c7b4d35211b4b369fe98e39ef1a5|
| bedtools | community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4 | community.wave.seqera.io/library/bedtools@sha256:a49419dd283ad6a2e7bb00cab8f9d939b73581f1fe45fdbe774fd6cd9e652f48 |


A python and an R image were needed for visualization and differential gene analysis, after a few smaller images i concentrated all needed functions in the following two images. Below are the Dockerfiles for their creation.

python-image: python_vis
R-image: r_tricho

```bash
FROM python:3.10-slim
RUN pip install --no-cache-dir numpy==1.26.0
RUN pip install --no-cache-dir pandas==2.2.3
RUN pip install --no-cache-dir matplotlib==3.7.1
RUN pip install --no-cache-dir seaborn==0.13.2
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*
```

```bash
FROM r-base@sha256:089317f336a61255bb35f1efd799820cef37136d2acf2a76ba4abb74af51d4a3
RUN apt-get update && apt-get install -y libcurl4-openssl-dev libxml2-dev libssl-dev
RUN R -e "install.packages(c('BiocManager', 'tidyverse'))"
RUN R -e "BiocManager::install(c('edgeR', 'tximport', 'limma', 'DESeq2'))"
```

Command used to create these images (has to be run in the same directory as the Dockerfile is located in, with no other Dockerfile being in it):

```bash
docker build -t name_of_your_new_image .
```


### testing and development

In order to test the functionality of the pipeline a smaller dataset was chosen to ensure short runtimes and to keep the hardware requirements in check.

The used dataset is a dedicated training set and can be found [here](https://bioinfogp.cnb.csic.es/files/samples/rnaseq/).

### Directory Structure

The pipelines necessary directory structure is outlined here, relative to the main directory of the project (from here on called "rna/").
In "rna/" nextflow.config, params.yml and rna_pipe.nf are necessary, the pipeline is run in this directory so all traces and reports will end up here as well.
Nextflow will create a work directory inside rna/ in which past runs are cached, it is advised to specify a large volume in the config file as the work directory, since nextflow
will otherwise quickly fill the entire disk. Using the default is usually NOT a great idea.
Here in "rna/" the directories "scripts" and "modules" are located as well as a .nextflow directory (created by nextflow).
Further all experiment directories must be subdirectories of "rna/", these then must have the subdirectories "results" and "material".

## Pipeline

The following part describes all parts of the pipeline as well as their purpose.

### Config

The file **nextflow.config** is the configuration file for the pipeline, the configfile in this repo has all settings tuned for this pipeline.
Here it is especially important to set the work directory to some directory with a lot of available space (>100 GB).
Further memory and cpus allocations for each process can be made here (i only allocated memory and cpus to processes that are especially ressource intensive and that benefit greatly from multi-threading). Parameters can be declared here, but i used a dedicated parameter file for this.

### Parameters

Most parameters at the start of the main "rna_pipe.nf" script are nonsensical as they are just placeholder. The actual parameters can be found in the paramter files "params.yml", "params2.yml". This may seem odd, but the intention behind this is that a user only needs to change the parameters in a dedicated parameter file and maybe the config file depending on resource availability.


### input files

All input files must be inside a subdirectory of the experiments directory called "material", further a csv file of the following layout is required:

```csv
name,forward,reverse,condition
```
The "name" contains the name of the sample without any file-type.
Second and third columns contain the relative paths to each sample (usually like this: name_of_experiment/material/sample_name.fq.gz), paired-end reads ar used.
The last column specifies what condition each sample has (control, more / less oxygen,...)

### Workflow

The file **rna_pipe.nf** defines the worklfow of the pipeline, it is the centerpiece tying it all together.
The pipeline is run by executing this file, no other commands are necessary, **rna_pipe.nf** will call all other scripts and handle their in -and output.
All parameters are saved in the file **params.yml**, trace and report of the run are given a timestamp at the end.
Consequently this command runs the pipeline:

```bash
nextflow run rna_pipe.nf -params-file params.yml -with-trace -with-report
```

The workflow is roughly split in 2 parts depending on the format of the annotation file: gff and gff3 set annotation_type = "gff", while gtf sets annotation_type = "gtf". This leads to many "duplicated" process calls, but is is a lot easier this way instead of having multiple checks throughout the worklflow.


### Modules

All processes have their own module, **rna_pipe.nf** does not contain any processes.
Each process itself is just a wrapper for whatever happens in the script block of the process, it defines the environment (here the docker container),
all inputs and outputs and where and how to save the results.
Nextflow handles parallelization itself, processes are therefore not necessarily executed in the order in which they are written in the workflow.
When the necessary input for a process exists Nextflow may start this process.
All processes receive the name of the experiment as an input, therefore it will not be repeated in all the modules inputs. This experiment name is used to distinguish between experiments and dynamically create the necessary directory structures so that the pipeline can be used for different experiments without having to touch any code.
All process names are the same as the name of the modules they can be found in, for example the process "fastqc" is found in the module "fastqc.nf".

### fastqc

purpose: get an initial quality report of the reads
input: csv file containing name, forwards and reverse reads as well as the condition of each sample
output: report in html, and a zip-folder containing the reports

### trim_galore

purpose: trim the reads received after sequencing
input: forward and reverse reads
output: trimmed reads, trimming report, reports for each read in html and a corresponding zip-folder

### multiqc

purpose: bundle trimming reports and inital quality reports into one report
input: fastqc report, trimming report, individual reports
output: combined report in html (and additional data file)

### create_gtf

this process is only called when the annotation_type is gff
purpose: create agtf file from the gff or gff3 file
input: annotation file (gff)
output: gtf file

### create_ref_transcriptome

purpose: creating a reference transcriptome for salmon and kallisto
input: annotation file (either gtf or gff), reference genome in fasta format
output: reference transcriptome in fasta format

### make_tx2gene_from_gff

purpose: create a tx2gene file that maps gene IDs to transcript IDs --> important for differential gene expression
input: annotation file (gtf works as well), python script doing the mapping
output: csv file with 2 columns (gene id and transcript id)

### kallisto_index

purpose: creates an index for the quantification with kallisto, without an index the quantification does not work
input: reference transcriptome (fasta) --> NOT the reference genome
output: kallisto index file (.idx)

### kallisto

purpose: quantification with the kallisto pseudo-aligner
input: forward and reverse reads (pairwise), kallisto index
output: abundance file (tsv) for each pair of reads

### clean_abundance

purpose: clean the generated abundance files, remove unnecessary strings
input: abundance file (tsv)
output: cleaned abundance (tsv)
**clean_quant** works in the exact same way

### dge_edger_kallisto

purpose: calculate differential gene expression using EdgeR for Kallisto. For this the R script **edger_v2.r** is called.
input: all generated abundance.tsv files, edger_v2.r script, csv-file containing paths to fastq files (also contains naming patterns and condition), tx2gene file
output: csv files for all pairwise comparisons between abundance files.
note: this is process is mirrored by dge_edger_salmon, dge_edger_bowtie2 and dge_edger_star, they do exactly the same but became their own process in order to keep each process small so that every process only handles exactly one thing, nothing more.

### dge_limma_kallisto

purpose: calculate differential gene expression using Limma-voom for Kallisto. For this the R script **limma_v2.r** is called.
input: all generated abundance.tsv files, edger_v2.r script, csv-file containing paths to fastq files (also contains naming patterns and condition), tx2gene file
output: csv files for all pairwise comparisons between abundance files, a single voom-plot (png).

### dge_deseq2_kallisto

purpose: calculate differential gene expression using DESeq2 for Kallisto. For this the R script **deseq2_v2.r** is called.
input: all generated abundance.tsv files, edger_v2.r script, csv-file containing paths to fastq files (also contains naming patterns and condition), tx2gene file
output: csv files for all pairwise comparisons between abundance files

### plot_volcano_edger_kallisto

purpose: create a volcano and an MA plot to visualize the differential gene expression results.
input: csv file produced by differential gene expressiona analysis, **volcano_edger.r** R script
output: two png files, one volcano plot, one MA plot

### plot_volcano_limma_kallisto

purpose: create a volcano and an MA plot to visualize the differential gene expression results.
input: csv file produced by differential gene expressiona analysis, **volcano_limma.r** R script
output: two png files, one volcano plot, one MA plot

### plot_volcano_deseq2_kallisto

purpose: create a volcano and an MA plot to visualize the differential gene expression results.
input: csv file produced by differential gene expressiona analysis, **volcano_deseq2.r** R script
output: two png files, one volcano plot, one MA plot


### salmon_index

purpose: create the index files for salmon
input: reference transcriptome
output: directory containing all necessary index files

### salmon

purpose: quantification using salmon
input: paired end reads, index directory
output: quant.sf file for every pair of reads (renamed to reflect the inputs name)

### clean_quant

purpose: clean up unnecessary strings
input: quantifcation files (.sf)
output: cleaned quantification files

### dge_edger_salmon

purpose: calculate differential gene expression using EdgeR for Salmon. For this the R script **edger_v2.r** is called.
input: all generated abundance.tsv files, edger_v2.r script, csv-file containing paths to fastq files (also contains naming patterns and condition), tx2gene file
output: csv files for all pairwise comparisons between abundance files.

### dge_limma_salmon

purpose: calculate differential gene expression using Limma-voom for Salmon. For this the R script **limma_v2.r** is called.
input: all generated abundance.tsv files, edger_v2.r script, csv-file containing paths to fastq files (also contains naming patterns and condition), tx2gene file
output: csv files for all pairwise comparisons between abundance files, a single voom-plot (png).

### dge_deseq2_salmon

purpose: calculate differential gene expression using DESeq2 for Salmon. For this the R script **deseq2_v2.r** is called.
input: all generated abundance.tsv files, edger_v2.r script, csv-file containing paths to fastq files (also contains naming patterns and condition), tx2gene file
output: csv files for all pairwise comparisons between abundance files

### plot_volcano_edger_salmon

purpose: create a volcano and an MA plot to visualize the differential gene expression results.
input: csv file produced by differential gene expressiona analysis, **volcano_edger.r** R script
output: two png files, one volcano plot, one MA plot

### plot_volcano_limma_salmon

purpose: create a volcano and an MA plot to visualize the differential gene expression results.
input: csv file produced by differential gene expressiona analysis, **volcano_limma.r** R script
output: two png files, one volcano plot, one MA plot

### plot_volcano_deseq2_salmon

purpose: create a volcano and an MA plot to visualize the differential gene expression results.
input: csv file produced by differential gene expressiona analysis, **volcano_deseq2.r** R script
output: two png files, one volcano plot, one MA plot


### bowtie2_index

purpose: build the necessary index files for the alignment with Bowtie2
input: reference genome
output: six .bt2 files in their own directory

### bowtie2

purpose: alignment with Bowtie2
input: paired end reads, index files
output: one sam file for each pair of reads

### sam_to_bam

purpose: convert sam files to their binary form
input: .sam files
output: .bam files

### samtools_bowtie2

purpose: sort bam files by query and coordinates, index them and finally mark the duplicates
input: bam-files
output: bam-files that are sorted, indexed and have the duplicates marked

### picard_alignment_summary_bowtie2

purpose: get a report / summary of the alignment
input: bam-files, reference genome (fasta format)
output: alignment summary for each bam-file (txt)

### featurecounts_bowtie2

purpose: get some actual (human-readable) counts from the bam-files
input: bam-files, reference gtf (if annotation is gff --> created version)
output: counts in txt format, one for each bam-file

### calculate_tpms_bowtie2

purpose: calculate the tpm values from the counts found in the counts textfiles
input: counts files from featureCounts (txt), python script(do_calculate_tpms.py)
output: counts textfiles, with two additional columns: RPK and TPM

### clean_features_bowtie2

purpose: clean up unnecessary strings
input: counts-textfiles from calculate_tpms_bowtie2
output: cleaned counts-textfiles

### dge_edger_bowtie2

purpose: calculate differential gene expression using EdgeR for Bowtie2. For this the R script **edger_v2.r** is called.
input: all generated abundance.tsv files, edger_v2.r script, csv-file containing paths to fastq files (also contains naming patterns and condition), tx2gene file
output: csv files for all pairwise comparisons between abundance files.

### dge_limma_bowtie2

purpose: calculate differential gene expression using Limma-voom for Bowtie2. For this the R script **limma_v2.r** is called.
input: all generated abundance.tsv files, edger_v2.r script, csv-file containing paths to fastq files (also contains naming patterns and condition), tx2gene file
output: csv files for all pairwise comparisons between abundance files, a single voom-plot (png).

### dge_deseq2_bowtie2

purpose: calculate differential gene expression using DESeq2 for Bowtie2. For this the R script **deseq2_v2.r** is called.
input: all generated abundance.tsv files, edger_v2.r script, csv-file containing paths to fastq files (also contains naming patterns and condition), tx2gene file
output: csv files for all pairwise comparisons between abundance files

### plot_volcano_edger_bowtie2

purpose: create a volcano and an MA plot to visualize the differential gene expression results.
input: csv file produced by differential gene expressiona analysis, **volcano_edger.r** R script
output: two png files, one volcano plot, one MA plot

### plot_volcano_limma_bowtie2

purpose: create a volcano and an MA plot to visualize the differential gene expression results.
input: csv file produced by differential gene expressiona analysis, **volcano_limma.r** R script
output: two png files, one volcano plot, one MA plot

### plot_volcano_deseq2_bowtie2

purpose: create a volcano and an MA plot to visualize the differential gene expression results.
input: csv file produced by differential gene expressiona analysis, **volcano_deseq2.r** R script
output: two png files, one volcano plot, one MA plot


### star_index

purpose: create index files for star
input: reference genome (fasta), reference gtf
output: multiple files in their own directory, these files only need to be created once

### star

purpose: alignment with the STAR aligner
input: paired-end reads, created star index directory
output: all files have characteristic file endings, the prefix is determined by the inputs name. Log.out, Log.progress.out, Log.final.out are all lofgiles, Aligned.sortedByCoord.out.bam is the counts file, ReadsPerGene.out.tab is a table containing gene counts (not used here for the sake of a uniform analysis with bowtie2)

### samtools_star

purpose: sort bam files by query and coordinates, index them and finally mark the duplicates
input: bam-files
output: bam-files that are sorted, indexed and have the duplicates marked

### picard_alignment_summary_star

purpose: get a report / summary of the alignment
input: bam-files, reference genome (fasta format)
output: alignment summary for each bam-file (txt)

### featurecounts_star

purpose: get some actual (human-readable) counts from the bam-files
input: bam-files, reference gtf (if annotation is gff --> created version)
output: counts in txt format, one for each bam-file

### star_stats

purpose: get some stats about each .bam file that STAR produces, used to generate info about the mapping quality
input: one bam-file
output: text file containing all information obtained by samtools stats.

### star_calc_genome_cov

purpose: calculate a bedgraph from a bam-file in order to get information about the genome coverage
input: one bam-file, gtf-file corresponding to the organism
output: bedgraph for the bam-file, this serves as the input for the following visualization

### star_genome_coverage

purpose: plot genome coverage based on the bedgraph
input: bedgraph, "genome_coverage_vis.py" python script
output: png containing all plots

### get_scaffold_order

purpose: determine the order in which the scaffolds are found in a bam-file, additionally create a genome.txt file which is sorted in the order in which the scaffolds are
input: one bam-file
output: genome.txt, scaffold-order text file

### bedtools_sort_gtf

purpose: sort the gtf file, so that it matches the order in which the genome.txt and the bam file are
input: unsorted gtf file, genome.txt file (serves as the template)
output: sorted gtf file

### star_calc_gene_cov

purpose: calculate the coverage per gene for one bam-file
input: bam-file, sorted gtf-file, genome.txt
output: single text file containing all information produced by bedtools coverage

### star_gene_coverage

purpose: plot the findings of **star_calc_gene_cov**
input: text file with the data from the coverage per gene calculation, "gene_overage_vis.py" python script
output: single png file with multiple subplots

### star_mq

purpose: plot mapping quality as well as some general stats about a bam-file
input: stas text file generated by **star_stats**, python script "mq_vis.py"
output: one png file showcasing the findings of samtools stats

### calculate_tpms_star

purpose: calculate the tpm values from the counts found in the counts textfiles
input: counts files from featureCounts (txt), python script(do_calculate_tpms.py)
output: counts textfiles, with two additional columns: RPK and TPM

### clean_features_star

purpose: clean up unnecessary strings
input: counts-textfiles from calculate_tpms_star
output: cleaned counts-textfiles

### dge_edger_bowtie2

purpose: calculate differential gene expression using EdgeR for STAR. For this the R script **edger_v2.r** is called.
input: all generated abundance.tsv files, edger_v2.r script, csv-file containing paths to fastq files (also contains naming patterns and condition), tx2gene file
output: csv files for all pairwise comparisons between abundance files.

### dge_limma_bowtie2

purpose: calculate differential gene expression using Limma-voom for STAR. For this the R script **limma_v2.r** is called.
input: all generated abundance.tsv files, edger_v2.r script, csv-file containing paths to fastq files (also contains naming patterns and condition), tx2gene file
output: csv files for all pairwise comparisons between abundance files, a single voom-plot (png).

### dge_deseq2_bowtie2

purpose: calculate differential gene expression using DESeq2 for STAR. For this the R script **deseq2_v2.r** is called.
input: all generated abundance.tsv files, edger_v2.r script, csv-file containing paths to fastq files (also contains naming patterns and condition), tx2gene file
output: csv files for all pairwise comparisons between abundance files

### plot_volcano_edger_star

purpose: create a volcano and an MA plot to visualize the differential gene expression results.
input: csv file produced by differential gene expressiona analysis, **volcano_edger.r** R script
output: two png files, one volcano plot, one MA plot

### plot_volcano_limma_star

purpose: create a volcano and an MA plot to visualize the differential gene expression results.
input: csv file produced by differential gene expressiona analysis, **volcano_limma.r** R script
output: two png files, one volcano plot, one MA plot

### plot_volcano_deseq2_star

purpose: create a volcano and an MA plot to visualize the differential gene expression results.
input: csv file produced by differential gene expressiona analysis, **volcano_deseq2.r** R script
output: two png files, one volcano plot, one MA plot



### workflow.onComplete

This part of the workflow is executed whenn everything else is done. Here the trace and report are given a timestamp so that multiple runs with different parameters can be easily compared


### Transfer of experimental data

Since the small disk space of the VM is insufficient a mounted volume was used to store experimental data (as well as nextflows work directory).
This volume has different access permissions than the VM, so in order to make rsync work I forced it to give me full permissions:

```bash
sudo chown -R ubuntu:ubuntu /home/ubuntu/point/tricho_unperturbed/material
```

In order to upload the data from the SSD to the server this command was used (exact path depending on current experiment)

```bash
rsync -avu --progress -e "ssh -i ~/Documents/uni/semester8/praxis/privKey.org_ecdsa -p 30167" /media/simon/T9/NGS/Metagenomics/Biofilm_202505 ubuntu@129.70.51.6:point/tricho_unperturbed/material/
```



## tmux sessions for longer runs

```bash
# start session
tmux new -s mysession

# detach
ctrl + b, after that press d

# reconnect to the session
tmux attach -t mysession
```

### transfer results to local machine

```bash
rsync -avu --progress -e "ssh -i ~/Documents/uni/semester8/praxis/privKey.org_ecdsa -p 30167" ubuntu@129.70.51.6:project/rna/rna_report.html ~/Documents/uni/semester8/praxis/rna/rna_report.html
```

### creation of flowcharts

To visualize the workflow of the pipeline Nextflows own visualization was used.

```bash
nextflow run rna_pipe.nf -preview -with-dag rna_mermaid.mmd
```

This does NOT produce a png file, but rather a mermaid file which was then converted (and customized) using the [Mermaid Live Editor](mermaid.live).