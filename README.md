# Tricho_RNA


**Objective:** Creation of an RNAseq pipeline with Nextflow for _Trichoplax adhaerens_

## Setup of Execution Environment

The Pipelines were executed on a virtual machine running Ubuntu 22.04, the VM itself ran on SimpleVMs Servers.

### Requirements

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
```

The commands to install nextflow follow [Nextflows Instructions](https://training.nextflow.io/latest/envsetup/02_local/#download-nextflow) for the local environment (2025-05-07).

### Transfer of Files to the VM

To keep a local directory (where the code is written) in sync with a directory on the VM **rsync** is used. The following example shows how this was done:

```b̀ash
rsync -avu --progress -e "ssh -i ~/path/to/my/private/key.org_ecdsa -p 12345" ~/path/tomy/local/directory/ ubuntu@123.45.67.8:path/to/vm/directory/
```

## Preparation Pipeline

Before running the pipeline some preparations have to take place.

### Nextflow

Nextflow itself requires the following the following software to run properly as is stated in [Nextflows Instructions](https://training.nextflow.io/latest/envsetup/02_local/#download-nextflow) (can be seen in the **Customisation of the VM** as well).

 - Docker
 - git
 - Java (versions 11 to 21)

### Docker

Docker of course requires images to run in its containers, the following table shows the program ran in a container, its cooresponding image and the source for the image.
All of these images can be installed using this command:

```b̀ash
docker pull name_docker_image
```

| program | image | source |
| ----- | ----- | ---- |
| fastqc & trim_galore | community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18 | d |
| multiqc | community.wave.seqera.io/library/pip_multiqc:ad8f247edb55897c | d |
| kallisto | quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1 | d |
| salmon | community.wave.seqera.io/library/salmon:1.10.3--fcd0755dd8abb423 | d |
| picard | community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6 | d |

### Directory Structure

In order to run the code without having to make any changes the directory structure should look like this:

_hier mit tree command graphik einfügen_

## Pipeline

The following part describes all parts of the pipeline as well as their purpose.

### Workflow

The file **rna_pipe.nf** defines the worklfow of the pipeline, it is the centerpiece tying it all together.
The pipeline is run by executing this file, no other commands are necessary, **rna_pipe.nf** will call all other scripts and handle their in -and output.
Consequently this command runs the pipeline:

```b̀ash
nextflow run rna_pipe.nf
```

### Config

The file **nextflow.config** is the configuration file for the pipeline, the configfile in this repo has all settings tuned for this pipeline.

### Parameters

There were no parameters given  when running the Pipeline, this is because all parameters come with default values. Furthermore the file **params.yml** specify all used parameters,
values given in params.yml overwrite the defaults set in the code, therefore the pipelines parameters can be altered without touching the code itself.

### Modules

The modules contain the processes, here the actual work is done. All used tools have at least one Module dedicated to them. Each module has one process which runs a docker container with the necessary tool for the step. Modules and processes share the same name which usually resembles the name of the used tool or purpose of the tool. The results are always saved in directories named after the used tools. The file docker_images.csv containes the docker images of the used tools.

#### fastqc

This module contains a single process that creates a docker container and runs [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) inside it. The purpose of this is to create an inital quality reports of the reads.

**Input:** single gzipped fastq file (.fq.gz)

**output:** zip folder, html-file

#### trim_galore

This module runs [TrimGalore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) and consequently trims the reads and produces trimmed reads, trimming reports of the paired-end reads as well as a general trimming report.

**input:** paired-end reads (.fq.gz), one forward and one reverse read

**output:** trimmed reads (fq.gz), trimming_report (.txt), fastqc_reports for each read (.zip, .html)

#### multiqc

This module creates a comprehensive quality report by combining all other quality reports into one single report.

**input:** all quality reports, parameter with the intended name for the report

**output:** comprehensive report (html)

#### kallisto_index

kallisto_index creates the index for kallisto, this functionality is separated in its own modules because the creation of the index only needs to happen once.

**input:** reference genome (.fasta / .fa)

**output:** kallisto index (.idx) 

#### Kallisto

[Kallisto](https://github.com/pachterlab/kallisto) is the first pseudo-aligner in this pipeline and is used to primarily get TPM values for the transcripts.

**input:** paired-end reads (fq.gz), one forward and one reverse, index for kallisto (.idx)

**output:** abundance file containing the TPM values (.tsv)

#### salmon_index

This module creates the salmon index for the following salmon module, unlike kallisto_index multiple files are generated all of which are needed for the following quantification

**input:** reference genome in fasta format (.fasta / .fa)

**output:** directory containing multiple files salmon uses as the index

#### salmon

[Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) is the second pseudo-aligner used here, as with kallisto it is used to generate TPM values for the transcripts

**input:** paired-end reads (.fq.gz), path to the directory containing the index files

**output:** single quant.sf file

#### bowtie2_index

Bowtie2, the first true aligner in this pipeline, needs an index for the alignment. This module created this index.

**input:** reference genome is fasta format (.fasta / .fa)

**output:** directory containing six .bt2 files, all of which Bowtie2 needs for the following quantification.

#### bowtie2

This module does the alignment using [Bowtie2](https://github.com/BenLangmead/bowtie2).

**input:** paired-end reads (.fq.gz), path to the directory containing the .bt2 files

**output:** single .sam file (one .sam file from 2 reads (forward and reverse))

#### samtools

Samtools is used to convert .sam to .bam and to sort the bam-files by coordinates

**input:** single .sam file

**output:** single .bam file

#### stringtie

Stringtie is used to get the expression levels from the .bam files --> TPM values

**input:** .bam file, gene annotaion file (.gtf)

**output:** tsv file containing the expression levels, another gtf file (that is not the main objective here)

#### star_index

As with the other (pseudo) aligners STAR also needs an index for the alignment to take place

**input:** reference genome in fasta format (.fasta / .fa), reference gene annotation file (.gtf)

**output:** directory containing all the index files STAR generates

#### star

This module runs STAR and performs the quantification, all resulting files have predetermined names but can be given prefixes, only the predetermined names are listed here.

**input:** paired-end reads (.fq.gz), apth to the index-directory

**output:** Log.out, Log.progress.out, Log.final.out, Aligned.sortedByCoord.out.bam, ReadsPerGene.out.tab

#### picard_alignment_summary

