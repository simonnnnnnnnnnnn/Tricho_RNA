library(DESeq2)
library(tximport)
library(dplyr)



# most of it is very similar to EdgeR --> at least the data input and handling


# get usefule error messages, --> more details so one actually sees what broke
options(error=function() {traceback(2); quit("now", 1, FALSE)})

# read CL arguments --> file paths
args <- commandArgs(trailingOnly = TRUE)
path_to_mapping <- args[1]
samples_amount <- as.numeric(args[2])
aligner <- args[3]

# depending on which aligner is used the tximport method is chosen
# STARs output can technically be parsed with import type "rsem", but since STARs output goes into featureCounts it must be set to "none"
import_type <- switch(aligner, "kallisto"="kallisto", "salmon"="salmon", "bowtie2"="none", "star"="none", "unknown") # i hope unknown never triggers

# just to make sure...
if (import_type == "unknown"){
    stop("aligner type not supported, how did you even do that?")
}
# data preparation
#final_path <- paste(path_to_all, aligner, "/")
file_path <- args[3:(length(args))]



# sadly they dont come in in the expected order.... fix for later: pass the pattern via CL as well
control <- grep("ph5_co004", file_path, value = TRUE)
condition1 <- grep("ph5_co50", file_path, value = TRUE)
condition2 <- grep("ph3_co004", file_path, value = TRUE)

# file mapping transcripts to genes
# not needed for Bowtie and STAR --> they are already at gene level
if(import_type %in% c("kallisto", "salmon")){
    gene_to_transcript <- read.csv(path_to_mapping, header = TRUE)
}
#gene_to_transcript <- read.table(path_to_mapping, header = FALSE)


# now all the listed files must be imported, the counts of the transcripts are converted to gene counts
# for this the gene_to_transcript file is used
if(import_type %in% c("kallisto", "salmon")){
    gene_counts_control <- tximport(control, type=import_type, tx2gene = gene_to_transcript, ignoreTxVersion = TRUE)#ignoreTxVersion = TRUE
    gene_counts_condition_1 <- tximport(condition1, type=import_type, tx2gene = gene_to_transcript, ignoreTxVersion = TRUE)
    gene_counts_condition_2 <- tximport(condition2, type=import_type, tx2gene = gene_to_transcript, ignoreTxVersion = TRUE)

    # deseq only likes integers so some rounding must take place
    gene_counts_combined <- list(control=gene_counts_control, condition1 = gene_counts_condition_1, condition2 = gene_counts_condition_2)

    counts <- cbind(gene_counts_control$counts, gene_counts_condition_1$counts, gene_counts_condition_2$counts)
    counts <- round(counts) # this is the crucial step for deseq

}else {
    gene_counts_control <- tximport(control, type = import_type, txIn = FALSE, txOut = FALSE, geneIdCol = "Geneid", countsCol = "Counts", lengthCol = "Length", abundanceCol = "TPM")
    gene_counts_condition_1 <- tximport(condition1, type = import_type, txIn = FALSE, txOut = FALSE, geneIdCol = "Geneid", countsCol = "Counts", lengthCol = "Length", abundanceCol = "TPM")
    gene_counts_condition_2 <- tximport(condition2, type = import_type, txIn = FALSE, txOut = FALSE, geneIdCol = "Geneid", countsCol = "Counts", lengthCol = "Length", abundanceCol = "TPM")

    gene_counts_combined <- list(control=gene_counts_control, condition1 = gene_counts_condition_1, condition2 = gene_counts_condition_2)

    counts <- cbind(gene_counts_control$counts, gene_counts_condition_1$counts, gene_counts_condition_2$counts)

}

# all files per condition are combined into one list
#gene_counts_combined <- list(control=gene_counts_control, condition1 = gene_counts_condition_1, condition2 = gene_counts_condition_2)

# after combining into the list, now the actual counts are merged
#counts <- cbind(gene_counts_control$counts, gene_counts_condition_1$counts, gene_counts_condition_2$counts)
# quick check
#head(counts)


# Now the divergence from the EdgeR script

experimental_groups <- factor(c(rep("control", 3), rep("condition1", 3), rep("condition2", 3)))# 3 = replicates per group

# if anything was messed up: correct it here
experimental_groups <- factor(experimental_groups, levels = c("control", "condition1", "condition2"))


# for DESeq2 the data needs to be formatted a bit differently:
col_data <- data.frame(row.names = colnames(counts), condition = experimental_groups)

# now that a suitable DF exists the DESeq dataset can be created
dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~ condition)

# run it
dds_run <- DESeq(dds)

# get some results
res_control_vs_condition1 <- results(dds_run, contrast = c("condition", "control", "condition1"))
write.csv(res_control_vs_condition1, "DESeq2_control_condition1.csv")

res_control_vs_condition2 <- results(dds_run, contrast = c("condition", "control", "condition2"))
write.csv(res_control_vs_condition2, "DESeq2_control_condition2.csv")
