library(tximport)
library(edgeR)
library(dplyr)

options(error=function() {traceback(2); quit("no", 1, FALSE)})

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

# this must be adjusted so all types of files can be handled properly -> different alinger = different names
# read all files --> patterns must be adjusted according to whatever the naming convention for the condition is
#if(import_type == "kallisto"){
#    #final_path = path_to_all + aligner + "/"
#    # since kallisto always writes files called "abundance.tsv" every replicate has its own folder --> the provided path (from CL) is then just results/kallisto/
#    control <- list.files(path=path_to_all, pattern="*ph5_co004*abundance.tsv", full.names = TRUE)
#    condition_1 <- list.files(path=path_to_all, pattern="*ph5_co50*abundance.tsv", full.names = TRUE)
#    condition_2 <- list.files(path=path_to_all, pattern="*ph3_co004*abundance.tsv", full.names = TRUE)

#} else if (import_type == "salmon") {
#    # same story as with kallisto, only that its "quant.sf" instead of "abundance.tsv"
#    control <- list.files(path=path_to_all, pattern="*ph5_co004*/quant.sf", full.names = TRUE)
#    condition_1 <- list.files(path=path_to_all, pattern="*ph5_co50*/quant.sf", full.names = TRUE)
#    condition_2 <- list.files(path=path_to_all, pattern="*ph3_co004*/quant.sf", full.names = TRUE)

#}else if (import_type == "none" && aligner == "bowtie2") {
#    # here the path to control, condition1, condition2 are exactly the same --> its all just in results/bowtie2
#    control <- list.files(path=path_to_all, pattern="*ph5_co004*counts.txt", full.names = TRUE)
#    condition_1 <- list.files(path=path_to_all, pattern="*ph5_co50*counts.txt", full.names = TRUE)
#    condition_2 <- list.files(path=path_to_all, pattern="*ph3_co004*counts.txt", full.names = TRUE)

#}else if (import_type == "none" && aligner == "star") {
    # the same as with bowtie2
#    control <- list.files(path=path_to_all, pattern="*ph5_co004*counts.txt", full.names = TRUE)
#    condition_1 <- list.files(path=path_to_all, pattern="*ph5_co50*counts.txt", full.names = TRUE)
#    condition_2 <- list.files(path=path_to_all, pattern="*ph3_co004*counts.txt", full.names = TRUE)
#}


# assuming three groups of course...which means 9 samples overall --> future fix: define samples per group and number of groups in CL
n <- samples_amount
#control <- file_paths[1:n]
#condition1 <- file_path[(n+1):(2*n)]
#condition2 <- file_path[(2*n+1):(3*n)]

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
}else {
    gene_counts_control <- tximport(control, type = import_type, txIn = FALSE, txOut = FALSE, geneIdCol = "Geneid", countsCol = "Counts", lengthCol = "Length", abundanceCol = "TPM")
    gene_counts_condition_1 <- tximport(condition1, type = import_type, txIn = FALSE, txOut = FALSE, geneIdCol = "Geneid", countsCol = "Counts", lengthCol = "Length", abundanceCol = "TPM")
    gene_counts_condition_2 <- tximport(condition2, type = import_type, txIn = FALSE, txOut = FALSE, geneIdCol = "Geneid", countsCol = "Counts", lengthCol = "Length", abundanceCol = "TPM")
}

# all files per condition are combined into one list
gene_counts_combined <- list(control=gene_counts_control, condition1 = gene_counts_condition_1, condition2 = gene_counts_condition_2)

# after combining into the list, now the actual counts are merged
counts <- cbind(gene_counts_control$counts, gene_counts_condition_1$counts, gene_counts_condition_2$counts)
# quick check
head(counts)

# ------------------------------------------------------------------------------


# analysis of differential gene expression


# depending on number of groups: define them
experimental_groups <- factor(c(rep("control", 3), rep("condition1", 3), rep("condition2", 3)))# 3 = replicates per group

# if anything was messed up: correct it here
experimental_groups <- factor(experimental_groups, levels = c("control", "condition1", "condition2"))


# why on earth is always y used...
# here the DGEList object is defined
y <- DGEList(counts = counts, group = experimental_groups)


# filtering for low expression --> only want reads with statistical power
good <- filterByExpr(y)
y <- y[good, ]#, good.lib.sizes=FALSE

# some normalization
y <- normLibSizes(y)




# design matrix, linear model is fitted --> thats what glmQLFit does
design_matrix <- model.matrix(~0+experimental_groups)
fit <- glmQLFit(y, design_matrix)

# now the actual comparisons
# with testing --> glmGLFTest

# control vs condition1
test_control_vs_condition1 <- glmQLFTest(fit, contrast = c(-1, 1, 0))
res_control_vs_condition1 <- topTags(test_control_vs_condition1, n = Inf) # toptags are there to get the results
# now write it all to csv --> easier for further analysis --> can be easily passed to other modules
write.csv(res_control_vs_condition1, "EdgeR_control_condition1.csv")

# control vs condition2
test_control_vs_condition2 <- glmQLFTest(fit, contrast = c(-1, 1, 0))
res_control_vs_condition2 <- topTags(test_control_vs_condition2, n = Inf)
write.csv(res_control_vs_condition2, "EdgeR_control_condition2.csv")