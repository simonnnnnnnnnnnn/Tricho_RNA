library(tximport)
library(edgeR)
library(dplyr)

options(error=function() {traceback(2); quit("no", 1, FALSE)})


args <- commandArgs(trailingOnly = TRUE)
path_to_mapping <- args[1]
#samples_amount <- as.numeric(args[2])
aligner <- args[2]
overview_file <- args[3]

# now all the files coming in
file_path <- args[4:(length(args))]

# define type according to aligner
import_type <- switch(aligner, "kallisto"="kallisto", "salmon"="salmon", "bowtie2"="none", "star"="none", "unknown")

if (import_type == "unknown"){
    stop("aligner type invalid, arguments probably in wrong order")
}

# read the overview file
infile <- read.csv(overview_file, sep = ",", header = TRUE)

# get distinct conditions --> this will define the new df
distinct_conditions <- unique(infile$condition)

# new list to match incoming files to their respective conditions --> all dynamically
temp <- setNames(vector("list", length(distinct_conditions)), distinct_conditions)

for (i in 1:nrow(infile)){
    h <- grep(infile$name[i], file_path, value = TRUE)

    place <- infile$condition[i]
    temp[[place]] <- append(temp[[place]], h)
}

# now they must go through tximport
if (import_type %in% c("kallisto", "salmon")){
    gene_to_transcript <- read.csv(path_to_mapping, header = TRUE)

    gene_counts_combined <- list()
    for(i in 1:length(distinct_conditions)){
        #gene_counts_combined <- append(gene_counts_combined, tximport(temp$distinct_conditions[i], type = import_type, tx2gene = gene_to_transcript, ignoreTxVersion = TRUE))

        # less cluttered --> that is a lot more readable, should do the same :D
        condition_files <- temp[[distinct_conditions[i]]]
        txi <- tximport(condition_files, type = import_type, tx2gene = gene_to_transcript, ignoreTxVersion = TRUE)
        gene_counts_combined[[distinct_conditions[i]]] <- txi
    }

    # ensure order --> bind all columns, but preserve the names
    counts <- do.call(cbind, lapply(gene_counts_combined, function(x) x[["counts"]]))
    # deseq2 doesnt like floats...
    counts <- round(counts)

}else{
    # here bowtie2 and star
    gene_counts_combined <- list()
    for (i in 1:length(distinct_conditions)){
        condition_files <- temp[[distinct_conditions[i]]]
        txi <- tximport(condition_files, type = import_type, txIn = FALSE, txOut = FALSE, geneIdCol = "Geneid", countsCol = "Counts", lengthCol = "Length", abundanceCol = "TPM")
        gene_counts_combined[[distinct_conditions[i]]] <- txi
    }
    # ensure order --> bind all columns, but preserve the names
    counts <- do.call(cbind, lapply(gene_counts_combined, function(x) x[["counts"]]))
}

# group samples to conditions
samples_per_condition <- sapply(temp, length)
experimental_groups <- factor(unlist(mapply(rep, distinct_conditions, samples_per_condition)))

# if anything was messed up: correct it here
experimental_groups <- factor(experimental_groups, levels = distinct_conditions)# use ditinct_conditions instead


# dge with edger
y <- DGEList(counts = counts, group = experimental_groups)

# filtering for low expression --> only want reads with statistical power
good <- filterByExpr(y)
y <- y[good, ]#, good.lib.sizes=FALSE

# some normalization
y <- normLibSizes(y)


# design matrix, linear model is fitted --> thats what glmQLFit does
design_matrix <- model.matrix(~0+experimental_groups)
fit <- glmQLFit(y, design_matrix)

# results and tests
pairwise_results <- list()

for (i in 1:(length(distinct_conditions)-1)){
    for (j in (i+1):length(distinct_conditions)){
        num <- distinct_conditions[i]
        denom <- distinct_conditions[j]

        contrast_dynamic <- rep(0, length(distinct_conditions))
        contrast_dynamic[i] <- -1
        contrast_dynamic[j] <- 1

        test <- glmQLFTest(fit, contrast = contrast_dynamic)
        h <- topTags(test, n = Inf)
        name <- paste0("EdgeR_", num, "_", denom)
        pairwise_results[[name]] <- h
    }
}

# save it --> generate files nextflow expects
for (comparison in names(pairwise_results)){
    tempres <- as.data.frame(pairwise_results[[comparison]])
    write.csv(tempres, paste0(comparison, ".csv"), row.names = TRUE)# here row.names must be true --> keep IDs
}