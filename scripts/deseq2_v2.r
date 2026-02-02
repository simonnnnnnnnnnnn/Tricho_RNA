library(tximport)
library(dplyr)
library(DESeq2)

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

# new df with colnames being named after the conditions --> dynamically decide number of conditions --> makes the script a lot more flexible
#temp <- data.frame(matrix(ncol=length(distinct_conditions), nrow = 0))
#colnames(temp) <- distinct_conditions

# iterate over infile: use name pattern from infile to get all files from file_path that match --> the append them to the corresponding conditions in the them df
# this way all files coming from featureCounts or pseudo alignment are sorted by their experimental condition

#for (i in 1:nrow(infile)){
#    h <- grep(infile$name[i], file_path, value = TRUE)
#    place <- infile$condition[i]
#    temp[[place]] <- append(temp[[place]], h)
#}

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


# now the deseq part --> conditionsnamen dynamisch # samples_maount benutzen
#experimental_groups <- factor(c(rep("control", 3), rep("condition1", 3), rep("condition2", 3)))# 3 = replicates per group --> use samples_amount
samples_per_condition <- sapply(temp, length)
experimental_groups <- factor(unlist(mapply(rep, distinct_conditions, samples_per_condition)))

# if anything was messed up: correct it here
experimental_groups <- factor(experimental_groups, levels = distinct_conditions)# use ditinct_conditions instead


# for DESeq2 the data needs to be formatted a bit differently:
col_data <- data.frame(row.names = colnames(counts), condition = experimental_groups)

# now that a suitable DF exists the DESeq dataset can be created
dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~ condition)

# run it
dds_run <- DESeq(dds)

# get some results
#res_control_vs_condition1 <- results(dds_run, contrast = distinct_conditions)# same story here
#write.csv(res_control_vs_condition1, "DESeq2_control_condition1_v2.csv")

#res_control_vs_condition2 <- results(dds_run, contrast = c("condition", "control", "condition2"))# c("control", "condition1", "condition2")
#write.csv(res_control_vs_condition2, "DESeq2_control_condition2_v2.csv")


# all pairwise results
pairwise_results <- list()

for (i in 1:(length(distinct_conditions)-1)){
    for (j in (i+1):length(distinct_conditions)){
        num <- distinct_conditions[i]
        denom <- distinct_conditions[j]

        h <- results(dds_run, contrast = c("condition", num, denom))

        name <- paste0("DESeq2_", num, "_", denom)
        pairwise_results[[name]] <- h # put the results in  the list --> later iterate over it and save all files
    }
}

# save it all
for (comparison in names(pairwise_results)){
    tempres <- as.data.frame(pairwise_results[[comparison]])
    write.csv(tempres, paste0(comparison, ".csv"), row.names = TRUE)# here row.names must be true --> keep IDs
}