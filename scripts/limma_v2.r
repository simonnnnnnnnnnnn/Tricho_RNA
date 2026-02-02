# libs
library(tximport)
library(limma)
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

# match infiles with patterns to conditions
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

# group by condition
samples_per_condition <- sapply(temp, length)
experimental_groups <- factor(unlist(mapply(rep, distinct_conditions, samples_per_condition)))

# if anything was messed up: correct it here
experimental_groups <- factor(experimental_groups, levels = distinct_conditions)# use ditinct_conditions instead



# now the limma
# a dgelist must be created
dge <- DGEList(counts = counts, group = experimental_groups)

design <- model.matrix(~0+experimental_groups)

# low expression not wanted --> kick out
good <- filterByExpr(dge, design)
dge <- dge[good, ]

dge <- calcNormFactors(dge) # this is the normalization

# V O O M
png("voomplot.png")
my_voom <- voom(dge, design, plot = TRUE)
dev.off()

limma_fit <- lmFit(my_voom, design)
limma_fit <- eBayes(limma_fit)

# now make comparisons for all pairs
pairwise_results <- list()

pairs <- combn(distinct_conditions, 2, simplify = FALSE)

design_colnames <- colnames(design)

for (pair in pairs){
    name <- paste0("Limma-voom_", pair[1], "_", pair[2])
    con1 <- design_colnames[grep(pair[1], design_colnames)]
    con2 <- design_colnames[grep(pair[2], design_colnames)]
    contrast_form <- paste0(con1, "-", con2)


    contrast_matrix <- makeContrasts(contrasts = contrast_form, levels = design)

    # now fit with contrast
    fit2 <- contrasts.fit(limma_fit, contrast_matrix)
    fit2 <- eBayes(fit2)

    output <- topTable(fit2, coef = 1, number = Inf)
    write.csv(output, paste0(name, ".csv"), row.names = TRUE)
}