library(tximport)
library(ggplot2)

options(error=function() {traceback(2); quit("no", 1, FALSE)})


args <- commandArgs(trailingOnly = TRUE)
samplename <- args[1]
#samples_amount <- as.numeric(args[2])
data <- args[2]

finalname <- paste0(samplename, "_volcano.png")
titlename <- paste0(samplename, "volcano plot")
data1 <- read.csv(data, header = TRUE, stringsAsFactors = FALSE)
data1 <- data1[!is.na(data1$padj) & !is.na(data1$baseMean) & data1$baseMean > 0, ]

print(head(data1))

#print(head(data$logFC))


# classify genes into significantly differentially expressed --> FDR < 0.05 and non-significant genes
# colors should mathc that

data1$differential <- ifelse(data1$padj < 0.05, "significant", "non-significant")
# add adjusted pvalues
data1$log10padj <- -log10(data1$padj)

print(head(data1))


volcano1 <- ggplot(data1, aes(x = log2FoldChange, y = log10padj, color = differential)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("significant" = "red", "non-significant" = "grey")) +
    theme_minimal() +
    labs(title = titlename, x = "logFC", y = "-log10(adj P-value)", color = "significance")

ggsave(finalname, plot = volcano1, width = 8, height = 6, dpi = 300, bg = "white")
#print(volcano1)

data1$log10baseMean <- log10(data1$baseMean +1)
# now the MA plot part
ma1 <- ggplot(data1, aes(x = log10baseMean, y = log2FoldChange, color = differential)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("significant" = "red", "non-significant" = "grey")) +
    theme_minimal() +
    labs(title = paste0(samplename, "_MA_Plot"),
    x = "log10 (basemean +1)",
    y = "log Fold Change",
    color = "Significance")
ggsave(paste0(samplename, "_MA_Plot.png"), plot = ma1, width = 8, height = 6, dpi = 300, bg = "white")