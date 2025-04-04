library(DESeq2)
library(EnhancedVolcano)

#example of 28 vs 20 celcius

dat <- read.csv("raw.csv", header = T, row.names = 1) #put only 2 conditions
info <- read.table("colData.txt", header = T, sep = '\t')

dds <- DESeqDataSetFromMatrix(dat, info, ~condition)

keep <- rowSums(counts(dds)) >= 10 # low-read genes filtering
dds <- dds[keep,]

ddsDE <- DESeq(dds)

res <- results(ddsDE, alpha = 0.05)
res.df <- as.data.frame(res)
res.dfO <- res.df[order(res.df$padj),]
write.csv(res.dfO, "DESeq.csv")

# Create a lab vector where genes starting with "LOC" have NA as their label
lab_vector <- row.names(res.df)
lab_vector[grepl("^LOC", row.names(res.df))] <- NA  # Set labels to NA for genes starting with "LOC"

b<-EnhancedVolcano(res.df,
                x = "log2FoldChange",
                y = "padj",
                lab = lab_vector,
                pCutoff = 0.00001,
                FCcutoff = 1.5,
                legendPosition = "none",
                title = NULL,
                ylab = "-log10(padj)",
                xlab = "log2(foldChange)",
                subtitle = NULL,
                caption = NULL,
                col = c('gray', 'gray', 'gray', 'blue'),
                drawConnectors = TRUE,
                widthConnectors = 0.1,
                cutoffLineType = "blank",
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                max.overlaps = 35)
#adjust color according to conditions

b <- b + ggtitle("28°C vs 20°C") + theme(plot.title = element_text(hjust = 0.5))
ggsave("axis.png", plot = b, width = 20, height = 13, dpi = 300, units = "cm", limitsize = FALSE)
