library(DESeq2)
library(ggplot2)

dat <- read.csv("raw.csv", header = T, row.names = 1) # after attached Entrez gene id using DAVID (filtered unmatched symbols)
info <- read.table("colData.txt", header = T, sep = '\t') # metaData

dds <- DESeqDataSetFromMatrix(dat, info, ~condition)

keep <- rowSums(counts(dds)) >= 10 # Low read filtering
dds <- dds[keep,]

ddsDE <- DESeq(dds)
vsd <- vst(ddsDE, blind=FALSE)

a <- plotPCA(vsd, intgroup=c("condition"))

b <- a + 
  theme_bw() +  # White background with no grid lines
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Set margins (top, right, bottom, left) to 0
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    text = element_text(size = 20, family = "Arial"),  # Set font to Arial
    axis.title = element_text(size = 19),  # Increase axis title size
    axis.text = element_text(size = 15),   # Increase axis labels size
    legend.text = element_text(size = 20),   # Increase legend text size
    legend.title = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.1),  # Outline the plot with black lines
  ) +
  geom_point(aes(color = condition, fill = condition), size = 8, stroke = 0.5, shape=21, color = "black") +
  scale_color_manual(values = c("20°C" = "blue", "28°C" = "orange", "37°C" = "red")) +
  scale_fill_manual(values = c("20°C" = "blue", "28°C" = "orange", "37°C" = "red")) +
  guides(
    fill = "none",  # Remove the fill legend (white-filled dots)
    color = guide_legend(override.aes = list(shape = 21, size = 8, stroke = 0.5, color = "black",fill = c("blue", "orange", "red")))
  )

print(b)

ggsave("PCAdeseq2013.png", plot = b, width = 20, height = 13, dpi = 300, units = "cm", limitsize = FALSE)

normCounts <- counts(ddsDE, normalized = T)
normCounts.df <- as.data.frame(normCounts)
write.csv(normCounts.df, "NormCounts.csv")
