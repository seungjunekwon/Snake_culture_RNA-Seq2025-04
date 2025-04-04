library(clusterProfiler)
library(org.Pguttatus.eg.db) # self-established (check 'establish-corn-snake-go.R')
library(ggplot2)
library(stringr)

# example of 20 vs 28

all <- as.character(unlist(read.table("28-20-all.txt"))) #txt files containing entrez ids
deg <- as.character(unlist(read.table("20-28-005-1.txt"))) #DEG selected with padj & foldchange threshold 0.05 & +-1

res <- enrichGO(deg, 'org.Pguttatus.eg.db', ont="BP",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                pAdjust = "BH",
                readable = TRUE,
                universe = all,
                keyType = "GID")

res_df <- as.data.frame(res)
res_df$GeneRatio <- as.numeric(str_split_fixed(res_df$GeneRatio, "/", 2)[,1]) / 
  as.numeric(str_split_fixed(res_df$GeneRatio, "/", 2)[,2])
res_df$BgRatio <- as.numeric(str_split_fixed(res_df$BgRatio, "/", 2)[,1]) / 
  as.numeric(str_split_fixed(res_df$BgRatio, "/", 2)[,2])

res_df$foldEnrichment <- res_df$GeneRatio / res_df$BgRatio

res_df <- res_df[, -1]
res_df <- res_df[, c(setdiff(names(res_df), "geneID"), "geneID")]

# Capitalize only the first letter of the first word in the 'Description' column
res_df$Description <- sapply(res_df$Description, function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
})

res_df$Description <- str_wrap(res_df$Description, width = 27)

# Keep only the top 7 rows if there are more than 10
if (nrow(res_df) > 7) {
  res_df <- head(res_df, 7)
}

p <- ggplot(res_df, aes(x = -log10(p.adjust), y = reorder(Description, -p.adjust), size = foldEnrichment)) +
  geom_point(color = "blue", shape = 16) +  # Dot color and shape
  theme_minimal() +  # Use a minimal theme
  theme(
    text = element_text(family = "Arial"),  # Set font to Arial
    axis.text.y = element_text(size = 15, color = "black"),  # Adjust y-axis text size
    axis.title = element_text(size = 12, color = "black"),  # Adjust axis title size
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    panel.grid = element_blank(),  # Remove background grid
    legend.position = "right",  # Place legend on the right
    plot.title = element_text(size = 15, hjust = 0.67, face="bold"),
    legend.title = element_text(hjust = 0.7),  # Move title to the right
    legend.title.align = 1,  # Align legend title to the right
    legend.margin = margin(r = 20)) +  # Add space on the right side of the legend
  labs(
    x = "-log10(padj)",  # Label for x-axis
    y = NULL,  # No label for y-axis
    title = "GO Biological Processes enriched at 20°C (vs 28°C)"
  ) +
  xlim(1.4,2.5) +
  scale_size_continuous(
    range = c(7, 15),
    guide = guide_legend(
      title = "foldEnrichment",
      override.aes = list(
        shape = 21,
        fill = NA,
        color = "black",
        stroke = 0.5)))

#adjust color accordingly

ggsave("20-28.png", plot = p, width = 14, height = 10, dpi = 300, units = "cm", limitsize = FALSE)
