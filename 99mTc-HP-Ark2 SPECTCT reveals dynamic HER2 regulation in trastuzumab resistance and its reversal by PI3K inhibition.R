read_gz_file <- function(filename) {
  data <- read.table(gzfile(filename), header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
  return(data)
}
merge_data_frames <- function(file_list) {
  merged_data <- read_gz_file(file_list[1])
  
  for (i in 2:length(file_list)) {
    next_data <- read_gz_file(file_list[i])
    merged_data <- merge(merged_data, next_data, by = "row.names", all = TRUE)
    rownames(merged_data) <- merged_data$Row.names
    merged_data$Row.names <- NULL
  }
  
  return(merged_data)
}
file_list <- list.files(pattern = "*.txt.gz")
combined_data <- merge_data_frames(file_list)
library(tibble)
combined_data <- rownames_to_column(combined_data, var = "ENSG")
library(AnnotationDbi)
library(org.Hs.eg.db)

gene_ids <- combined_data[,1]


gene_symbols <- mapIds(
  org.Hs.eg.db,                   
  keys = gene_ids,                
  column = "SYMBOL",              
  keytype = "ENSEMBL",            
  multiVals = "first"             
)
gene_conversion <- data.frame(
  ensembl_gene_id = gene_ids,  
  hgnc_symbol = gene_symbols,  
  stringsAsFactors = FALSE
)
combined_data_with_symbols <- merge(
  combined_data,                
  gene_conversion,               
  by.x = "ENSG",            
  by.y = "ensembl_gene_id",      
  all.x = TRUE                  
)
colnames <- names(combined_data_with_symbols)
df <- combined_data_with_symbols[, c("hgnc_symbol", colnames[colnames != "hgnc_symbol"])]
df <- df[,-2]
keep <- rowSums(df >= 10) >= 3
filtered_data <- df[keep, ]
gene_names <- df[, 1]
expression_data <- df[, -1]
total_expression <- rowSums(expression_data)
df_with_total <- data.frame(Gene = gene_names, TotalExpression = total_expression, expression_data)
library(dplyr)
df_dedup <- df_with_total %>%
  group_by(Gene) %>%
  filter(TotalExpression == max(TotalExpression)) %>%
  slice_head(n = 1) %>%  
  ungroup()
df_dedup <- df_dedup[, -2]
keep_genes <- rowSums(df_dedup[, -1] >= 10) >= 3
df <- df_dedup[keep_genes, ]
print(df)
saveRDS(df, file = "df.rds")
write.table(df, file = "df.txt", sep = "\t", quote = FALSE, col.names = NA)
df <- df[complete.cases(df$Gene), ]
df <- column_to_rownames(df, var = "Gene")
sample_info <- data.frame(
  condition = factor(c(rep("SKBR3", 3), rep("SKBR3_HR", 3)))
)
rownames(sample_info) <- colnames(df)
sum(is.na(df))
df <- df[complete.cases(df), ]
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = df,
  colData = sample_info,
  design = ~ condition
)
dds <- DESeq(dds)
results_table <- results(dds)
results_table_df <- as.data.frame(results_table)
write.csv(results_table_df, file = "DESeq2_results.csv", row.names = TRUE)
significant_genes <- subset(results_table_df, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(significant_genes, file = "DESeq2_significant_genes.csv", row.names = TRUE)
library(EnhancedVolcano)
results_table_df <- as.data.frame(results_table)
library(ggplot2)
highlight <- data.frame(
    gene = c("ERBB2","PGP","ABCC1","CASP3","PCNA","VEGFA"),
    log2FC = c(0.6113155, 0.6009124, 0.5650866, -0.7745263, 0.7221752, 0.2331025),
    log10_pvalue = c(4.3018660, 26.0978812, 0.5650866, 13.8821536, 0.2056663, 0.5966332)
)
plot2 <- ggplot(data, aes(x = log2FC, y = log10_pvalue)) +
    geom_point(color = "blue") +
    geom_segment(data = highlight,aes(x = log2FC, y = log10_pvalue, 
                     xend = ifelse(highlight$log2FC < 0, log2FC - 1, log2FC + 1), 
                     yend = log10_pvalue + 5),
                 color = ifelse(highlight$log2FC < 0, "red", "green"),
                 linetype = "dashed") +
    geom_text(data = highlight,aes(label = gene, x = ifelse(highlight$log2FC < 0, log2FC - 1.5, log2FC + 1.5), 
                  y = log10_pvalue + 5),
              vjust = ifelse(highlight$log2FC < 0, 1.5, -0.5),
              color = ifelse(highlight$log2FC < 0, "red", "green")) +
    theme_minimal() +
    labs(title = "Volcano plot", x = "log2FC", y = "-log10 p-value")

print(plot2)
ggsave("volcano_plot.pdf", plot = plot2, width = 8, height = 6)
ggsave("volcano_plot_cell_style.pdf", plot = volcano_plot, device = "pdf", width = 8, height = 6)
row_names <- rownames(results_table_df)
second_column <- results_table_df[, 2]
new_df <- data.frame(RowNames = row_names, SecondColumn = second_column)
gsym.fc <- new_df
colnames(gsym.fc) <-c("SYMBOL","logFC")
dim(gsym.fc)
head(gsym.fc)
gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
gsym.fc.id <- merge(gsym.fc, gsym.id, by="SYMBOL", all=F)
head(gsym.fc.id)
gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T),]
head(gsym.fc.id.sorted)
id.fc <- gsym.fc.id.sorted$logFC
names(id.fc) <- gsym.fc.id.sorted$ENTREZID
head(id.fc)
kk <- gseKEGG(id.fc, organism = "hsa")
dim(kk)
kk.gsym <- setReadable(kk, 'org.Hs.eg.db', 
                       'ENTREZID')
sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T),]
head(sortkk[1:5])
write.csv(sortkk, file = "GSEA_gene.csv", row.names = TRUE)
geneSetID <- c("hsa04668")
selectedGeneID <- c("TNF", "AKT3")
mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
base_size = 12
rankmode <- "sep" 
data <- read.csv("DESeq2_results.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
results_table_df <- data
results_table_df$log2BaseMean <- log2(results_table_df$baseMean + 1)
gene_list <- c("ERBB2","PGP","ABCC1","CASP3","PCNA","VEGFA")
subset_data <- results_table_df[rownames(results_table_df) %in% gene_list, ]
ma_plot <- ggplot(results_table_df, aes(x = log2BaseMean, y = log2FoldChange)) +
    geom_point(alpha = 0.5, color = "blue") +  
    geom_point(data = subset_data, aes(x = log2BaseMean, y = log2FoldChange), 
               color = "red", size = 3) +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
    geom_text_repel(data = subset_data,
                    aes(label = rownames(subset_data)),
                    color = "darkgreen",  
                    size = 5,
                    box.padding = 0.5,  
                    point.padding = 0.5,  
                    segment.color = "darkred",  
                    segment.size = 0.5) +  
    theme_minimal() +
    labs(title = "MA plot", x = "log2(Mean Expression)", y = "log2(Fold Change)") + 
    theme(legend.position = "none")
print(ma_plot)
ggsave("ma_plot_with_labels.pdf", plot = ma_plot, device = "pdf", width = 8, height = 6)
library(ggplot2)
library(ggpubr)
gene_name <- "ERBB2"
normalized_counts <- counts(dds, normalized=TRUE)
gene_expression <- normalized_counts[rownames(normalized_counts) == gene_name, ]
plot_data <- data.frame(
  expression = as.numeric(gene_expression),
  condition = colData(dds)$condition
)
custom_colors <- c("#E69F00", "#56B4E9")
p <- ggplot(plot_data, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +  
  geom_jitter(width = 0.2, size = 1.5, color = "black", alpha = 0.7) +  
  scale_fill_manual(values = custom_colors) +  
  labs(title = paste("Differential Expression of", gene_name),
       x = "Condition",
       y = "Normalized Expression") +
  theme_minimal(base_size = 14) +  
  theme(
    legend.position = "none",  
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    axis.title = element_text(face = "bold"),  
    axis.text = element_text(color = "black") 
  ) +
  stat_compare_means(method = "t.test", label = "p.format", label.y = max(plot_data$expression) * 1.1, 
                     size = 5, label.sep = "\n", color = "red") +  
  geom_signif(comparisons = list(c("Condition1", "Condition2")), 
              map_signif_level=TRUE, textsize=4)  
ggsave("Differential_Expression_Boxplot.pdf", plot = p, width = 8, height = 6)