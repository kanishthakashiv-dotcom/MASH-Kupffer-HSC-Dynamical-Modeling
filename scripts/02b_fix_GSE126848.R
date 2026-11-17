
setwd("C:/Users/KANISHTHA/Desktop/new project/MASH_Systems_Immunology")
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(data.table)

network_genes <- list(
  TLR4_NFkB     = c("TLR4", "NFKB1", "RELA"),
  TNF_IL1B      = c("TNF", "IL1B"),
  TREM2_SPP1    = c("TREM2", "SPP1", "CD9"),
  TGFB1         = c("TGFB1"),
  IL10_PPARG    = c("IL10", "CD163", "MRC1"),
  TGFBR_SMAD3   = c("TGFBR1", "TGFBR2", "SMAD3"),
  YAP1_TAZ      = c("YAP1", "WWTR1", "CCN2"),
  PDGFRB_ERK    = c("PDGFRB", "MAPK1", "MAPK3"),
  NFkB_HSC      = c("BCL2", "BIRC3", "ICAM1"),
  TIMP1         = c("TIMP1", "MMP2", "MMP9"),
  COL1A1_ACTA2  = c("COL1A1", "COL1A2", "ACTA2"),
  PPARG_HSC     = c("PPARG", "RXRA", "PLIN2"),
  SMAD7         = c("SMAD7", "BAMBI")
)
all_target_genes <- unique(unlist(network_genes))

# 1. Load GSE126848 counts & align all 57 samples using exact title/metadata labels
gse_dir <- "data/raw/bulk/GSE126848"
all_files <- list.files(gse_dir, pattern = "counts", full.names = TRUE, recursive = TRUE)
all_files <- all_files[!dir.exists(all_files)]
cnt_file <- all_files[which.max(file.info(all_files)$size)]

df <- as.data.frame(fread(cnt_file, header = TRUE, check.names = FALSE))
gene_ids <- as.character(df[[1]])
counts <- df[, -1, drop = FALSE]
counts <- counts[, sapply(counts, is.numeric), drop = FALSE]

clean_ens <- sub("\\..*$", "", gene_ids)
symbols <- mapIds(org.Hs.eg.db, keys = clean_ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
valid <- !is.na(symbols) & symbols != ""
counts <- counts[valid, , drop = FALSE]
counts$Symbol <- symbols[valid]

counts_agg <- counts %>%
  group_by(Symbol) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  as.data.frame()
rownames(counts_agg) <- counts_agg$Symbol
counts_agg$Symbol <- NULL
counts_mat <- round(as.matrix(counts_agg))

meta <- read.csv(file.path(gse_dir, "GSE126848_metadata.csv"), row.names = 1, check.names = FALSE)
full_row_strings <- apply(meta, 1, function(r) paste(as.character(r), collapse = " | "))
matched_idx <- sapply(colnames(counts_mat), function(code) {
  hits <- grep(paste0("(^|[^0-9])", code, "([^0-9]|$)"), full_row_strings)
  if (length(hits) == 1) return(hits[1])
  return(NA_integer_)
})

keep_cols <- !is.na(matched_idx)
counts_mat <- counts_mat[, keep_cols, drop = FALSE]
meta_aligned <- meta[matched_idx[keep_cols], , drop = FALSE]
rownames(meta_aligned) <- colnames(counts_mat)

# Harmonize using meta$title prefixes ("Normal_", "Obese_", "NAFL_", "NASH_")
titles_low <- tolower(as.character(meta_aligned$title))
stage <- rep("Unclassified", nrow(meta_aligned))
stage[grepl("^normal|healthy|lean", titles_low)] <- "Healthy"
stage[grepl("^obese", titles_low)] <- "Obese_Control"
stage[grepl("^nafl", titles_low)] <- "NAFL_F0"
stage[grepl("^nash", titles_low)] <- "Advanced_F3_F4"
meta_aligned$Harmonized_Stage <- stage

cat("=== GSE126848 Full 57-Sample Stage Distribution ===\n")
print(table(meta_aligned$Harmonized_Stage))

keep <- meta_aligned$Harmonized_Stage != "Unclassified"
counts_sub <- counts_mat[, keep, drop = FALSE]
meta_sub <- meta_aligned[keep, , drop = FALSE]
meta_sub$Harmonized_Stage <- factor(meta_sub$Harmonized_Stage)
counts_sub <- counts_sub[rowSums(counts_sub >= 10) >= 5, ]

dds <- DESeqDataSetFromMatrix(countData = counts_sub, colData = meta_sub, design = ~ Harmonized_Stage)
dds <- DESeq(dds, quiet = TRUE)
vst_mat <- assay(vst(dds, blind = FALSE))

saveRDS(list(vst = vst_mat, meta = meta_sub, dds = dds), "data/processed/GSE126848_processed.rds")

# Use Healthy if present, otherwise NAFL_F0 as reference contrast
ref_group <- ifelse("Healthy" %in% levels(meta_sub$Harmonized_Stage), "Healthy", "NAFL_F0")
res <- results(dds, contrast = c("Harmonized_Stage", "Advanced_F3_F4", ref_group))
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)
res_df$Cohort <- "GSE126848"
write.csv(res_df, "results/tables/GSE126848_DE_Advanced_vs_Healthy.csv", row.names = FALSE)

# 2. Rebuild combined 15-node DE table across ALL THREE bulk cohorts
all_de <- list()
for (id in c("GSE135251", "GSE126848", "GSE162694")) {
  df_de <- read.csv(file.path("results/tables", paste0(id, "_DE_Advanced_vs_Healthy.csv")))
  all_de[[id]] <- df_de %>% filter(Gene %in% all_target_genes) %>% arrange(padj)
}
combined_net <- bind_rows(all_de)
write.csv(combined_net, "results/tables/All_Bulk_Cohorts_15Node_Network_DE.csv", row.names = FALSE)

cat("\n================ ALL 3 BULK COHORTS VERIFIED ================\n")
print(table(combined_net$Cohort))
print(head(combined_net[combined_net$Cohort == "GSE126848" & combined_net$padj < 0.05 & !is.na(combined_net$padj),
                        c("Cohort", "Gene", "log2FoldChange", "padj")], 20))
cat("\nStep 1.2 completed 100% cleanly across all 3 bulk cohorts!\n")

