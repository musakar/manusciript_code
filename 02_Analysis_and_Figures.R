# -------------------------------------------------- #
# Script 02: Main Analysis and Figure Generation
# Purpose: Perform data pre-processing, statistical analysis, and
#          generate all figures and tables for the manuscript.
# -------------------------------------------------- #

# SECTION 0: LOAD LIBRARIES
# --------------------------------------------------
cat("Loading required packages...\n")
library(affy)
library(Biobase)
library(sva)
library(limma)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(UpSetR)
library(clusterProfiler)
library(org.At.tair.db)
library(ath1121501.db)
library(AnnotationDbi)
library(knitr)

# SECTION 1: DEFINE PATHS & LOAD DATA
# --------------------------------------------------
cat("\nSECTION 1: Loading and Preparing Data\n")

# Define project folder structure
base_path <- "C:/Users/Musa KAR/Desktop/yeniarabophis"
data_path <- file.path(base_path, "data")
results_folder <- file.path(base_path, "results")

# Create results directory if it doesn't exist
if (!dir.exists(results_folder)) {
  dir.create(results_folder)
}

# Read raw .CEL files
raw_data <- ReadAffy(celfile.path = data_path)

# Define master metadata for all 54 samples
master_metadata <- tribble(
  ~filename, ~GSE, ~pesticide, ~condition,
  "GSM264761.CEL", "GSE10464", "Control", "Control", "GSM264762.CEL", "GSE10464", "Control", "Control", "GSM264763.CEL", "GSE10464", "Control", "Control", "GSM264764.CEL", "GSE10464", "PQT", "Treated", "GSM264765.CEL", "GSE10464", "PQT", "Treated", "GSM264766.CEL", "GSE10464", "PQT", "Treated",
  "GSM226328.CEL", "GSE8927", "Control", "Control", "GSM226329.CEL", "GSE8927", "Control", "Control", "GSM226330.CEL", "GSE8927", "Control", "Control", "GSM226331.CEL", "GSE8927", "GLPS", "Treated", "GSM226332.CEL", "GSE8927", "GLPS", "Treated", "GSM226333.CEL", "GSE8927", "GLPS", "Treated",
  "GSM226322.CEL", "GSE8926", "Control", "Control", "GSM226323.CEL", "GSE8926", "Control", "Control", "GSM226324.CEL", "GSE8926", "Control", "Control", "GSM226325.CEL", "GSE8926", "TYZ", "Treated", "GSM226326.CEL", "GSE8926", "TYZ", "Treated", "GSM226327.CEL", "GSE8926", "TYZ", "Treated",
  "GSM226316.CEL", "GSE8925", "Control", "Control", "GSM226317.CEL", "GSE8925", "Control", "Control", "GSM226318.CEL", "GSE8925", "Control", "Control", "GSM226319.CEL", "GSE8925", "IDZ", "Treated", "GSM226314.CEL", "GSE8925", "IDZ", "Treated", "GSM226315.CEL", "GSE8925", "IDZ", "Treated",
  "GSM225909.CEL", "GSE8912", "Control", "Control", "GSM225910.CEL", "GSE8912", "Control", "Control", "GSM225911.CEL", "GSE8912", "Control", "Control", "GSM225912.CEL", "GSE8912", "SMZ", "Treated", "GSM225913.CEL", "GSE8912", "SMZ", "Treated", "GSM225914.CEL", "GSE8912", "SMZ", "Treated",
  "GSM225903.CEL", "GSE8913", "Control", "Control", "GSM225904.CEL", "GSE8913", "Control", "Control", "GSM225905.CEL", "GSE8913", "Control", "Control", "GSM225906.CEL", "GSE8913", "PSF", "Treated", "GSM225907.CEL", "GSE8913", "PSF", "Treated", "GSM225908.CEL", "GSE8913", "PSF", "Treated",
  "GSM591834.CEL", "GSE24052", "Control", "Control", "GSM591835.CEL", "GSE24052", "Control", "Control", "GSM591836.CEL", "GSE24052", "Control", "Control", "GSM591837.CEL", "GSE24052", "DCB", "Treated", "GSM591838.CEL", "GSE24052", "DCB", "Treated", "GSM591839.CEL", "GSE24052", "DCB", "Treated",
  "GSM1259918_R1_AVA098_wtC_NP.CEL", "GSE52117", "Control", "Control", "GSM1259930_R2_AVA172_wtC_NP.CEL", "GSE52117", "Control", "Control", "GSM1259939_R3_AVA194_wtC_NP.CEL", "GSE52117", "Control", "Control", "GSM1259952_T_R1_AVA097_wtC.CEL",   "GSE52117", "CHS", "Treated", "GSM1259956_T_R2_AVA171_wtC.CEL",   "GSE52117", "CHS", "Treated", "GSM1259960_T_R3_AVA197_wtC.CEL",   "GSE52117", "CHS", "Treated",
  "GSM1257538_DMSO_1.CEL", "GSE52013", "Control", "Control", "GSM1257539_DMSO_2.CEL", "GSE52013", "Control", "Control", "GSM1257540_DMSO_3.CEL", "GSE52013", "Control", "Control", "GSM1257541_Treated_1.CEL", "GSE52013", "auxg", "Treated", "GSM1257542_Treated_2.CEL", "GSE52013", "auxg", "Treated", "GSM1257543_Treated_3.CEL", "GSE52013", "auxg", "Treated"
)

# Attach phenoData to AffyBatch object
actual_sample_names <- sampleNames(raw_data)
phenoData_ordered <- master_metadata %>% filter(filename %in% actual_sample_names) %>% slice(match(actual_sample_names, filename))
phenoData_ordered <- as.data.frame(phenoData_ordered)
rownames(phenoData_ordered) <- phenoData_ordered$filename
pdata_final <- AnnotatedDataFrame(data = phenoData_ordered)
phenoData(raw_data) <- pdata_final


# SECTION 2: NORMALIZATION & BATCH CORRECTION
# --------------------------------------------------
cat("\nSECTION 2: Normalizing Data and Correcting Batch Effects\n")
# Perform RMA normalization
eset <- rma(raw_data)
# Correct for batch effects using ComBat
batch_info <- pData(eset)$GSE
batch_corrected_matrix <- sva::ComBat(dat = exprs(eset), batch = batch_info)


# SECTION 3: GENERATE FIGURE 1 (PCA PLOTS)
# --------------------------------------------------
cat("\nSECTION 3: Generating Figure 1 (PCA Plots)...\n")
# Prepare data for Fig 1A
pca_data_before <- prcomp(t(exprs(eset)))
pca_df_before <- as.data.frame(pca_data_before$x) %>% rownames_to_column(var = "filename") %>% left_join(pData(eset), by = "filename")
percent_variance_before <- round(100 * pca_data_before$sdev^2 / sum(pca_data_before$sdev^2), 1)
# Create Fig 1A
pca_plot_before <- ggplot(pca_df_before, aes(x = PC1, y = PC2, color = GSE, shape = condition)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(title = "Before Batch Correction", x = paste0("PC1: ", percent_variance_before[1], "% variance"), y = paste0("PC2: ", percent_variance_before[2], "% variance")) +
  theme_bw(base_size = 11) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") 
# Prepare data for Fig 1B
pca_data_corrected <- prcomp(t(batch_corrected_matrix))
pca_df_corrected <- as.data.frame(pca_data_corrected$x) %>% rownames_to_column(var = "filename") %>% left_join(pData(eset), by = "filename")
percent_variance_corrected <- round(100 * pca_data_corrected$sdev^2 / sum(pca_data_corrected$sdev^2), 1)
my_colors <- c("Control" = "#000000", "PQT" = "#E69F00", "GLPS" = "#56B4E9", "TYZ" = "#009E73", "IDZ" = "#F0E442", "SMZ" = "#0072B2", "PSF" = "#D55E00", "DCB" = "#CC79A7", "CHS" = "red", "auxg" = "purple")
# Create Fig 1B
pca_plot_corrected <- ggplot(pca_df_corrected, aes(x = PC1, y = PC2, color = pesticide, shape = condition)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = my_colors) +
  labs(title = "After Batch Correction", x = paste0("PC1: ", percent_variance_corrected[1], "% variance"), y = paste0("PC2: ", percent_variance_corrected[2], "% variance"), color = "Pesticide", shape = "Condition") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# Combine and save Figure 1
combined_pca_plot <- pca_plot_before + pca_plot_corrected + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom')
ggsave(filename = "Figure1_PCA_Combined.png", plot = combined_pca_plot, path = results_folder, width = 14, height = 7, dpi = 300)
cat("Figure 1 saved to 'results' folder.\n")


# SECTION 4: LIMMA DEG ANALYSIS & TABLE 1
# --------------------------------------------------
cat("\nSECTION 4: Performing DEG Analysis and Generating Table 1...\n")
# Define design and contrast matrices
group <- factor(pData(eset)$pesticide)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
contrast_matrix <- makeContrasts(PQT_vs_Control = PQT - Control, GLPS_vs_Control = GLPS - Control, TYZ_vs_Control = TYZ - Control, IDZ_vs_Control = IDZ - Control, SMZ_vs_Control = SMZ - Control, PSF_vs_Control = PSF - Control, DCB_vs_Control = DCB - Control, CHS_vs_Control = CHS - Control, auxg_vs_Control = auxg - Control, levels = design)
# Fit linear model and apply contrasts
fit <- lmFit(batch_corrected_matrix, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
# Generate and save DEG summary table (Table 1)
all_results_summary <- summary(decideTests(fit2))
summary_table <- as.data.frame(t(all_results_summary))
colnames(summary_table) <- c("Downregulated", "Not Significant", "Upregulated")
write.csv(summary_table, file = file.path(results_folder, "Table1_DEG_Summary.csv"), row.names = TRUE)
cat("Table 1 saved to 'results' folder.\n")


# SECTION 5: VOLCANO PLOTS (FIGURE 2 & SUPP. FIG S1)
# --------------------------------------------------
cat("\nSECTION 5: Generating Volcano Plots...\n")
# Define a reusable function to generate volcano plots
create_volcano_plot <- function(fit_object, coef_name, pesticide_name) {
  results <- topTable(fit_object, coef = coef_name, n = Inf) %>% mutate(significance = case_when(adj.P.Val < 0.05 & logFC > 1  ~ "Upregulated", adj.P.Val < 0.05 & logFC < -1 ~ "Downregulated", TRUE ~ "Not Significant")) %>% rownames_to_column(var = "gene_id")
  top_genes_to_label <- results %>% filter(significance != "Not Significant") %>% slice_max(order_by = B, n = 10)
  volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = significance), alpha = 0.7, size = 2) +
    scale_color_manual(name = "Significance", values = c("Downregulated" = "#0072B2", "Upregulated" = "#D55E00", "Not Significant" = "grey80")) +
    geom_text_repel(data = top_genes_to_label, aes(label = gene_id), size = 3.5, box.padding = 0.5, max.overlaps = Inf, color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    labs(title = paste("Volcano Plot of DEGs for", pesticide_name, "Treatment"), subtitle = paste(pesticide_name, "vs. Control"), x = expression(Log[2]*" (Fold Change)"), y = expression(-Log[10]*" (Adjusted P-Value)")) +
    theme_bw(base_size = 14) + theme(legend.position = "top", plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5))
  return(volcano_plot)
}
# Generate and save Figure 2 (SMZ)
volcano_plot_smz <- create_volcano_plot(fit_object = fit2, coef_name = "SMZ_vs_Control", pesticide_name = "SMZ")
ggsave(filename = "Figure2_Volcano_SMZ.png", plot = volcano_plot_smz, path = results_folder, width = 10, height = 8, dpi = 300)
# Generate and save Supplementary Figure S1
volcano_plot_psf <- create_volcano_plot(fit_object = fit2, coef_name = "PSF_vs_Control", pesticide_name = "PSF")
ggsave(filename = "Supplementary_Figure_S1A_Volcano_PSF.png", plot = volcano_plot_psf, path = results_folder, width = 10, height = 8, dpi = 300)
volcano_plot_dcb <- create_volcano_plot(fit_object = fit2, coef_name = "DCB_vs_Control", pesticide_name = "DCB")
ggsave(filename = "Supplementary_Figure_S1B_Volcano_DCB.png", plot = volcano_plot_dcb, path = results_folder, width = 10, height = 8, dpi = 300)
cat("Figure 2 and Supplementary Figure S1 saved to 'results' folder.\n")


# SECTION 6: UPSET PLOT (FIGURE 3)
# --------------------------------------------------
cat("\nSECTION 6: Generating Figure 3 (UpSet Plot)...\n")
# Get DEG lists for intersection
degs_smz <- topTable(fit2, coef = "SMZ_vs_Control", n = Inf, p.value = 0.05)
degs_psf <- topTable(fit2, coef = "PSF_vs_Control", n = Inf, p.value = 0.05)
degs_tyz <- topTable(fit2, coef = "TYZ_vs_Control", n = Inf, p.value = 0.05)
deg_list_for_upset <- list(SMZ = rownames(degs_smz), PSF = rownames(degs_psf), TYZ = rownames(degs_tyz))
# Generate and save UpSet plot to a PNG file
png(file.path(results_folder, "Figure3_UpSet_Plot.png"), width = 1000, height = 700)
upset(fromList(deg_list_for_upset), nsets = 3, order.by = "freq", text.scale = 1.5, point.size = 3.5, line.size = 2)
dev.off()
cat("Figure 3 saved to 'results' folder.\n")


# SECTION 7: FUNCTIONAL ENRICHMENT (FIGURE 4, TABLE 2, SUPP. FIG S2)
# --------------------------------------------------
cat("\nSECTION 7: Performing Functional Enrichment Analysis...\n")
# Define core gene set and map IDs
core_response_genes_probes <- intersect(rownames(degs_smz), intersect(rownames(degs_psf), rownames(degs_tyz)))
tair_ids <- mapIds(ath1121501.db, keys = core_response_genes_probes, column = "TAIR", keytype = "PROBEID", multiVals = "first")
tair_ids <- na.omit(tair_ids)
# GO Analysis and Figure 4
go_enrichment <- enrichGO(gene = tair_ids, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
if (!is.null(go_enrichment) && nrow(go_enrichment) > 0) {
  go_dot_plot <- dotplot(go_enrichment, showCategory = 20) + labs(title = "GO Enrichment for Core Response Genes (BP)")
  ggsave(filename = "Figure4_GO_DotPlot.png", plot = go_dot_plot, path = results_folder, width = 10, height = 8, dpi = 300)
  cat("Figure 4 saved to 'results' folder.\n")
  # Table 2
  table2_data <- as.data.frame(go_enrichment) %>% slice_head(n = 15) %>% select(ID, Description, GeneRatio, p.adjust, Count) %>% rename("GO Term ID" = ID, "Biological Process" = Description, "Gene Ratio" = GeneRatio, "Adjusted P-Value" = p.adjust, "Gene Count" = Count)
  write.csv(table2_data, file = file.path(results_folder, "Table2_Top_GO_Terms.csv"), row.names = FALSE)
  cat("Table 2 saved to 'results' folder.\n")
} else {cat("No significant GO terms found.\n")}
# KEGG Analysis and Supplementary Figure S2
kegg_entrez_ids <- mapIds(org.At.tair.db, keys = tair_ids, column = "ENTREZID", keytype = "TAIR", multiVals = "first")
kegg_entrez_ids <- na.omit(kegg_entrez_ids)
kegg_enrichment <- enrichKEGG(gene = kegg_entrez_ids, organism = 'ath', pvalueCutoff = 0.05)
if (!is.null(kegg_enrichment) && nrow(kegg_enrichment) > 0) {
  kegg_cnet_plot <- cnetplot(kegg_enrichment, showCategory = 10, node_label = "gene") + labs(title = "KEGG Pathway Enrichment Network")
  ggsave(filename = "Supplementary_Figure_S2_KEGG_Network.png", plot = kegg_cnet_plot, path = results_folder, width = 12, height = 10, dpi = 300)
  cat("Supplementary Figure S2 saved to 'results' folder.\n")
} else {cat("No significant KEGG pathways found.\n")}


# SECTION 8: FINAL WORKSPACE (OPTIONAL)
# --------------------------------------------------
cat("\nSaving the final R workspace image...\n")
save.image(file = file.path(results_folder, "Full_Analysis_Workspace.RData"))
cat("All processes are complete! Check the 'results' folder for your outputs.\n")
