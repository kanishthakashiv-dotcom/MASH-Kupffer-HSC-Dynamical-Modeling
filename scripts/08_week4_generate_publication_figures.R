    setwd("C:/Users/KANISHTHA/Desktop/new project/MASH_Systems_Immunology")
    pkgs <- c("ggplot2", "dplyr", "tidyr", "patchwork", "igraph", 
        "scales")
    for (p in pkgs) {
        if (!requireNamespace(p, quietly = TRUE)) 
            install.packages(p, repos = "https://cloud.r-project.org")
        library(p, character.only = TRUE)
    }
    dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
    theme_npj <- function(base_size = 10) {
        theme_classic(base_size = base_size) + theme(plot.title = element_text(face = "bold", 
            size = base_size + 1, hjust = 0), plot.subtitle = element_text(size = base_size - 
            1, color = "grey30"), axis.title = element_text(face = "bold", 
            size = base_size), axis.text = element_text(color = "black", 
            size = base_size - 1), legend.title = element_text(face = "bold", 
            size = base_size - 1), legend.text = element_text(size = base_size - 
            1), strip.background = element_rect(fill = "grey94", 
            color = "black", linewidth = 0.5), strip.text = element_text(face = "bold", 
            size = base_size - 1))
    }
    save_pub_fig <- function(plot_obj, filename_base, w = 11, 
        h = 5.5) {
        ggsave(file.path("results/figures", paste0(filename_base, 
            ".png")), plot_obj, width = w, height = h, dpi = 300)
        ggsave(file.path("results/figures", paste0(filename_base, 
            ".pdf")), plot_obj, width = w, height = h)
        cat("Exported:", filename_base, "(300 DPI PNG + Vector PDF)\n")
    }
    u4 <- read.csv("results/tables/All_4_Cohorts_Unified_15Node_Evidence.csv")
    key_genes <- c("TREM2", "SPP1", "TGFB1", "IL1B", "TNF", "IL10", 
        "TLR4", "COL1A1", "COL1A2", "ACTA2", "TIMP1", "MMP2", 
        "PDGFRB", "CCN2", "YAP1", "SMAD3", "PPARG")
    u4_sub <- u4 %>% filter(Gene %in% key_genes)
    f1_long <- u4_sub %>% dplyr::select(Gene, `Bulk: GSE135251` = log2FoldChange_GSE135251, 
        `Bulk: GSE126848` = log2FoldChange_GSE126848, `Bulk: GSE162694` = log2FoldChange_GSE162694, 
        `scRNA Mac (GSE136103)` = Mac_Log2Ratio, `scRNA HSC (GSE136103)` = HSC_Log2Ratio) %>% 
        pivot_longer(-Gene, names_to = "Dataset", values_to = "Log2FC") %>% 
        mutate(Gene = factor(Gene, levels = rev(key_genes)), 
            Dataset = factor(Dataset, levels = c("Bulk: GSE135251", 
                "Bulk: GSE126848", "Bulk: GSE162694", "scRNA Mac (GSE136103)", 
                "scRNA HSC (GSE136103)")))
    p1a <- ggplot(f1_long, aes(x = Dataset, y = Gene, fill = pmax(-3.5, 
        pmin(3.5, Log2FC)))) + geom_tile(color = "white", linewidth = 0.4) + 
        geom_text(aes(label = ifelse(is.na(Log2FC), "NA", sprintf("%.2f", 
            Log2FC))), size = 2.7) + scale_fill_gradient2(low = "#2166AC", 
        mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey90", 
        name = "Log2 FC /\nLog2 Ratio") + labs(title = "A. Cross-Cohort Conserved Dysregulation of Kupffer/SAMac-HSC Circuit", 
        subtitle = "Consistent induction of TREM2/SPP1/TGFB1 & ECM-YAP1 axis alongside IL10 repression across 4 GEO cohorts", 
        x = NULL, y = "Core Network Gene") + theme_npj() + theme(axis.text.x = element_text(angle = 25, 
        hjust = 1))
    sc_val <- read.csv("results/tables/GSE136103_SingleCell_Mac_HSC_Network_Validation.csv") %>% 
        filter(Gene %in% c("TREM2", "SPP1", "TGFB1", "IL10", 
            "COL1A1", "ACTA2", "TIMP1", "MMP2", "PDGFRB", "PPARG"))
    p1b <- ggplot(sc_val, aes(x = Mac_Log2Ratio, y = HSC_Log2Ratio, 
        label = Gene)) + geom_hline(yintercept = 0, linetype = "dashed", 
        color = "grey60") + geom_vline(xintercept = 0, linetype = "dashed", 
        color = "grey60") + geom_point(aes(color = ifelse(Mac_Log2Ratio > 
        0 & HSC_Log2Ratio > 0, "Co-Induced", "Compartment-Specific")), 
        size = 3.5) + geom_text(vjust = -0.8, size = 3.1, fontface = "bold") + 
        scale_color_manual(values = c(`Co-Induced` = "#B2182B", 
            `Compartment-Specific` = "#2166AC"), name = "Compartment Shift") + 
        labs(title = "B. Single-Cell Compartment Partitioning (GSE136103)", 
            subtitle = "Cirrhotic vs. Healthy Log2 expression shift in Macrophages vs. HSCs", 
            x = "Macrophage / SAMac Log2 Shift", y = "HSC / Myofibroblast Log2 Shift") + 
        theme_npj()
    fig1 <- p1a + p1b + plot_layout(widths = c(1.3, 1))
    save_pub_fig(fig1, "Fig1_MultiCohort_Network_Validation", 
        w = 12.5, h = 5.8)
    topo <- read.table("network/MASH_Kupffer_HSC_15Node.topo", 
        header = TRUE, stringsAsFactors = FALSE)
    attr_mat <- read.csv("results/tables/BoolNet_SteadyState_Attractors.csv", 
        row.names = 1)
    png("results/figures/Fig2A_Network_Topology_Raw.png", width = 1800, 
        height = 1600, res = 300)
    par(mar = c(1, 1, 2.5, 1))
    g <- graph_from_data_frame(topo[, c("Source", "Target")], 
        directed = TRUE)
    mac_nodes <- c("TLR4_NFkB", "TNF_IL1B", "TREM2_SPP1", "TGFB1", 
        "IL10_PPARG")
    V(g)$color <- ifelse(V(g)$name %in% mac_nodes, "#4393C3", 
        "#D6604D")
    E(g)$color <- ifelse(topo$Type == 1, "#1B7837", "#B2182B")
    E(g)$lty <- ifelse(topo$Type == 1, 1, 2)
    set.seed(42)
    plot(g, layout = layout_in_circle(g), vertex.size = 24, vertex.label.cex = 0.72, 
        vertex.label.color = "black", vertex.label.font = 2, 
        vertex.frame.color = "black", edge.arrow.size = 0.5, 
        edge.width = 1.6, edge.curved = 0.15, main = "A. Signed Directed Kupffer/SAMac - HSC Regulatory Circuit")
    legend("bottomleft", legend = c("Macrophage/SAMac Module", 
        "HSC/Stellate Module", "Activation (->)", "Inhibition (--|)"), 
        pch = c(19, 19, NA, NA), lty = c(NA, NA, 1, 2), lwd = 2, 
        col = c("#4393C3", "#D6604D", "#1B7837", "#B2182B"), 
        cex = 0.72, bty = "n")
    dev.off()
    attr_df <- as.data.frame(attr_mat)
    attr_df$Attractor <- rownames(attr_df)
    attr_long <- attr_df %>% pivot_longer(-Attractor, names_to = "Node", 
        values_to = "State")
    p2b <- ggplot(attr_long, aes(x = Node, y = Attractor, fill = factor(round(State)))) + 
        geom_tile(color = "white", linewidth = 0.5) + scale_fill_manual(values = c(`0` = "#F7F7F7", 
        `1` = "#B2182B"), labels = c("OFF (0)", "ON (1)"), name = "Boolean State") + 
        labs(title = "B. BoolNet Steady-State Attractors & Basin Sizes", 
            subtitle = "99.5% of state space converges to locked TREM2_SPP1=1 / YAP1_TAZ=1 fibrotic attractors", 
            x = NULL, y = NULL) + theme_npj() + theme(axis.text.x = element_text(angle = 35, 
        hjust = 1))
    save_pub_fig(p2b, "Fig2B_Boolean_Attractors", w = 9.5, h = 4.2)
    ens <- readRDS("data/processed/sRACIPE_MASH_Ensemble.rds")
    ss_mat <- ens$ss_mat
    st_labels <- ens$state_labels
    set.seed(2026)
    sub_idx <- sample(seq_len(ncol(ss_mat)), min(3000, ncol(ss_mat)))
    pca_res <- prcomp(t(ss_mat[, sub_idx]), scale. = FALSE)
    pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 
        2], Attractor = st_labels[sub_idx])
    var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 
        1)
    p3a <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Attractor)) + 
        geom_point(alpha = 0.55, size = 1.3) + scale_color_manual(values = c(State1_Quiescent_Homeostatic = "#2166AC", 
        State2_Reversible_Inflammatory = "#FDB863", State3_Irreversible_Fibrotic = "#B2182B")) + 
        labs(title = "A. Continuous ODE Ensemble PCA (sRACIPE, n = 37,500 trajectories)", 
            subtitle = "Emergence of three discrete phenotypic basins across randomized kinetic parameters", 
            x = paste0("PC1 (", var_exp[1], "% variance)"), y = paste0("PC2 (", 
                var_exp[2], "% variance)")) + theme_npj() + theme(legend.position = "bottom")
    centroids <- read.csv("results/tables/sRACIPE_Attractor_Centroids.csv")
    cent_long <- centroids %>% pivot_longer(-Node, names_to = "Attractor", 
        values_to = "ZScore")
    p3b <- ggplot(cent_long, aes(x = Attractor, y = Node, fill = ZScore)) + 
        geom_tile(color = "white", linewidth = 0.4) + geom_text(aes(label = sprintf("%.2f", 
        ZScore)), size = 2.8) + scale_fill_gradient2(low = "#2166AC", 
        mid = "white", high = "#B2182B", midpoint = 0, name = "Mean Z-Score") + 
        labs(title = "B. Attractor Node Expression Centroids", 
            subtitle = "Antagonistic switch between PPARG/IL10 & YAP1/COL1A1/TREM2", 
            x = NULL, y = NULL) + theme_npj() + theme(axis.text.x = element_text(angle = 20, 
        hjust = 1))
    fig3 <- p3a + p3b + plot_layout(widths = c(1.25, 1))
    save_pub_fig(fig3, "Fig3_sRACIPE_Continuous_Landscape", w = 12.5, 
        h = 5.6)
    hyst <- read.csv("results/tables/Hysteresis_Trigger_Removal_Trajectory.csv")
    hyst_long <- hyst %>% pivot_longer(cols = c(COL1A1_ACTA2, 
        YAP1_TAZ, TREM2_SPP1, PPARG_HSC), names_to = "Module", 
        values_to = "Activity")
    p4 <- ggplot(hyst_long, aes(x = TLR4_Drive, y = Activity, 
        color = Direction, linetype = Direction)) + geom_line(linewidth = 1.1) + 
        geom_point(size = 1.8) + facet_wrap(~Module, scales = "free_y", 
        ncol = 2) + scale_color_manual(values = c(Forward_Injury_Escalation = "#2166AC", 
        Reverse_Trigger_Removal = "#B2182B"), labels = c("Forward Injury Escalation (0 -> 2.0)", 
        "Reverse Trigger Withdrawal (2.0 -> 0)")) + scale_linetype_manual(values = c(Forward_Injury_Escalation = "solid", 
        Reverse_Trigger_Removal = "dashed"), labels = c("Forward Injury Escalation (0 -> 2.0)", 
        "Reverse Trigger Withdrawal (2.0 -> 0)")) + labs(title = "Irreversible Hysteresis Across the Kupffer/SAMac - HSC Regulatory Circuit", 
        subtitle = "Complete withdrawal of upstream TLR4/NFkB drive fails to restore quiescent PPARG_HSC or collapse COL1A1_ACTA2", 
        x = "Upstream Metabolic / Inflammatory Drive (TLR4_NFkB Input)", 
        y = "Steady-State Module Activity (a.u.)") + theme_npj() + 
        theme(legend.position = "bottom")
    save_pub_fig(p4, "Fig4_Hysteresis_Bifurcation_Proof", w = 10.5, 
        h = 6.5)
    pat <- read.csv("results/tables/All_416_Patients_sRACIPE_Attractor_Projection.csv")
    stage_order <- c("Healthy", "Obese_Control", "NAFL_F0", "Early_F1_F2", 
        "Advanced_F3_F4")
    pat$Stage <- factor(pat$Stage, levels = intersect(stage_order, 
        unique(pat$Stage)))
    p5a <- ggplot(pat, aes(x = Stage, y = Net_Fibrotic_Attractor_Score, 
        fill = Stage)) + geom_boxplot(outlier.shape = NA, alpha = 0.85, 
        width = 0.6) + geom_jitter(width = 0.15, size = 1, alpha = 0.45) + 
        facet_wrap(~Cohort, scales = "free_x") + scale_fill_brewer(palette = "RdYlBu", 
        direction = -1) + labs(title = "A. Clinical Attractor Trajectory Across 416 Human Liver Biopsies", 
        subtitle = "GSE135251 (p = 1.86e-11), GSE162694 (p = 1.19e-05), GSE126848 (p = 0.053)", 
        x = "Harmonized Histological Stage", y = "Net Fibrotic Attractor Score") + 
        theme_npj() + theme(axis.text.x = element_text(angle = 30, 
        hjust = 1), legend.position = "none")
    sc_proj <- read.csv("results/tables/GSE136103_SingleCell_Attractor_Scores.csv")
    p5b <- ggplot(sc_proj, aes(x = Cell_State, y = SingleCell_Fibrotic_Attractor_Score, 
        fill = Condition)) + geom_boxplot(outlier.shape = NA, 
        alpha = 0.85) + scale_fill_manual(values = c(Healthy = "#4393C3", 
        Cirrhotic_MASH = "#B2182B")) + labs(title = "B. Single-Cell Niche Shift (GSE136103)", 
        subtitle = "SAMac_TREM2_SPP1 & Kupffer activation in Cirrhotic MASH", 
        x = NULL, y = "Single-Cell Attractor Score") + theme_npj() + 
        theme(axis.text.x = element_text(angle = 25, hjust = 1), 
            legend.position = "bottom")
    fig5 <- p5a + p5b + plot_layout(widths = c(1.65, 1))
    save_pub_fig(fig5, "Fig5_Clinical_and_SingleCell_Projection", 
        w = 13, h = 5.5)
    s_pert <- read.csv("results/tables/InSilico_Single_Node_Perturbations.csv")
    d_pert <- read.csv("results/tables/InSilico_Dual_Node_Synergy_Screen.csv")
    top_single <- s_pert %>% head(10) %>% mutate(Category = "Single-Node Monotherapy")
    top_dual <- d_pert %>% head(10) %>% mutate(Category = "Dual-Node Combination")
    comb_pert <- bind_rows(top_single %>% dplyr::select(Perturbation, 
        Category, Final_COL1A1_ACTA2, Final_PPARG_HSC, Fibrotic_Collapse_Pct), 
        top_dual %>% dplyr::select(Perturbation, Category, Final_COL1A1_ACTA2, 
            Final_PPARG_HSC, Fibrotic_Collapse_Pct))
    p6a <- ggplot(comb_pert, aes(x = reorder(Perturbation, Fibrotic_Collapse_Pct), 
        y = Fibrotic_Collapse_Pct, fill = Category)) + geom_col(width = 0.72) + 
        coord_flip() + scale_fill_manual(values = c(`Single-Node Monotherapy` = "#878787", 
        `Dual-Node Combination` = "#B2182B")) + labs(title = "A. In-Silico Reversal of Locked Fibrotic Attractor (Single vs. Dual Targets)", 
        subtitle = "Comparison of 13 single-node perturbations vs. top 10 dual-node combinations", 
        x = NULL, y = "Fibrotic Attractor Collapse (%)") + theme_npj() + 
        theme(legend.position = "bottom")
    p6b <- ggplot(d_pert, aes(x = Final_COL1A1_ACTA2, y = Final_PPARG_HSC, 
        color = Attractor_Reverted)) + geom_point(size = 2.6, 
        alpha = 0.8) + scale_color_manual(values = c(`FALSE` = "grey60", 
        `TRUE` = "#1B7837"), labels = c("Locked Fibrotic", "Reverted to Quiescent"), 
        name = "Final Attractor") + labs(title = "B. Phase-Plane Attractor Reversion Across All 78 Dual Combinations", 
        subtitle = "Simultaneous suppression of COL1A1/YAP1/ECM & restoration of PPARG_HSC", 
        x = "Post-Treatment COL1A1_ACTA2 Activity", y = "Post-Treatment PPARG_HSC Activity") + 
        theme_npj() + theme(legend.position = "bottom")
    fig6 <- p6a + p6b + plot_layout(widths = c(1.25, 1))
    save_pub_fig(fig6, "Fig6_Perturbation_Synergy_Screen", w = 13, 
        h = 6)
    cat("\nStep 4.1 completed 100% cleanly! All 6 publication figures saved in results/figures/\n")
