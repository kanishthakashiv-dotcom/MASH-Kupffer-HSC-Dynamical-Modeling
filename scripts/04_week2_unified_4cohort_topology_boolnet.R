    setwd("C:/Users/KANISHTHA/Desktop/new project/MASH_Systems_Immunology")
    if (!requireNamespace("BoolNet", quietly = TRUE)) 
        install.packages("BoolNet", repos = "https://cloud.r-project.org")
    library(BoolNet)
    library(dplyr)
    library(tidyr)
    bulk_de <- read.csv("results/tables/All_Bulk_Cohorts_15Node_Network_DE.csv")
    sc_val <- read.csv("results/tables/GSE136103_SingleCell_Mac_HSC_Network_Validation.csv")
    bulk_wide <- bulk_de %>% dplyr::select(Cohort, Gene, log2FoldChange, 
        padj) %>% pivot_wider(names_from = Cohort, values_from = c(log2FoldChange, 
        padj))
    unified_4cohort <- full_join(bulk_wide, sc_val, by = "Gene") %>% 
        arrange(Gene)
    write.csv(unified_4cohort, "results/tables/All_4_Cohorts_Unified_15Node_Evidence.csv", 
        row.names = FALSE)
    cat("=== Unified 4-Cohort Evidence Table Saved (n =", nrow(unified_4cohort), 
        "genes across 4 GEO cohorts) ===\n")
    print(head(unified_4cohort[, c("Gene", "log2FoldChange_GSE135251", 
        "log2FoldChange_GSE126848", "log2FoldChange_GSE162694", 
        "Mac_Log2Ratio", "HSC_Log2Ratio")], 15))
    topo_edges <- data.frame(Source = c("TLR4_NFkB", "TLR4_NFkB", 
        "TNF_IL1B", "TNF_IL1B", "TREM2_SPP1", "IL10_PPARG", "IL10_PPARG", 
        "TGFB1", "TREM2_SPP1", "TGFB1", "TREM2_SPP1", "TNF_IL1B", 
        "TNF_IL1B", "TGFBR_SMAD3", "TGFBR_SMAD3", "TGFBR_SMAD3", 
        "YAP1_TAZ", "YAP1_TAZ", "PDGFRB_ERK", "PDGFRB_ERK", "NFkB_HSC", 
        "TIMP1", "COL1A1_ACTA2", "PPARG_HSC", "PPARG_HSC", "TGFBR_SMAD3", 
        "YAP1_TAZ", "SMAD7", "TGFBR_SMAD3", "IL10_PPARG", "COL1A1_ACTA2"), 
        Target = c("TNF_IL1B", "NFkB_HSC", "TREM2_SPP1", "TLR4_NFkB", 
            "TGFB1", "TLR4_NFkB", "TNF_IL1B", "IL10_PPARG", "IL10_PPARG", 
            "TGFBR_SMAD3", "YAP1_TAZ", "NFkB_HSC", "PDGFRB_ERK", 
            "COL1A1_ACTA2", "TIMP1", "YAP1_TAZ", "COL1A1_ACTA2", 
            "PDGFRB_ERK", "COL1A1_ACTA2", "YAP1_TAZ", "TIMP1", 
            "COL1A1_ACTA2", "YAP1_TAZ", "TGFBR_SMAD3", "YAP1_TAZ", 
            "PPARG_HSC", "PPARG_HSC", "TGFBR_SMAD3", "SMAD7", 
            "PPARG_HSC", "TREM2_SPP1"), Type = c(1, 1, 1, 1, 
            1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
            1, 1, 2, 2, 2, 2, 2, 1, 1, 1), stringsAsFactors = FALSE)
    dir.create("network", showWarnings = FALSE)
    write.table(topo_edges, "network/MASH_Kupffer_HSC_15Node.topo", 
        sep = "\t", row.names = FALSE, quote = FALSE)
    cat("\nSaved signed regulatory topology (.topo) with", nrow(topo_edges), 
        "edges across", length(unique(c(topo_edges$Source, topo_edges$Target))), 
        "nodes.\n")
    all_nodes <- sort(unique(c(topo_edges$Source, topo_edges$Target)))
    rules <- c("targets, factors")
    for (nd in all_nodes) {
        acts <- topo_edges$Source[topo_edges$Target == nd & topo_edges$Type == 
            1]
        inhs <- topo_edges$Source[topo_edges$Target == nd & topo_edges$Type == 
            2]
        act_str <- if (length(acts) > 0) 
            paste0("(", paste(acts, collapse = " | "), ")")
        else ""
        inh_str <- if (length(inhs) > 0) 
            paste0("!(", paste(inhs, collapse = " | "), ")")
        else ""
        rhs <- paste(c(act_str, inh_str)[c(act_str, inh_str) != 
            ""], collapse = " & ")
        if (rhs == "") 
            rhs <- nd
        rules <- c(rules, paste0(nd, ", ", rhs))
    }
    writeLines(rules, "network/MASH_Kupffer_HSC_15Node.boolnet")
    net <- loadNetwork("network/MASH_Kupffer_HSC_15Node.boolnet")
    attr_sync <- getAttractors(net, type = "synchronous")
    attr_df <- getAttractorSequence(attr_sync, 1)
    attr_list <- lapply(seq_along(attr_sync$attractors), function(i) {
        seq_mat <- getAttractorSequence(attr_sync, i)
        colMeans(seq_mat)
    })
    attr_mat <- do.call(rbind, attr_list)
    rownames(attr_mat) <- paste0("Attractor_", seq_len(nrow(attr_mat)), 
        " (Basin=", round(sapply(attr_sync$attractors, function(a) a$basinSize)/(2^length(all_nodes)) * 
            100, 1), "%)")
    write.csv(attr_mat, "results/tables/BoolNet_SteadyState_Attractors.csv", 
        row.names = TRUE)
    saveRDS(list(net = net, attractors = attr_sync, topo = topo_edges), 
        "data/processed/BoolNet_MASH_Model.rds")
    cat("\n=== Boolean Steady-State Attractors Discovered ===\n")
    print(round(attr_mat, 2))
    cat("\nStep 2.1 completed 100% cleanly! All 4 cohorts merged and Boolean attractors identified.\n")
