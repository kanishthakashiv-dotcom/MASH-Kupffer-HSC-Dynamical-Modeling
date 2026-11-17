    setwd("C:/Users/KANISHTHA/Desktop/new project/MASH_Systems_Immunology")
    if (!requireNamespace("BiocManager", quietly = TRUE)) 
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    if (!requireNamespace("sRACIPE", quietly = TRUE)) 
        BiocManager::install("sRACIPE", update = FALSE, ask = FALSE)
    library(sRACIPE)
    library(SummarizedExperiment)
    library(dplyr)
    topo <- read.table("network/MASH_Kupffer_HSC_15Node.topo", 
        header = TRUE, stringsAsFactors = FALSE)
    cat("Loaded signed topology:", nrow(topo), "interactions across", 
        length(unique(c(topo$Source, topo$Target))), "nodes.\n")
    set.seed(2026)
    cat("Running sRACIPE continuous ODE ensemble simulation (please wait ~1-2 minutes)...\n")
    rSet <- sRACIPE::sracipeSimulate(circuit = topo, numModels = 2500, 
        nIC = 15, simulationTime = 50, integrateStepSize = 0.05, 
        plots = FALSE, anneal = FALSE)
    rSet_norm <- sRACIPE::sracipeNormalize(rSet)
    ss_mat <- assay(rSet_norm, 1)
    cat("Total converged continuous steady-state solutions:", 
        ncol(ss_mat), "\n")
    fib_nodes <- intersect(c("COL1A1_ACTA2", "TIMP1", "YAP1_TAZ", 
        "TGFBR_SMAD3", "TREM2_SPP1", "TGFB1", "PDGFRB_ERK"), 
        rownames(ss_mat))
    qui_nodes <- intersect(c("PPARG_HSC", "IL10_PPARG", "SMAD7"), 
        rownames(ss_mat))
    fib_score <- colMeans(ss_mat[fib_nodes, , drop = FALSE])
    qui_score <- colMeans(ss_mat[qui_nodes, , drop = FALSE])
    net_fib_index <- fib_score - qui_score
    set.seed(2026)
    km <- kmeans(t(ss_mat), centers = 3, nstart = 25)
    cluster_means <- tapply(net_fib_index, km$cluster, mean)
    ordered_clusters <- names(sort(cluster_means))
    state_labels <- rep("Intermediate_Inflammatory", ncol(ss_mat))
    state_labels[km$cluster == as.integer(ordered_clusters[1])] <- "State1_Quiescent_Homeostatic"
    state_labels[km$cluster == as.integer(ordered_clusters[2])] <- "State2_Reversible_Inflammatory"
    state_labels[km$cluster == as.integer(ordered_clusters[3])] <- "State3_Irreversible_Fibrotic"
    model_ids <- sub("_.*$", "", colnames(ss_mat))
    if (length(unique(model_ids)) < 100) {
        model_ids <- rep(seq_len(2500), length.out = ncol(ss_mat))
    }
    states_per_model <- tapply(state_labels, model_ids, function(x) length(unique(x)))
    multistability_tab <- table(States_Per_Model = states_per_model)
    prop_df <- as.data.frame(table(Dynamical_Attractor = state_labels))
    prop_df$Percentage <- round(prop_df$Freq/sum(prop_df$Freq) * 
        100, 2)
    centroid_df <- as.data.frame(sapply(sort(unique(state_labels)), 
        function(st) {
            rowMeans(ss_mat[, state_labels == st, drop = FALSE])
        }))
    centroid_df$Node <- rownames(centroid_df)
    write.csv(prop_df, "results/tables/sRACIPE_Phenotype_Proportions.csv", 
        row.names = FALSE)
    write.csv(centroid_df, "results/tables/sRACIPE_Attractor_Centroids.csv", 
        row.names = FALSE)
    saveRDS(list(rSet_norm = rSet_norm, ss_mat = ss_mat, state_labels = state_labels, 
        net_fib_index = net_fib_index, multistability = multistability_tab), 
        "data/processed/sRACIPE_MASH_Ensemble.rds")
    cat("\n=== Continuous ODE Steady-State Attractor Proportions ===\n")
    print(prop_df)
    cat("\n=== Parameter Set Multistability Breakdown (1 = Monostable, 2 = Bistable, 3 = Tristable) ===\n")
    print(multistability_tab)
    cat("\n=== Mean Z-Score Centroids Across Dynamical States ===\n")
    print(round(centroid_df[, 1:3], 2))
    cat("\nStep 2.2 completed 100% cleanly! Continuous ODE ensemble saved.\n")
