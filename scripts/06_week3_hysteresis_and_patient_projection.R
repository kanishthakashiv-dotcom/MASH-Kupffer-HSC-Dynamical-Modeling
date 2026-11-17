    setwd("C:/Users/KANISHTHA/Desktop/new project/MASH_Systems_Immunology")
    library(dplyr)
    library(tidyr)
    ens <- readRDS("data/processed/sRACIPE_MASH_Ensemble.rds")
    topo <- read.table("network/MASH_Kupffer_HSC_15Node.topo", 
        header = TRUE, stringsAsFactors = FALSE)
    centroids <- read.csv("results/tables/sRACIPE_Attractor_Centroids.csv", 
        row.names = "Node")
    clean_name <- function(x) gsub("_", "", x)
    topo$Src_Clean <- clean_name(topo$Source)
    topo$Tgt_Clean <- clean_name(topo$Target)
    nodes <- rownames(centroids)
    n_nodes <- length(nodes)
    hill_act <- function(x, thr = 0.5, n = 4, lambda = 2) 1 + 
        (lambda - 1) * (x^n)/(thr^n + x^n)
    hill_inh <- function(x, thr = 0.5, n = 4, lambda = 0.1) lambda + 
        (1 - lambda) * (thr^n)/(thr^n + x^n)
    simulate_steady_state <- function(drive_val, init_state, 
        steps = 400, dt = 0.05) {
        state <- init_state
        names(state) <- nodes
        for (s in seq_len(steps)) {
            dstate <- numeric(n_nodes)
            names(dstate) <- nodes
            for (nd in nodes) {
                acts <- topo$Src_Clean[topo$Tgt_Clean == nd & 
                  topo$Type == 1]
                inhs <- topo$Src_Clean[topo$Tgt_Clean == nd & 
                  topo$Type == 2]
                prod_rate <- 1
                if (nd == "TLR4NFkB") 
                  prod_rate <- prod_rate + drive_val
                for (a in acts) prod_rate <- prod_rate * hill_act(state[a])
                for (inh in inhs) prod_rate <- prod_rate * hill_inh(state[inh])
                dstate[nd] <- prod_rate - 1 * state[nd]
            }
            state <- pmax(0, state + dt * dstate)
        }
        return(state)
    }
    drive_forward <- seq(0, 2.5, by = 0.1)
    drive_reverse <- seq(2.5, 0, by = -0.1)
    curr <- rep(0.2, n_nodes)
    names(curr) <- nodes
    curr[c("PPARGHSC", "IL10PPARG")] <- 1.5
    fwd_res <- lapply(drive_forward, function(dv) {
        curr <<- simulate_steady_state(dv, curr)
        data.frame(Direction = "1_Forward_Injury_Escalation", 
            TLR4_Drive = dv, COL1A1_ACTA2 = curr["COL1A1ACTA2"], 
            YAP1_TAZ = curr["YAP1TAZ"], TREM2_SPP1 = curr["TREM2SPP1"], 
            PPARG_HSC = curr["PPARGHSC"])
    })
    rev_res <- lapply(drive_reverse, function(dv) {
        curr <<- simulate_steady_state(dv, curr)
        data.frame(Direction = "2_Reverse_Trigger_Removal", TLR4_Drive = dv, 
            COL1A1_ACTA2 = curr["COL1A1ACTA2"], YAP1_TAZ = curr["YAP1TAZ"], 
            TREM2_SPP1 = curr["TREM2SPP1"], PPARG_HSC = curr["PPARGHSC"])
    })
    hysteresis_df <- bind_rows(c(fwd_res, rev_res))
    write.csv(hysteresis_df, "results/tables/Hysteresis_Trigger_Removal_Trajectory.csv", 
        row.names = FALSE)
    cat("=== Hysteresis Proof: Baseline (Drive=0) vs Post-Injury Trigger Removal (Drive=0) ===\n")
    print(hysteresis_df[hysteresis_df$TLR4_Drive %in% c(0, 1, 
        2, 2.5), ])
    gene_map <- list(TLR4NFkB = c("TLR4", "NFKB1", "RELA"), TNFIL1B = c("TNF", 
        "IL1B"), TREM2SPP1 = c("TREM2", "SPP1", "CD9"), TGFB1 = c("TGFB1"), 
        IL10PPARG = c("IL10", "CD163", "MRC1"), TGFBRSMAD3 = c("TGFBR1", 
            "TGFBR2", "SMAD3"), YAP1TAZ = c("YAP1", "WWTR1", 
            "CCN2"), PDGFRBERK = c("PDGFRB", "MAPK1", "MAPK3"), 
        NFkBHSC = c("BCL2", "BIRC3", "ICAM1"), TIMP1 = c("TIMP1", 
            "MMP2", "MMP9"), COL1A1ACTA2 = c("COL1A1", "COL1A2", 
            "ACTA2"), PPARGHSC = c("PPARG", "RXRA", "PLIN2"), 
        SMAD7 = c("SMAD7", "BAMBI"))
    patient_proj_list <- list()
    for (id in c("GSE135251", "GSE126848", "GSE162694")) {
        obj <- readRDS(file.path("data/processed", paste0(id, 
            "_processed.rds")))
        vst_z <- t(scale(t(obj$vst)))
        node_mat <- sapply(names(gene_map), function(nd) {
            g <- intersect(gene_map[[nd]], rownames(vst_z))
            if (length(g) == 1) 
                return(as.numeric(vst_z[g, ]))
            return(colMeans(vst_z[g, , drop = FALSE], na.rm = TRUE))
        })
        node_mat <- t(node_mat)
        cor_mat <- cor(node_mat, centroids[rownames(node_mat), 
            1:3], use = "pairwise.complete.obs")
        assigned_state <- colnames(cor_mat)[max.col(cor_mat, 
            ties.method = "first")]
        fib_idx <- colMeans(node_mat[c("COL1A1ACTA2", "TIMP1", 
            "YAP1TAZ", "TGFBRSMAD3", "TREM2SPP1", "TGFB1"), , 
            drop = FALSE]) - colMeans(node_mat[c("PPARGHSC", 
            "IL10PPARG"), , drop = FALSE])
        df_p <- data.frame(Cohort = id, Sample = colnames(obj$vst), 
            Stage = as.character(obj$meta$Harmonized_Stage), 
            Net_Fibrotic_Attractor_Score = fib_idx, Assigned_Dynamical_Attractor = assigned_state, 
            Cor_State1_Quiescent = cor_mat[, "State1_Quiescent_Homeostatic"], 
            Cor_State2_Inflammatory = cor_mat[, "State2_Reversible_Inflammatory"], 
            Cor_State3_Fibrotic = cor_mat[, "State3_Irreversible_Fibrotic"])
        patient_proj_list[[id]] <- df_p
    }
    all_patients <- bind_rows(patient_proj_list)
    write.csv(all_patients, "results/tables/All_416_Patients_sRACIPE_Attractor_Projection.csv", 
        row.names = FALSE)
    cat("\n=== Patient-Level Attractor Assignment Across Clinical Stages (All 3 Bulk Cohorts) ===\n")
    print(table(all_patients$Stage, all_patients$Assigned_Dynamical_Attractor))
    cat("\n=== Mean Net Fibrotic Attractor Score by Stage & Cohort ===\n")
    stage_summary <- all_patients %>% group_by(Cohort, Stage) %>% 
        summarise(N = n(), Mean_Attractor_Score = round(mean(Net_Fibrotic_Attractor_Score), 
            3), SD = round(sd(Net_Fibrotic_Attractor_Score), 
            3), .groups = "drop")
    print(as.data.frame(stage_summary))
    for (id in c("GSE135251", "GSE126848", "GSE162694")) {
        sub_p <- all_patients[all_patients$Cohort == id, ]
        kw <- kruskal.test(Net_Fibrotic_Attractor_Score ~ Stage, 
            data = sub_p)
        cat(sprintf("Kruskal-Wallis p-value for %s: %.3e\n", 
            id, kw$p.value))
    }
    cat("\nStep 3.1 completed 100% cleanly! Hysteresis & patient attractor projections saved.\n")
