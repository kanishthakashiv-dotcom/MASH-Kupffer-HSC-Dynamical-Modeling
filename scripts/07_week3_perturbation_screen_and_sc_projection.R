    setwd("C:/Users/KANISHTHA/Desktop/new project/MASH_Systems_Immunology")
    library(Seurat)
    library(Matrix)
    library(dplyr)
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
    simulate_perturbed <- function(init_state, pert_mult, drive_val = 0, 
        steps = 450, dt = 0.05) {
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
                prod_rate <- prod_rate * pert_mult[nd]
                dstate[nd] <- prod_rate - 1 * state[nd]
            }
            state <- pmax(0, state + dt * dstate)
        }
        return(state)
    }
    base_mult <- rep(1, n_nodes)
    names(base_mult) <- nodes
    s0 <- rep(0.2, n_nodes)
    names(s0) <- nodes
    s0[c("PPARGHSC", "IL10PPARG")] <- 1.5
    s_high <- simulate_perturbed(s0, base_mult, drive_val = 2.5)
    locked_fibrotic <- simulate_perturbed(s_high, base_mult, 
        drive_val = 0)
    fib_readout <- function(st) mean(st[c("COL1A1ACTA2", "TIMP1", 
        "YAP1TAZ")]) - mean(st[c("PPARGHSC", "IL10PPARG")])
    base_fib <- fib_readout(locked_fibrotic)
    quiescent_ref <- fib_readout(s0)
    homeostatic_nodes <- c("PPARGHSC", "IL10PPARG", "SMAD7")
    get_node_mult <- function(nd) ifelse(nd %in% homeostatic_nodes, 
        3, 0.1)
    get_node_label <- function(nd) ifelse(nd %in% homeostatic_nodes, 
        paste0(nd, "_ACT"), paste0(nd, "_KD"))
    single_res <- lapply(nodes, function(nd) {
        pm <- base_mult
        pm[nd] <- get_node_mult(nd)
        st_end <- simulate_perturbed(locked_fibrotic, pm, drive_val = 0)
        f_end <- fib_readout(st_end)
        collapse_pct <- round((base_fib - f_end)/(base_fib - 
            quiescent_ref) * 100, 2)
        data.frame(Type = "Single_Node", Perturbation = get_node_label(nd), 
            Node1 = nd, Node2 = "None", Final_COL1A1_ACTA2 = round(st_end["COL1A1ACTA2"], 
                3), Final_YAP1_TAZ = round(st_end["YAP1TAZ"], 
                3), Final_PPARG_HSC = round(st_end["PPARGHSC"], 
                3), Fibrotic_Collapse_Pct = collapse_pct, Attractor_Reverted = (st_end["COL1A1ACTA2"] < 
                0.8 & st_end["PPARGHSC"] > 1))
    })
    single_df <- bind_rows(single_res) %>% arrange(desc(Fibrotic_Collapse_Pct))
    single_lookup <- setNames(single_df$Fibrotic_Collapse_Pct, 
        single_df$Node1)
    pairs <- combn(nodes, 2, simplify = FALSE)
    dual_res <- lapply(pairs, function(pr) {
        n1 <- pr[1]
        n2 <- pr[2]
        pm <- base_mult
        pm[n1] <- get_node_mult(n1)
        pm[n2] <- get_node_mult(n2)
        st_end <- simulate_perturbed(locked_fibrotic, pm, drive_val = 0)
        f_end <- fib_readout(st_end)
        collapse_pct <- round((base_fib - f_end)/(base_fib - 
            quiescent_ref) * 100, 2)
        e1 <- pmax(0, pmin(1, single_lookup[n1]/100))
        e2 <- pmax(0, pmin(1, single_lookup[n2]/100))
        bliss_expected <- round((e1 + e2 - e1 * e2) * 100, 2)
        synergy_score <- round(collapse_pct - bliss_expected, 
            2)
        data.frame(Type = "Dual_Node", Perturbation = paste0(get_node_label(n1), 
            " + ", get_node_label(n2)), Node1 = n1, Node2 = n2, 
            Final_COL1A1_ACTA2 = round(st_end["COL1A1ACTA2"], 
                3), Final_YAP1_TAZ = round(st_end["YAP1TAZ"], 
                3), Final_PPARG_HSC = round(st_end["PPARGHSC"], 
                3), Fibrotic_Collapse_Pct = collapse_pct, Bliss_Synergy_Score = synergy_score, 
            Attractor_Reverted = (st_end["COL1A1ACTA2"] < 0.8 & 
                st_end["PPARGHSC"] > 1))
    })
    dual_df <- bind_rows(dual_res) %>% arrange(desc(Fibrotic_Collapse_Pct), 
        desc(Bliss_Synergy_Score))
    write.csv(single_df, "results/tables/InSilico_Single_Node_Perturbations.csv", 
        row.names = FALSE)
    write.csv(dual_df, "results/tables/InSilico_Dual_Node_Synergy_Screen.csv", 
        row.names = FALSE)
    cat("=== Single-Node Perturbation Screen (Locked Fibrotic State) ===\n")
    print(single_df[, c("Perturbation", "Final_COL1A1_ACTA2", 
        "Final_PPARG_HSC", "Fibrotic_Collapse_Pct", "Attractor_Reverted")])
    cat("\n=== Top 12 Synergistic Dual-Node Perturbations (Reverting Irreversible Fibrosis) ===\n")
    print(head(dual_df[, c("Perturbation", "Final_COL1A1_ACTA2", 
        "Final_PPARG_HSC", "Fibrotic_Collapse_Pct", "Bliss_Synergy_Score", 
        "Attractor_Reverted")], 12))
    seu <- readRDS("data/processed/GSE136103_seurat_processed.rds")
    expr <- tryCatch(GetAssayData(seu, layer = "data"), error = function(e) GetAssayData(seu, 
        slot = "data"))
    gene_map <- list(TLR4NFkB = c("TLR4", "NFKB1", "RELA"), TNFIL1B = c("TNF", 
        "IL1B"), TREM2SPP1 = c("TREM2", "SPP1", "CD9"), TGFB1 = c("TGFB1"), 
        IL10PPARG = c("IL10", "CD163", "MRC1"), TGFBRSMAD3 = c("TGFBR1", 
            "TGFBR2", "SMAD3"), YAP1TAZ = c("YAP1", "WWTR1", 
            "CCN2"), PDGFRBERK = c("PDGFRB", "MAPK1", "MAPK3"), 
        NFkBHSC = c("BCL2", "BIRC3", "ICAM1"), TIMP1 = c("TIMP1", 
            "MMP2", "MMP9"), COL1A1ACTA2 = c("COL1A1", "COL1A2", 
            "ACTA2"), PPARGHSC = c("PPARG", "RXRA", "PLIN2"), 
        SMAD7 = c("SMAD7", "BAMBI"))
    niche_idx <- which(seu$Cell_State %in% c("Macrophage_Kupffer", 
        "SAMac_TREM2_SPP1", "HSC_Myofibroblast"))
    sub_expr <- expr[, niche_idx, drop = FALSE]
    sc_node_mat <- sapply(names(gene_map), function(nd) {
        g <- intersect(gene_map[[nd]], rownames(sub_expr))
        if (length(g) == 1) 
            return(as.numeric(sub_expr[g, ]))
        return(Matrix::colMeans(sub_expr[g, , drop = FALSE]))
    })
    sc_node_mat <- t(sc_node_mat)
    sc_fib_score <- colMeans(sc_node_mat[c("COL1A1ACTA2", "TIMP1", 
        "YAP1TAZ", "TGFBRSMAD3", "TREM2SPP1", "TGFB1"), , drop = FALSE]) - 
        colMeans(sc_node_mat[c("PPARGHSC", "IL10PPARG"), , drop = FALSE])
    sc_proj_df <- data.frame(Cell_Barcode = colnames(sub_expr), 
        Condition = seu$Condition[niche_idx], Cell_State = seu$Cell_State[niche_idx], 
        SingleCell_Fibrotic_Attractor_Score = sc_fib_score)
    write.csv(sc_proj_df, "results/tables/GSE136103_SingleCell_Attractor_Scores.csv", 
        row.names = FALSE)
    cat("\n=== Single-Cell (GSE136103) Fibrotic Attractor Score by Cell State & Condition ===\n")
    print(sc_proj_df %>% group_by(Cell_State, Condition) %>% 
        summarise(Cells = n(), Mean_Score = round(mean(SingleCell_Fibrotic_Attractor_Score), 
            3), .groups = "drop"))
    cat("\nStep 3.2 and ALL OF WEEK 3 completed 100% cleanly!\n")
