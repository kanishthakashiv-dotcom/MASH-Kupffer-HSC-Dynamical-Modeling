    setwd("C:/Users/KANISHTHA/Desktop/new project/MASH_Systems_Immunology")
    library(dplyr)
    topo <- read.table("network/MASH_Kupffer_HSC_15Node.topo", 
        header = TRUE, stringsAsFactors = FALSE)
    nodes <- sort(unique(c(topo$Source, topo$Target)))
    n_nodes <- length(nodes)
    hill_act <- function(x, thr = 0.85, n = 4) 1/(1 + (thr/pmax(x, 
        1e-12))^n)
    hill_inh <- function(x, thr = 0.85, n = 4) 1 - 1/(1 + (thr/pmax(x, 
        1e-12))^n)
    simulate_ode <- function(init_state, pert_mult, drive_val = 0, 
        steps = 600, dt = 0.05) {
        state <- as.numeric(init_state)
        names(state) <- nodes
        for (s in seq_len(steps)) {
            dstate <- numeric(n_nodes)
            names(dstate) <- nodes
            for (nd in nodes) {
                acts <- topo$Source[topo$Target == nd & topo$Type == 
                  1]
                inhs <- topo$Source[topo$Target == nd & topo$Type == 
                  2]
                act_term <- if (length(acts) > 0) 
                  mean(sapply(acts, function(a) hill_act(state[a])))
                else 0.1
                inh_term <- if (length(inhs) > 0) 
                  prod(sapply(inhs, function(inh) hill_inh(state[inh])))
                else 1
                base_basal <- if (nd %in% c("PPARG_HSC", "IL10_PPARG")) 
                  0.9
                else 0.08
                drive_term <- if (nd == "TLR4_NFkB") 
                  drive_val
                else 0
                prod_rate <- (base_basal + 1.8 * act_term + drive_term) * 
                  inh_term * pert_mult[nd]
                dstate[nd] <- prod_rate - 1 * state[nd]
            }
            state <- pmax(0, state + dt * dstate)
            names(state) <- nodes
        }
        return(state)
    }
    base_mult <- rep(1, n_nodes)
    names(base_mult) <- nodes
    s0 <- rep(0.05, n_nodes)
    names(s0) <- nodes
    s0[c("PPARG_HSC", "IL10_PPARG")] <- 1.6
    quiescent_steady <- simulate_ode(s0, base_mult, drive_val = 0)
    curr <- quiescent_steady
    fwd_res <- lapply(seq(0, 2, by = 0.1), function(dv) {
        curr <<- simulate_ode(curr, base_mult, drive_val = dv)
        data.frame(Direction = "Forward_Injury_Escalation", TLR4_Drive = dv, 
            COL1A1_ACTA2 = round(curr["COL1A1_ACTA2"], 4), YAP1_TAZ = round(curr["YAP1_TAZ"], 
                4), TREM2_SPP1 = round(curr["TREM2_SPP1"], 4), 
            PPARG_HSC = round(curr["PPARG_HSC"], 4), row.names = NULL)
    })
    rev_res <- lapply(seq(2, 0, by = -0.1), function(dv) {
        curr <<- simulate_ode(curr, base_mult, drive_val = dv)
        data.frame(Direction = "Reverse_Trigger_Removal", TLR4_Drive = dv, 
            COL1A1_ACTA2 = round(curr["COL1A1_ACTA2"], 4), YAP1_TAZ = round(curr["YAP1_TAZ"], 
                4), TREM2_SPP1 = round(curr["TREM2_SPP1"], 4), 
            PPARG_HSC = round(curr["PPARG_HSC"], 4), row.names = NULL)
    })
    post_removal_steady <- curr
    hysteresis_df <- bind_rows(c(fwd_res, rev_res))
    write.csv(hysteresis_df, "results/tables/Hysteresis_Trigger_Removal_Trajectory.csv", 
        row.names = FALSE)
    cat("=== Calibrated Hysteresis Proof (Drive=0 Quiescent vs Drive=0 Post-Injury Removal) ===\n")
    print(hysteresis_df[hysteresis_df$TLR4_Drive %in% c(0, 0.4, 
        0.8, 1.2, 1.6, 2), ])
    fib_readout <- function(st) mean(st[c("COL1A1_ACTA2", "TIMP1", 
        "YAP1_TAZ", "TREM2_SPP1")]) - mean(st[c("PPARG_HSC", 
        "IL10_PPARG")])
    base_fib <- fib_readout(post_removal_steady)
    quiescent_ref <- fib_readout(quiescent_steady)
    homeostatic_nodes <- c("PPARG_HSC", "IL10_PPARG", "SMAD7")
    get_node_mult <- function(nd) ifelse(nd %in% homeostatic_nodes, 
        2.5, 0.1)
    get_node_label <- function(nd) ifelse(nd %in% homeostatic_nodes, 
        paste0(nd, "_ACT"), paste0(nd, "_KD"))
    single_res <- lapply(nodes, function(nd) {
        pm <- base_mult
        pm[nd] <- get_node_mult(nd)
        st_end <- simulate_ode(post_removal_steady, pm, drive_val = 0)
        f_end <- fib_readout(st_end)
        collapse_pct <- round(pmax(0, pmin(100, (base_fib - f_end)/(base_fib - 
            quiescent_ref) * 100)), 2)
        data.frame(Type = "Single_Node", Perturbation = get_node_label(nd), 
            Node1 = nd, Node2 = "None", Final_COL1A1_ACTA2 = round(st_end["COL1A1_ACTA2"], 
                3), Final_YAP1_TAZ = round(st_end["YAP1_TAZ"], 
                3), Final_PPARG_HSC = round(st_end["PPARG_HSC"], 
                3), Fibrotic_Collapse_Pct = collapse_pct, Attractor_Reverted = (st_end["PPARG_HSC"] > 
                st_end["COL1A1_ACTA2"]), row.names = NULL)
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
        st_end <- simulate_ode(post_removal_steady, pm, drive_val = 0)
        f_end <- fib_readout(st_end)
        collapse_pct <- round(pmax(0, pmin(100, (base_fib - f_end)/(base_fib - 
            quiescent_ref) * 100)), 2)
        e1 <- pmax(0, pmin(1, single_lookup[n1]/100))
        e2 <- pmax(0, pmin(1, single_lookup[n2]/100))
        bliss_expected <- round((e1 + e2 - e1 * e2) * 100, 2)
        synergy_score <- round(collapse_pct - bliss_expected, 
            2)
        data.frame(Type = "Dual_Node", Perturbation = paste0(get_node_label(n1), 
            " + ", get_node_label(n2)), Node1 = n1, Node2 = n2, 
            Final_COL1A1_ACTA2 = round(st_end["COL1A1_ACTA2"], 
                3), Final_YAP1_TAZ = round(st_end["YAP1_TAZ"], 
                3), Final_PPARG_HSC = round(st_end["PPARG_HSC"], 
                3), Fibrotic_Collapse_Pct = collapse_pct, Bliss_Synergy_Score = synergy_score, 
            Attractor_Reverted = (st_end["PPARG_HSC"] > st_end["COL1A1_ACTA2"]), 
            row.names = NULL)
    })
    dual_df <- bind_rows(dual_res) %>% arrange(desc(Fibrotic_Collapse_Pct), 
        desc(Bliss_Synergy_Score))
    write.csv(single_df, "results/tables/InSilico_Single_Node_Perturbations.csv", 
        row.names = FALSE)
    write.csv(dual_df, "results/tables/InSilico_Dual_Node_Synergy_Screen.csv", 
        row.names = FALSE)
    cat("\n=== Calibrated Single-Node Perturbation Screen (0-100% Scale) ===\n")
    print(single_df[, c("Perturbation", "Final_COL1A1_ACTA2", 
        "Final_PPARG_HSC", "Fibrotic_Collapse_Pct", "Attractor_Reverted")])
    cat("\n=== Top 12 Synergistic Dual-Node Perturbations (Reverting Irreversible Fibrosis) ===\n")
    print(head(dual_df[, c("Perturbation", "Final_COL1A1_ACTA2", 
        "Final_PPARG_HSC", "Fibrotic_Collapse_Pct", "Bliss_Synergy_Score", 
        "Attractor_Reverted")], 12))
