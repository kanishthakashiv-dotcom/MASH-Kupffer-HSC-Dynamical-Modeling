    setwd("C:/Users/KANISHTHA/Desktop/new project/MASH_Systems_Immunology")
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    topo <- read.table("network/MASH_Kupffer_HSC_15Node.topo", 
        header = TRUE, stringsAsFactors = FALSE)
    nodes <- sort(unique(c(topo$Source, topo$Target)))
    n_nodes <- length(nodes)
    hill_act <- function(x, thr = 0.65, n = 4) 1/(1 + (thr/pmax(x, 
        1e-12))^n)
    hill_inh <- function(x, thr = 0.65, n = 4) 1 - 1/(1 + (thr/pmax(x, 
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
                else 0.08
                inh_term <- if (length(inhs) > 0) 
                  prod(sapply(inhs, function(inh) hill_inh(state[inh])))
                else 1
                base_basal <- if (nd %in% c("PPARG_HSC", "IL10_PPARG")) 
                  0.85
                else 0.06
                drive_term <- if (nd == "TLR4_NFkB") 
                  drive_val
                else 0
                prod_rate <- (base_basal + 1.85 * act_term + 
                  drive_term) * inh_term * pert_mult[nd]
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
    s0[c("PPARG_HSC", "IL10_PPARG")] <- 1.65
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
    fib_readout <- function(st) mean(st[c("COL1A1_ACTA2", "TIMP1", 
        "YAP1_TAZ", "TREM2_SPP1")]) - mean(st[c("PPARG_HSC", 
        "IL10_PPARG")])
    base_fib <- fib_readout(post_removal_steady)
    quiescent_ref <- fib_readout(quiescent_steady)
    homeostatic_nodes <- c("PPARG_HSC", "IL10_PPARG", "SMAD7")
    get_node_mult <- function(nd) ifelse(nd %in% homeostatic_nodes, 
        2.2, 0.25)
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
                st_end["COL1A1_ACTA2"] & st_end["COL1A1_ACTA2"] < 
                0.5), row.names = NULL)
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
            Attractor_Reverted = (st_end["PPARG_HSC"] > st_end["COL1A1_ACTA2"] & 
                st_end["COL1A1_ACTA2"] < 0.5), row.names = NULL)
    })
    dual_df <- bind_rows(dual_res) %>% arrange(desc(Fibrotic_Collapse_Pct), 
        desc(Bliss_Synergy_Score))
    write.csv(single_df, "results/tables/InSilico_Single_Node_Perturbations.csv", 
        row.names = FALSE)
    write.csv(dual_df, "results/tables/InSilico_Dual_Node_Synergy_Screen.csv", 
        row.names = FALSE)
    cat("=== Single-Node Summary (Reverted =", sum(single_df$Attractor_Reverted), 
        "of 13) ===\n")
    print(head(single_df[, c("Perturbation", "Final_COL1A1_ACTA2", 
        "Final_PPARG_HSC", "Fibrotic_Collapse_Pct", "Attractor_Reverted")], 
        6))
    cat("\n=== Top Dual-Node Summary (Reverted =", sum(dual_df$Attractor_Reverted), 
        "of 78) ===\n")
    print(head(dual_df[, c("Perturbation", "Final_COL1A1_ACTA2", 
        "Final_PPARG_HSC", "Fibrotic_Collapse_Pct", "Bliss_Synergy_Score", 
        "Attractor_Reverted")], 10))
    theme_npj <- function(base_size = 10) {
        theme_classic(base_size = base_size) + theme(plot.title = element_text(face = "bold", 
            size = base_size + 0.5, hjust = 0), plot.subtitle = element_text(size = base_size - 
            1.2, color = "grey30"), axis.title = element_text(face = "bold", 
            size = base_size), axis.text = element_text(color = "black", 
            size = base_size - 1), legend.title = element_text(face = "bold", 
            size = base_size - 1), legend.text = element_text(size = base_size - 
            1), strip.background = element_rect(fill = "grey94", 
            color = "black", linewidth = 0.5), strip.text = element_text(face = "bold", 
            size = base_size - 1), plot.margin = margin(12, 22, 
            12, 12))
    }
    hyst_long <- hysteresis_df %>% pivot_longer(cols = c(COL1A1_ACTA2, 
        YAP1_TAZ, TREM2_SPP1, PPARG_HSC), names_to = "Module", 
        values_to = "Activity")
    p4 <- ggplot(hyst_long, aes(x = TLR4_Drive, y = Activity, 
        color = Direction, linetype = Direction)) + geom_line(linewidth = 1.15) + 
        geom_point(size = 1.9) + facet_wrap(~Module, scales = "free_y", 
        ncol = 2) + scale_color_manual(values = c(Forward_Injury_Escalation = "#2166AC", 
        Reverse_Trigger_Removal = "#B2182B"), labels = c("Forward Injury Escalation (0 -> 2.0)", 
        "Reverse Trigger Withdrawal (2.0 -> 0)")) + scale_linetype_manual(values = c(Forward_Injury_Escalation = "solid", 
        Reverse_Trigger_Removal = "dashed"), labels = c("Forward Injury Escalation (0 -> 2.0)", 
        "Reverse Trigger Withdrawal (2.0 -> 0)")) + labs(title = "Irreversible Hysteresis Across the Kupffer/SAMac - HSC Regulatory Circuit", 
        subtitle = "Complete withdrawal of upstream TLR4/NFkB drive fails to restore quiescent PPARG_HSC or collapse COL1A1_ACTA2", 
        x = "Upstream Metabolic / Inflammatory Drive (TLR4_NFkB Input)", 
        y = "Steady-State Module Activity (a.u.)") + theme_npj() + 
        theme(legend.position = "bottom")
    ggsave("results/figures/Fig4_Hysteresis_Bifurcation_Proof.png", 
        p4, width = 10.5, height = 6.5, dpi = 300)
    ggsave("results/figures/Fig4_Hysteresis_Bifurcation_Proof.pdf", 
        p4, width = 10.5, height = 6.5)
    top_single <- single_df %>% head(8) %>% mutate(Category = "Single-Node Monotherapy")
    top_dual <- dual_df %>% head(10) %>% mutate(Category = "Dual-Node Combination")
    comb_pert <- bind_rows(top_single %>% dplyr::select(Perturbation, 
        Category, Final_COL1A1_ACTA2, Final_PPARG_HSC, Fibrotic_Collapse_Pct), 
        top_dual %>% dplyr::select(Perturbation, Category, Final_COL1A1_ACTA2, 
            Final_PPARG_HSC, Fibrotic_Collapse_Pct))
    p6a <- ggplot(comb_pert, aes(x = reorder(Perturbation, Fibrotic_Collapse_Pct), 
        y = Fibrotic_Collapse_Pct, fill = Category)) + geom_col(width = 0.72) + 
        geom_text(aes(label = sprintf("%.1f%%", Fibrotic_Collapse_Pct)), 
            hjust = -0.15, size = 2.8, fontface = "bold") + coord_flip(ylim = c(0, 
        115)) + scale_fill_manual(values = c(`Single-Node Monotherapy` = "#7F7F7F", 
        `Dual-Node Combination` = "#B2182B")) + labs(title = "A. In-Silico Reversal of Locked Fibrotic Attractor", 
        subtitle = "Top single-node monotherapies vs. top dual-node combinations", 
        x = NULL, y = "Fibrotic Attractor Collapse (%)") + theme_npj(base_size = 9.5) + 
        theme(legend.position = "bottom")
    d_pert_plot <- dual_df %>% mutate(Status = factor(ifelse(Attractor_Reverted, 
        "Reverted to Quiescent", "Locked Fibrotic"), levels = c("Locked Fibrotic", 
        "Reverted to Quiescent")))
    top_labels <- d_pert_plot %>% filter(Attractor_Reverted) %>% 
        head(5)
    p6b <- ggplot(d_pert_plot, aes(x = Final_COL1A1_ACTA2, y = Final_PPARG_HSC, 
        color = Status)) + geom_point(size = 2.8, alpha = 0.85) + 
        geom_text(data = top_labels, aes(label = Perturbation), 
            size = 2.5, vjust = -0.8, hjust = 0.5, check_overlap = TRUE, 
            show.legend = FALSE) + scale_color_manual(values = c(`Locked Fibrotic` = "#7F7F7F", 
        `Reverted to Quiescent` = "#1B7837"), name = "Final Attractor") + 
        scale_x_continuous(expand = expansion(mult = c(0.08, 
            0.15))) + scale_y_continuous(expand = expansion(mult = c(0.08, 
        0.15))) + labs(title = "B. Phase-Plane Attractor Reversion (All 78 Pairs)", 
        subtitle = "Simultaneous COL1A1/ECM collapse & PPARG_HSC restoration", 
        x = "Post-Treatment COL1A1_ACTA2 Activity", y = "Post-Treatment PPARG_HSC Activity") + 
        theme_npj(base_size = 9.5) + theme(legend.position = "bottom")
    fig6 <- (p6a | p6b) + plot_layout(widths = c(1.3, 1.05)) + 
        plot_annotation(theme = theme(plot.margin = margin(8, 
            18, 8, 8)))
    ggsave("results/figures/Fig6_Perturbation_Synergy_Screen.png", 
        fig6, width = 14, height = 6.4, dpi = 300)
    ggsave("results/figures/Fig6_Perturbation_Synergy_Screen.pdf", 
        fig6, width = 14, height = 6.4)
    cat("\nExported unclipped Fig4 and Fig6 (300 DPI PNG + Vector PDF) cleanly!\n")
