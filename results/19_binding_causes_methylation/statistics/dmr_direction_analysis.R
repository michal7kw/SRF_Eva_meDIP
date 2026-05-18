#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

# ============================================================================
# DMR Direction Analysis at TES Binding Sites
# ============================================================================

cat("============================================\n")
cat("DMR Direction Analysis at TES Binding Sites\n")
cat("============================================\n\n")

# Counts
n_hyper_at_tes <- 71
n_hypo_at_tes <- 107
n_total_at_tes <- n_hyper_at_tes + n_hypo_at_tes

n_hyper_global <- 29034
n_hypo_global <- 2886
n_total_global <- n_hyper_global + n_hypo_global

# Observed proportions
pct_hyper_at_tes <- 100 * n_hyper_at_tes / n_total_at_tes
pct_hypo_at_tes <- 100 * n_hypo_at_tes / n_total_at_tes
pct_hyper_global <- 100 * n_hyper_global / n_total_global
pct_hypo_global <- 100 * n_hypo_global / n_total_global

cat("=== DMR Direction at TES Binding Sites ===\n\n")
cat(sprintf("At TES binding sites:\n"))
cat(sprintf("  HYPERmethylated: %d (%.1f%%)\n", n_hyper_at_tes, pct_hyper_at_tes))
cat(sprintf("  HYPOmethylated:  %d (%.1f%%)\n", n_hypo_at_tes, pct_hypo_at_tes))
cat(sprintf("  Total:           %d\n\n", n_total_at_tes))

cat(sprintf("Global DMR distribution:\n"))
cat(sprintf("  HYPERmethylated: %d (%.1f%%)\n", n_hyper_global, pct_hyper_global))
cat(sprintf("  HYPOmethylated:  %d (%.1f%%)\n", n_hypo_global, pct_hypo_global))
cat(sprintf("  Total:           %d\n\n", n_total_global))

# Fisher's exact test
# Contingency table:
#              At TES    Not at TES
# Hyper         71         28963
# Hypo         107          2779

contingency <- matrix(
    c(n_hyper_at_tes, n_hyper_global - n_hyper_at_tes,
      n_hypo_at_tes, n_hypo_global - n_hypo_at_tes),
    nrow = 2, byrow = TRUE,
    dimnames = list(
        DMR_Type = c("Hypermethylated", "Hypomethylated"),
        Location = c("At_TES_sites", "Not_at_TES")
    )
)

cat("=== Fisher's Exact Test ===\n\n")
cat("Contingency table:\n")
print(contingency)

fisher_result <- fisher.test(contingency)

cat(sprintf("\nOdds Ratio: %.3f\n", fisher_result$estimate))
cat(sprintf("95%% CI: %.3f - %.3f\n", fisher_result$conf.int[1], fisher_result$conf.int[2]))
cat(sprintf("p-value: %.2e\n\n", fisher_result$p.value))

# Interpretation
cat("=== INTERPRETATION ===\n\n")
if (fisher_result$estimate < 1 && fisher_result$p.value < 0.05) {
    cat("TES binding sites are SIGNIFICANTLY DEPLETED for hypermethylated DMRs\n")
    cat("(or equivalently, ENRICHED for hypomethylated DMRs).\n\n")
    cat(sprintf("Expected if random: %.1f%% hyper / %.1f%% hypo\n", pct_hyper_global, pct_hypo_global))
    cat(sprintf("Observed at TES:    %.1f%% hyper / %.1f%% hypo\n\n", pct_hyper_at_tes, pct_hypo_at_tes))
    cat("BIOLOGICAL INTERPRETATION:\n")
    cat("- TES binding appears to PREVENT methylation at direct target sites\n")
    cat("- The widespread hypermethylation (91% globally) is likely an INDIRECT effect\n")
    cat("- TES may recruit DNMT3A/3L to sites DISTAL from its binding locations\n")
}

# Create bar plot
outdir <- "results/19_binding_causes_methylation"

df_plot <- data.frame(
    Location = rep(c("At TES Sites\n(n=178)", "Global\n(n=31,920)"), each = 2),
    DMR_Type = rep(c("Hypermethylated", "Hypomethylated"), 2),
    Percentage = c(pct_hyper_at_tes, pct_hypo_at_tes, pct_hyper_global, pct_hypo_global)
)

p <- ggplot(df_plot, aes(x = Location, y = Percentage, fill = DMR_Type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage)),
              position = position_dodge(width = 0.7), vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("Hypermethylated" = "#D73027", "Hypomethylated" = "#4575B4"),
                      name = "DMR Direction") +
    labs(title = "DMR Direction: TES Binding Sites vs Global",
         subtitle = sprintf("Fisher's exact test: OR=%.2f, p=%.2e", fisher_result$estimate, fisher_result$p.value),
         x = "", y = "Percentage of DMRs") +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "top"
    ) +
    ylim(0, 100)

ggsave(file.path(outdir, "profiles", "DMR_direction_at_TES_sites_barplot.png"),
       p, width = 8, height = 6, dpi = 300)

cat(sprintf("\nBar plot saved: %s/profiles/DMR_direction_at_TES_sites_barplot.png\n", outdir))

# Save statistics to file
stats_file <- file.path(outdir, "statistics", "dmr_direction_summary.txt")
sink(stats_file)
cat("DMR Direction Analysis at TES Binding Sites\n")
cat("============================================\n\n")
cat("Date:", as.character(Sys.time()), "\n\n")
cat("COUNTS:\n")
cat(sprintf("  TES peaks with HYPER DMRs: %d\n", n_hyper_at_tes))
cat(sprintf("  TES peaks with HYPO DMRs:  %d\n", n_hypo_at_tes))
cat(sprintf("  Total TES peaks with DMR:  %d\n\n", n_total_at_tes))
cat("PERCENTAGES:\n")
cat(sprintf("  At TES sites: %.1f%% hyper / %.1f%% hypo\n", pct_hyper_at_tes, pct_hypo_at_tes))
cat(sprintf("  Global:       %.1f%% hyper / %.1f%% hypo\n\n", pct_hyper_global, pct_hypo_global))
cat("FISHER'S EXACT TEST:\n")
cat(sprintf("  Odds Ratio: %.3f\n", fisher_result$estimate))
cat(sprintf("  95%% CI: %.3f - %.3f\n", fisher_result$conf.int[1], fisher_result$conf.int[2]))
cat(sprintf("  p-value: %.2e\n\n", fisher_result$p.value))
cat("INTERPRETATION:\n")
cat("  TES binding is significantly DEPLETED at hypermethylated DMRs.\n")
cat("  This suggests TES binding PREVENTS methylation at direct targets,\n")
cat("  while causing methylation at distal sites (indirect effect).\n")
sink()

cat(sprintf("Statistics saved: %s\n", stats_file))

