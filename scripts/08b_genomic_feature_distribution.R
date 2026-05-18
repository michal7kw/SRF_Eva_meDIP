#!/usr/bin/env Rscript

################################################################################
# Script: 08b_genomic_feature_distribution.R
# Purpose: Enhanced genomic feature distribution plots for DMRs
#
# Description:
#   Creates publication-quality multi-panel visualizations showing the
#   distribution of DMRs across genomic features (promoter, exon, intron,
#   intergenic, UTR), with separate panels for hypermethylated vs
#   hypomethylated DMRs.
#
# Input:
#   - Annotated DMRs: ../results/08_annotation/*_annotated.csv
#
# Output:
#   - Multi-panel PDFs and PNGs (300 DPI)
#   - Summary statistics CSV
#
# Dependencies: ggplot2, dplyr, patchwork, scales, tidyr
#
# Runtime: ~10-15 minutes
################################################################################

cat("=================================================\n")
cat("DMR Genomic Feature Distribution Analysis\n")
cat("=================================================\n")
cat(paste("Start:", Sys.time(), "\n\n"))

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(scales)
    library(patchwork)
})

################################################################################
# Configuration
################################################################################

base_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/meDIP"
in_dir <- file.path(base_dir, "results/08_annotation")
out_dir <- file.path(base_dir, "results/08_annotation/feature_distribution")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

cat(paste("Input directory:", in_dir, "\n"))
cat(paste("Output directory:", out_dir, "\n\n"))

################################################################################
# Color Palettes
################################################################################

# Genomic feature colors (consistent across all panels)
feature_colors <- c(
    "Promoter" = "#E41A1C",
    "5' UTR" = "#377EB8",
    "Exon" = "#4DAF4A",
    "Intron" = "#984EA3",
    "3' UTR" = "#FF7F00",
    "Downstream" = "#FFFF33",
    "Distal Intergenic" = "#A65628",
    "Other" = "#999999"
)

feature_order <- c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR",
                   "Downstream", "Distal Intergenic", "Other")

# Direction colors
hyper_color <- "#D73027"  # Dark red for hypermethylated
hypo_color <- "#4575B4"   # Dark blue for hypomethylated
all_color <- "#756BB1"    # Purple for combined

################################################################################
# Helper Functions
################################################################################

#' Simplify ChIPseeker annotation categories
#' @param anno_df Data frame with 'annotation' column from ChIPseeker
#' @return Data frame with additional 'feature_simple' column
simplify_annotation <- function(anno_df) {
    anno_df %>%
        mutate(
            feature_simple = case_when(
                grepl("Promoter", annotation) ~ "Promoter",
                grepl("5' UTR", annotation) ~ "5' UTR",
                grepl("3' UTR", annotation) ~ "3' UTR",
                grepl("Exon", annotation) ~ "Exon",
                grepl("Intron", annotation) ~ "Intron",
                grepl("Downstream", annotation) ~ "Downstream",
                grepl("Intergenic", annotation) ~ "Distal Intergenic",
                TRUE ~ "Other"
            ),
            feature_simple = factor(feature_simple, levels = feature_order)
        )
}

#' Create pie chart showing genomic feature distribution
#' @param data Data frame with 'feature_simple' column
#' @param title Plot title
#' @param colors Named vector of colors
#' @return ggplot object
create_feature_pie <- function(data, title, colors = feature_colors) {
    # Calculate percentages
    feature_counts <- data %>%
        count(feature_simple, .drop = FALSE) %>%
        mutate(
            percentage = n / sum(n) * 100,
            label = ifelse(percentage >= 2,
                           paste0(round(percentage, 1), "%"),
                           "")
        )

    ggplot(feature_counts, aes(x = "", y = percentage, fill = feature_simple)) +
        geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
        coord_polar("y", start = 0) +
        geom_text(aes(label = label),
                  position = position_stack(vjust = 0.5),
                  size = 3, color = "white", fontface = "bold") +
        scale_fill_manual(values = colors, name = "Feature", drop = FALSE) +
        labs(title = title) +
        theme_void(base_size = 11) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            legend.position = "right",
            legend.title = element_text(face = "bold"),
            plot.margin = margin(5, 5, 5, 5)
        )
}

#' Create bar chart showing absolute feature counts
#' @param data Data frame with 'feature_simple' column
#' @param title Plot title
#' @param fill_color Single fill color for bars
#' @return ggplot object
create_feature_bar <- function(data, title, fill_color) {
    feature_counts <- data %>%
        count(feature_simple, .drop = FALSE) %>%
        mutate(feature_simple = factor(feature_simple, levels = feature_order))

    ggplot(feature_counts, aes(x = feature_simple, y = n)) +
        geom_bar(stat = "identity", fill = fill_color, color = "black", linewidth = 0.3) +
        geom_text(aes(label = comma(n)), vjust = -0.3, size = 3) +
        labs(title = title, x = NULL, y = "Number of DMRs") +
        theme_classic(base_size = 11) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
            axis.title.y = element_text(size = 10),
            plot.margin = margin(5, 10, 5, 5)
        ) +
        scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15)))
}

#' Create TSS distance histogram
#' @param data Data frame with 'distanceToTSS' column
#' @param title Plot title
#' @param fill_color Fill color for histogram
#' @return ggplot object
create_tss_distance_plot <- function(data, title, fill_color) {
    # Filter to reasonable range (within 100kb of TSS)
    plot_data <- data %>%
        filter(!is.na(distanceToTSS)) %>%
        filter(abs(distanceToTSS) <= 100000)

    if (nrow(plot_data) == 0) {
        # Return empty plot with message
        return(
            ggplot() +
                annotate("text", x = 0.5, y = 0.5, label = "No data", size = 5) +
                theme_void() +
                labs(title = title)
        )
    }

    ggplot(plot_data, aes(x = distanceToTSS / 1000)) +
        geom_histogram(bins = 50, fill = fill_color, color = "black", linewidth = 0.2) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
        labs(
            title = title,
            x = "Distance to TSS (kb)",
            y = "Number of DMRs"
        ) +
        theme_classic(base_size = 11) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            axis.title = element_text(size = 10),
            plot.margin = margin(5, 10, 5, 5)
        ) +
        scale_y_continuous(labels = comma) +
        scale_x_continuous(breaks = seq(-100, 100, 25))
}

#' Create combined multi-panel figure
#' @param all_data All DMRs data
#' @param hyper_data Hypermethylated DMRs data
#' @param hypo_data Hypomethylated DMRs data
#' @param contrast_name Name of the contrast for title
#' @return patchwork object
create_combined_figure <- function(all_data, hyper_data, hypo_data, contrast_name) {

    # Row 1: All DMRs
    n_all <- comma(nrow(all_data))
    p1_pie <- create_feature_pie(all_data, paste0("All DMRs (n=", n_all, ")"))
    p1_bar <- create_feature_bar(all_data, "All DMRs", all_color)
    p1_tss <- create_tss_distance_plot(all_data, "All DMRs", all_color)

    # Row 2: Hypermethylated
    n_hyper <- comma(nrow(hyper_data))
    p2_pie <- create_feature_pie(hyper_data, paste0("Hypermethylated (n=", n_hyper, ")"))
    p2_bar <- create_feature_bar(hyper_data, "Hypermethylated", hyper_color)
    p2_tss <- create_tss_distance_plot(hyper_data, "Hypermethylated", hyper_color)

    # Row 3: Hypomethylated
    n_hypo <- comma(nrow(hypo_data))
    p3_pie <- create_feature_pie(hypo_data, paste0("Hypomethylated (n=", n_hypo, ")"))
    p3_bar <- create_feature_bar(hypo_data, "Hypomethylated", hypo_color)
    p3_tss <- create_tss_distance_plot(hypo_data, "Hypomethylated", hypo_color)

    # Combine with patchwork (3 rows x 3 columns)
    combined <- (p1_pie | p1_bar | p1_tss) /
                (p2_pie | p2_bar | p2_tss) /
                (p3_pie | p3_bar | p3_tss) +
        plot_annotation(
            title = paste0("DMR Genomic Feature Distribution: ", contrast_name),
            subtitle = "Rows: All DMRs | Hypermethylated | Hypomethylated    |    Columns: Pie Chart | Bar Chart | TSS Distance",
            theme = theme(
                plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40")
            )
        ) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")

    return(combined)
}

#' Generate summary statistics
#' @param all_data All DMRs data
#' @param hyper_data Hypermethylated DMRs data
#' @param hypo_data Hypomethylated DMRs data
#' @param contrast_name Name of the contrast
#' @param out_dir Output directory
generate_summary_stats <- function(all_data, hyper_data, hypo_data, contrast_name, out_dir) {

    # Helper function to summarize features
    summarize_features <- function(data, label) {
        data %>%
            count(feature_simple, .drop = FALSE) %>%
            mutate(
                subset = label,
                percentage = round(n / sum(n) * 100, 2)
            ) %>%
            select(subset, feature_simple, n, percentage)
    }

    # Combine summaries
    summary_df <- bind_rows(
        summarize_features(all_data, "All_DMRs"),
        summarize_features(hyper_data, "Hypermethylated"),
        summarize_features(hypo_data, "Hypomethylated")
    )

    # Pivot to wide format for easier reading
    summary_wide <- summary_df %>%
        pivot_wider(
            id_cols = feature_simple,
            names_from = subset,
            values_from = c(n, percentage),
            names_glue = "{subset}_{.value}"
        ) %>%
        arrange(factor(feature_simple, levels = feature_order))

    # Add totals row
    totals <- data.frame(
        feature_simple = "TOTAL",
        All_DMRs_n = nrow(all_data),
        All_DMRs_percentage = 100,
        Hypermethylated_n = nrow(hyper_data),
        Hypermethylated_percentage = 100,
        Hypomethylated_n = nrow(hypo_data),
        Hypomethylated_percentage = 100
    )

    summary_wide <- bind_rows(summary_wide, totals)

    # Save
    out_file <- file.path(out_dir, paste0(contrast_name, "_feature_summary.csv"))
    write.csv(summary_wide, out_file, row.names = FALSE)
    cat(paste("  Summary saved:", basename(out_file), "\n"))

    return(summary_wide)
}

################################################################################
# Main Processing
################################################################################

# Find all annotated DMR files
anno_files <- list.files(in_dir, pattern = "_annotated\\.csv$", full.names = TRUE)

if (length(anno_files) == 0) {
    stop("ERROR: No annotated DMR files found in ", in_dir, "\n",
         "Expected files matching pattern: *_annotated.csv\n",
         "Please run 08_annotation.R first.")
}

cat(paste("Found", length(anno_files), "annotated DMR file(s)\n\n"))

for (anno_file in anno_files) {
    contrast_name <- gsub("_annotated\\.csv$", "", basename(anno_file))

    cat("==========================================\n")
    cat(paste("Processing:", contrast_name, "\n"))
    cat("==========================================\n")

    # Load annotated DMRs
    cat("Loading data...\n")
    dmr_data <- read.csv(anno_file, stringsAsFactors = FALSE)
    cat(paste("  Loaded", nrow(dmr_data), "DMRs\n"))

    # Simplify annotation categories
    dmr_data <- simplify_annotation(dmr_data)

    # Split by methylation direction
    if (!"direction" %in% colnames(dmr_data)) {
        # Add direction column if missing (based on logFC)
        dmr_data <- dmr_data %>%
            mutate(direction = ifelse(logFC > 0, "Hypermethylated", "Hypomethylated"))
    }

    hyper_data <- dmr_data %>% filter(direction == "Hypermethylated")
    hypo_data <- dmr_data %>% filter(direction == "Hypomethylated")

    cat(paste("  Hypermethylated:", nrow(hyper_data), "DMRs\n"))
    cat(paste("  Hypomethylated:", nrow(hypo_data), "DMRs\n\n"))

    # Generate summary statistics
    cat("Generating summary statistics...\n")
    summary_stats <- generate_summary_stats(dmr_data, hyper_data, hypo_data,
                                            contrast_name, out_dir)

    # Create combined multi-panel figure
    cat("Creating multi-panel figure...\n")
    combined_fig <- create_combined_figure(dmr_data, hyper_data, hypo_data, contrast_name)

    # Save combined figure
    pdf_file <- file.path(out_dir, paste0(contrast_name, "_feature_distribution_combined.pdf"))
    png_file <- file.path(out_dir, paste0(contrast_name, "_feature_distribution_combined.png"))

    ggsave(pdf_file, combined_fig, width = 18, height = 15, units = "in")
    cat(paste("  Saved:", basename(pdf_file), "\n"))

    ggsave(png_file, combined_fig, width = 18, height = 15, units = "in", dpi = 300)
    cat(paste("  Saved:", basename(png_file), "\n"))

    # Create individual subset plots (optional but useful)
    cat("\nCreating individual plots...\n")

    # Individual pie charts
    for (subset_name in c("All", "Hypermethylated", "Hypomethylated")) {
        subset_data <- switch(subset_name,
            "All" = dmr_data,
            "Hypermethylated" = hyper_data,
            "Hypomethylated" = hypo_data
        )
        fill_col <- switch(subset_name,
            "All" = all_color,
            "Hypermethylated" = hyper_color,
            "Hypomethylated" = hypo_color
        )

        # Pie chart
        pie_plot <- create_feature_pie(subset_data,
                                       paste0(contrast_name, " - ", subset_name, " DMRs"))
        pie_file <- file.path(out_dir, paste0(contrast_name, "_", tolower(subset_name), "_pie.png"))
        ggsave(pie_file, pie_plot, width = 8, height = 6, dpi = 300)

        # Bar chart
        bar_plot <- create_feature_bar(subset_data,
                                       paste0(contrast_name, " - ", subset_name, " DMRs"),
                                       fill_col)
        bar_file <- file.path(out_dir, paste0(contrast_name, "_", tolower(subset_name), "_bar.png"))
        ggsave(bar_file, bar_plot, width = 10, height = 6, dpi = 300)

        # TSS distance
        tss_plot <- create_tss_distance_plot(subset_data,
                                             paste0(contrast_name, " - ", subset_name, " DMRs"),
                                             fill_col)
        tss_file <- file.path(out_dir, paste0(contrast_name, "_", tolower(subset_name), "_tss_distance.png"))
        ggsave(tss_file, tss_plot, width = 10, height = 6, dpi = 300)
    }
    cat("  Individual plots saved.\n")

    cat("\n")
}

################################################################################
# Print Summary
################################################################################

cat("==========================================\n")
cat("Analysis Complete!\n")
cat("==========================================\n")
cat(paste("End:", Sys.time(), "\n\n"))

cat("Output files in:", out_dir, "\n\n")

cat("Key outputs:\n")
cat("1. *_feature_distribution_combined.pdf - Multi-panel figure (18x15 inches)\n")
cat("2. *_feature_distribution_combined.png - Multi-panel figure (300 DPI)\n")
cat("3. *_feature_summary.csv - Feature count statistics\n")
cat("4. *_{all,hyper,hypo}_{pie,bar,tss_distance}.png - Individual plots\n\n")

cat("Figure layout:\n")
cat("+-------------------+-------------------+-------------------+\n")
cat("|   All DMRs        |   All DMRs        |   All DMRs        |\n")
cat("|   Pie Chart       |   Bar Chart       |   TSS Distance    |\n")
cat("+-------------------+-------------------+-------------------+\n")
cat("|   Hypermethylated |   Hypermethylated |   Hypermethylated |\n")
cat("|   Pie Chart       |   Bar Chart       |   TSS Distance    |\n")
cat("+-------------------+-------------------+-------------------+\n")
cat("|   Hypomethylated  |   Hypomethylated  |   Hypomethylated  |\n")
cat("|   Pie Chart       |   Bar Chart       |   TSS Distance    |\n")
cat("+-------------------+-------------------+-------------------+\n\n")
