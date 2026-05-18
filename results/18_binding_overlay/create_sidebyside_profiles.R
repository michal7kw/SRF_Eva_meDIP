#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(gridExtra)
    library(grid)
    library(reshape2)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]

cat("Creating publication-quality side-by-side profile figures...\n")

# Publication-ready color palette (ColorBrewer)
COLORS_BINDING <- c("#2166AC", "#B2182B")  # Blue and Red
COLORS_MEDIP <- c("#A6D96A", "#1B7837", "#DFC27D", "#762A83")  # Light green, dark green, tan, purple

# Publication-ready theme
theme_publication <- function(base_size = 12) {
    theme_bw(base_size = base_size) +
    theme(
        # Title
        plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        # Axis
        axis.title = element_text(size = base_size, face = "bold"),
        axis.text = element_text(size = base_size - 1, color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        # Legend
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = base_size - 1),
        legend.key.size = unit(0.8, "lines"),
        legend.background = element_rect(fill = "white", color = NA),
        legend.margin = margin(t = 0, b = 5),
        # Panel
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        # Margins
        plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
    )
}

# Function to read deepTools matrix and create profile data
read_deeptools_matrix <- function(matrix_file) {
    if (!file.exists(matrix_file)) {
        return(NULL)
    }

    # Read the gzipped matrix
    con <- gzfile(matrix_file, "r")
    header <- readLines(con, n = 1)
    close(con)

    # Parse header for sample labels and parameters
    header_json <- gsub("^@", "", header)
    header_info <- jsonlite::fromJSON(header_json)

    # Read the data
    mat_data <- read.table(gzfile(matrix_file), skip = 1, header = FALSE)

    # Extract sample labels
    sample_labels <- header_info$sample_labels

    # Calculate number of bins
    n_bins <- (ncol(mat_data) - 6) / length(sample_labels)

    # Extract signal columns for each sample
    profiles <- list()
    for (i in seq_along(sample_labels)) {
        start_col <- 7 + (i - 1) * n_bins
        end_col <- start_col + n_bins - 1
        sample_mat <- mat_data[, start_col:end_col]
        profiles[[sample_labels[i]]] <- colMeans(sample_mat, na.rm = TRUE)
    }

    # Create position vector (assuming -5000 to +5000 bp)
    positions <- seq(-5000, 5000, length.out = n_bins)

    # Convert to data frame for ggplot
    df <- data.frame(Position = positions)
    for (name in names(profiles)) {
        df[[name]] <- profiles[[name]]
    }

    return(df)
}

# Function to create publication-quality profile plot
create_profile_plot <- function(df, samples, colors, title, ylab = "Signal") {
    if (is.null(df)) return(NULL)

    # Melt data for ggplot
    df_long <- melt(df, id.vars = "Position",
                    measure.vars = samples,
                    variable.name = "Sample",
                    value.name = "Signal")

    p <- ggplot(df_long, aes(x = Position, y = Signal, color = Sample)) +
        geom_line(linewidth = 1.2) +
        scale_color_manual(values = colors) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
        scale_x_continuous(
            breaks = c(-5000, -2500, 0, 2500, 5000),
            labels = c("-5kb", "-2.5kb", "DMR", "+2.5kb", "+5kb"),
            expand = c(0.02, 0)
        ) +
        scale_y_continuous(expand = c(0.02, 0)) +
        labs(title = title, x = NULL, y = ylab) +
        theme_publication(base_size = 11) +
        guides(color = guide_legend(nrow = 1))

    return(p)
}

# Process each DMR category
dmr_types <- c("hypermethylated", "hypomethylated", "all_dmrs", "stringent")

for (dmr_type in dmr_types) {
    binding_matrix <- file.path(outdir, "matrices", paste0(dmr_type, "_binding_matrix.gz"))
    medip_matrix <- file.path(outdir, "matrices", paste0(dmr_type, "_medip_matrix.gz"))

    if (!file.exists(binding_matrix) || !file.exists(medip_matrix)) {
        cat(paste("  Skipping", dmr_type, "- matrix files not found\n"))
        next
    }

    cat(paste("  Processing", dmr_type, "...\n"))

    tryCatch({
        # Read matrices
        binding_df <- read_deeptools_matrix(binding_matrix)
        medip_df <- read_deeptools_matrix(medip_matrix)

        if (is.null(binding_df) || is.null(medip_df)) {
            cat(paste("    Failed to read matrices for", dmr_type, "\n"))
            next
        }

        # Create binding plot
        binding_samples <- colnames(binding_df)[colnames(binding_df) != "Position"]
        p_binding <- create_profile_plot(
            binding_df,
            binding_samples,
            COLORS_BINDING,
            "TF Binding (Cut&Tag)",
            "Signal"
        )

        # Create meDIP plot
        medip_samples <- colnames(medip_df)[colnames(medip_df) != "Position"]
        p_medip <- create_profile_plot(
            medip_df,
            medip_samples,
            COLORS_MEDIP,
            "DNA Methylation (meDIP)",
            "Signal (RPKM)"
        )

        # Combine side by side - PNG output
        output_file <- file.path(outdir, "profiles", paste0(dmr_type, "_sidebyside_profile.png"))

        png(output_file, width = 14, height = 5, units = "in", res = 300)
        grid.arrange(
            p_binding, p_medip,
            ncol = 2,
            top = textGrob(
                paste0(tools::toTitleCase(gsub("_", " ", dmr_type)), " DMRs: Binding vs Methylation"),
                gp = gpar(fontsize = 14, fontface = "bold")
            )
        )
        dev.off()

        cat(paste("    Saved:", output_file, "\n"))

    }, error = function(e) {
        cat(paste("    Error processing", dmr_type, ":", e$message, "\n"))
    })
}

cat("Done creating side-by-side figures.\n")
