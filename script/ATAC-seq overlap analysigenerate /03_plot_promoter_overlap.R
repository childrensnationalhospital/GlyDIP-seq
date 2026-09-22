#!/usr/bin/env Rscript

# ==============================================================================
# Visualization and statistical comparison of promoter overlap
#
# Inputs:
#   promoter_overlap_summary.tsv
#   GlyDIP_ATAC_permutations.tsv
#   OxiDIP_ATAC_permutations.tsv
#
# Outputs:
#   Figure_ATAC_promoter_overlap.pdf
#   GlyDIP_vs_OxiDIP_statistics.tsv
# ==============================================================================


# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
})


# ==============================================================================
# Load summary data
# ==============================================================================

summary_df <- read.table(
  "promoter_overlap_summary.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

print(summary_df)


# ==============================================================================
# Load permutation distributions
# ==============================================================================

gly_perm <- read.table(
  "GlyDIP_ATAC_permutations.tsv",
  header = TRUE,
  sep = "\t"
)

oxi_perm <- read.table(
  "OxiDIP_ATAC_permutations.tsv",
  header = TRUE,
  sep = "\t"
)


# Add labels
gly_perm$dataset <- "GlyDIP-seq"
oxi_perm$dataset <- "OxiDIP-seq"

perm_df <- rbind(
  gly_perm,
  oxi_perm
)


# ==============================================================================
# Prepare observed values
# ==============================================================================

observed_df <- data.frame(
  dataset = c(
    "GlyDIP-seq",
    "OxiDIP-seq"
  ),
  overlap_percent = c(
    summary_df$observed_percent[
      summary_df$dataset == "GlyDIP"
    ],
    summary_df$observed_percent[
      summary_df$dataset == "OxiDIP"
    ]
  )
)


# ==============================================================================
# Figure:
# Null distribution + observed overlap
# ==============================================================================

perm_df$dataset <- factor(
  perm_df$dataset,
  levels = c(
    "GlyDIP-seq",
    "OxiDIP-seq"
  )
)

observed_df$dataset <- factor(
  observed_df$dataset,
  levels = c(
    "GlyDIP-seq",
    "OxiDIP-seq"
  )
)


p <- ggplot(
  perm_df,
  aes(
    x = dataset,
    y = overlap_percent
  )
) +

  # Null distribution from 1,000 permutations
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +

  # Individual permutation values
  geom_jitter(
    width = 0.12,
    alpha = 0.15,
    size = 0.8
  ) +

  # Observed values
  geom_point(
    data = observed_df,
    aes(
      x = dataset,
      y = overlap_percent
    ),
    inherit.aes = FALSE,
    size = 4,
    shape = 18
  ) +

  labs(
    x = NULL,
    y = "DIP-seq peaks overlapping ATAC-seq peaks (%)"
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    axis.text.x = element_text(
      size = 13
    ),
    axis.title.y = element_text(
      size = 13
    )
  )


ggsave(
  "Figure_ATAC_promoter_overlap.pdf",
  plot = p,
  width = 5,
  height = 5
)


# ==============================================================================
# Direct comparison:
# GlyDIP-seq vs OxiDIP-seq
#
# Construct 2 x 2 table:
#
#                    ATAC overlap    No overlap
# GlyDIP
# OxiDIP
#
# Fisher's exact test is used for comparison of overlap proportions.
# ==============================================================================


gly_total <- summary_df$total_peaks[
  summary_df$dataset == "GlyDIP"
]

gly_overlap <- summary_df$observed_count[
  summary_df$dataset == "GlyDIP"
]


oxi_total <- summary_df$total_peaks[
  summary_df$dataset == "OxiDIP"
]

oxi_overlap <- summary_df$observed_count[
  summary_df$dataset == "OxiDIP"
]


comparison_table <- matrix(
  c(
    gly_overlap,
    gly_total - gly_overlap,
    oxi_overlap,
    oxi_total - oxi_overlap
  ),
  nrow = 2,
  byrow = TRUE
)


rownames(comparison_table) <- c(
  "GlyDIP-seq",
  "OxiDIP-seq"
)

colnames(comparison_table) <- c(
  "ATAC_overlap",
  "No_ATAC_overlap"
)


print(comparison_table)


# ------------------------------------------------------------------------------
# Fisher's exact test
# ------------------------------------------------------------------------------

fisher_result <- fisher.test(
  comparison_table
)

print(fisher_result)


# ==============================================================================
# Save statistics
# ==============================================================================

statistics_df <- data.frame(
  GlyDIP_overlap_percent =
    100 * gly_overlap / gly_total,

  OxiDIP_overlap_percent =
    100 * oxi_overlap / oxi_total,

  odds_ratio =
    unname(fisher_result$estimate),

  p_value =
    fisher_result$p.value
)


write.table(
  statistics_df,
  file = "GlyDIP_vs_OxiDIP_statistics.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ==============================================================================
# Print summary
# ==============================================================================

cat("\n========================================\n")
cat("GlyDIP vs OxiDIP\n")
cat("========================================\n")

cat(
  "GlyDIP overlap:",
  round(
    100 * gly_overlap / gly_total,
    2
  ),
  "%\n"
)

cat(
  "OxiDIP overlap:",
  round(
    100 * oxi_overlap / oxi_total,
    2
  ),
  "%\n"
)

cat(
  "Odds ratio:",
  round(
    unname(fisher_result$estimate),
    3
  ),
  "\n"
)

cat(
  "Fisher's exact P:",
  format.pval(
    fisher_result$p.value,
    digits = 3
  ),
  "\n"
)
