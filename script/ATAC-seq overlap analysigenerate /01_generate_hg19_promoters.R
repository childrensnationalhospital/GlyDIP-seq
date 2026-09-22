#!/usr/bin/env Rscript

# ==============================================================================
# Generate hg19 promoter regions (TSS +/- 2 kb)
#
# Reference annotation:
#   TxDb.Hsapiens.UCSC.hg19.knownGene
#
# Output:
#   hg19_promoters.bed
#
# BED coordinates are 0-based, half-open.
# ==============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
})

# ------------------------------------------------------------------------------
# Load hg19 gene annotation
# ------------------------------------------------------------------------------

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

genes_gr <- genes(txdb)

# ------------------------------------------------------------------------------
# Define promoter regions as TSS +/- 2 kb
# ------------------------------------------------------------------------------

promoters_gr <- promoters(
  genes_gr,
  upstream = 2000,
  downstream = 2000
)

# ------------------------------------------------------------------------------
# Convert GRanges to BED format
#
# GRanges: 1-based
# BED:     0-based
# ------------------------------------------------------------------------------

promoter_df <- data.frame(
  chr   = as.character(seqnames(promoters_gr)),
  start = start(promoters_gr) - 1,
  end   = end(promoters_gr)
)

# ------------------------------------------------------------------------------
# Retain standard chromosomes only
# ------------------------------------------------------------------------------

standard_chr <- paste0(
  "chr",
  c(1:22, "X", "Y")
)

promoter_df <- promoter_df[
  promoter_df$chr %in% standard_chr &
  promoter_df$start >= 0,
]

# ------------------------------------------------------------------------------
# Sort genomic coordinates
# ------------------------------------------------------------------------------

promoter_df <- promoter_df[
  order(
    match(promoter_df$chr, standard_chr),
    promoter_df$start,
    promoter_df$end
  ),
]

# ------------------------------------------------------------------------------
# Export BED
# ------------------------------------------------------------------------------

write.table(
  promoter_df,
  file = "hg19_promoters.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
