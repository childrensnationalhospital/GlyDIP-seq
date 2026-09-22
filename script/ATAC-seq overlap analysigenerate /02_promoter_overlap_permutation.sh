#!/usr/bin/env bash

set -euo pipefail

# ==============================================================================
# Promoter-restricted permutation analysis
#
# Purpose:
#   Test whether GlyDIP-seq and OxiDIP-seq promoter peaks preferentially
#   overlap with ATAC-seq peaks beyond that expected by chance.
#
# Promoter definition:
#   TSS +/- 2 kb (hg19)
#
# Randomization:
#   ATAC-seq promoter peaks are randomly shuffled within promoter regions.
#
# Number of permutations:
#   1,000
#
# Required:
#   BEDTools
#
# Inputs:
#   hg19_promoters.bed
#   hg19.genome
#   ATAC_seq_GMMG1_promoter.bed
#   GlyDIP_GM_MG1_promoter.bed
#   OxiDIP_GMMG1_promoter.bed
#
# Outputs:
#   GlyDIP_ATAC_permutations.tsv
#   OxiDIP_ATAC_permutations.tsv
#   promoter_overlap_summary.tsv
# ==============================================================================


# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------

PROMOTERS="hg19_promoters.bed"
GENOME="hg19.genome"

ATAC="ATAC_seq_GMMG1_promoter.bed"

GLYDIP="GlyDIP_GM_MG1_promoter.bed"
OXIDIP="OxiDIP_GMMG1_promoter.bed"

N_PERM=1000


# ==============================================================================
# Prepare promoter genomic space
# ==============================================================================

bedtools sort \
    -i "${PROMOTERS}" \
    | bedtools merge -i - \
    > hg19_promoters_merged.bed

PROMOTER_SPACE="hg19_promoters_merged.bed"


# ==============================================================================
# Function for permutation analysis
# ==============================================================================

run_permutation () {

    DIP_FILE=$1
    LABEL=$2

    OUT="${LABEL}_ATAC_permutations.tsv"

    TMPDIR="tmp_${LABEL}_ATAC_shuffle"
    mkdir -p "${TMPDIR}"

    # --------------------------------------------------------------------------
    # Number of DIP peaks
    # --------------------------------------------------------------------------

    TOTAL_DIP=$(wc -l < "${DIP_FILE}")


    # --------------------------------------------------------------------------
    # Observed overlap
    #
    # -u reports each DIP peak once if it overlaps >=1 ATAC peak.
    # --------------------------------------------------------------------------

    OBS_COUNT=$(
        bedtools intersect \
            -u \
            -a "${DIP_FILE}" \
            -b "${ATAC}" \
        | wc -l
    )

    OBS_PERCENT=$(
        awk \
            -v x="${OBS_COUNT}" \
            -v n="${TOTAL_DIP}" \
            'BEGIN {
                printf "%.4f", 100*x/n
            }'
    )


    echo ""
    echo "========================================"
    echo "${LABEL}"
    echo "========================================"

    echo "Total ${LABEL} peaks : ${TOTAL_DIP}"
    echo "Observed overlap     : ${OBS_COUNT}"
    echo "Observed percentage  : ${OBS_PERCENT}%"


    # --------------------------------------------------------------------------
    # Permutations
    # --------------------------------------------------------------------------

    echo -e \
        "iteration\toverlap_count\toverlap_percent" \
        > "${OUT}"

    for i in $(seq 1 "${N_PERM}")
    do

        SHUFFLED="${TMPDIR}/ATAC_shuffle_${i}.bed"

        # Shuffle ATAC peaks within promoter space.
        # Peak lengths are retained.
        bedtools shuffle \
            -i "${ATAC}" \
            -g "${GENOME}" \
            -incl "${PROMOTER_SPACE}" \
            -seed "${i}" \
            > "${SHUFFLED}"


        # Count DIP peaks overlapping randomized ATAC peaks.
        RANDOM_COUNT=$(
            bedtools intersect \
                -u \
                -a "${DIP_FILE}" \
                -b "${SHUFFLED}" \
            | wc -l
        )


        RANDOM_PERCENT=$(
            awk \
                -v x="${RANDOM_COUNT}" \
                -v n="${TOTAL_DIP}" \
                'BEGIN {
                    printf "%.4f", 100*x/n
                }'
        )


        echo -e \
            "${i}\t${RANDOM_COUNT}\t${RANDOM_PERCENT}" \
            >> "${OUT}"


        if (( i % 100 == 0 ))
        then
            echo "${LABEL}: ${i}/${N_PERM} permutations completed"
        fi

    done


    # --------------------------------------------------------------------------
    # Null distribution
    # --------------------------------------------------------------------------

    RANDOM_MEAN=$(
        awk \
            'NR>1 {
                sum += $3
             }
             END {
                printf "%.4f", sum/(NR-1)
             }' \
            "${OUT}"
    )


    RANDOM_SD=$(
        awk \
            -v mean="${RANDOM_MEAN}" \
            'NR>1 {
                ss += ($3-mean)^2
             }
             END {
                printf "%.4f", sqrt(ss/(NR-2))
             }' \
            "${OUT}"
    )


    # --------------------------------------------------------------------------
    # Empirical P value
    #
    # P = (k + 1) / (N + 1)
    #
    # k = permutations with overlap >= observed overlap
    # --------------------------------------------------------------------------

    K=$(
        awk \
            -v obs="${OBS_PERCENT}" \
            'NR>1 && $3 >= obs {
                n++
             }
             END {
                print n+0
             }' \
            "${OUT}"
    )


    EMPIRICAL_P=$(
        awk \
            -v k="${K}" \
            -v n="${N_PERM}" \
            'BEGIN {
                printf "%.6f", (k+1)/(n+1)
             }'
    )


    # --------------------------------------------------------------------------
    # Save summary
    # --------------------------------------------------------------------------

    echo -e \
"${LABEL}\t${TOTAL_DIP}\t${OBS_COUNT}\t${OBS_PERCENT}\t${RANDOM_MEAN}\t${RANDOM_SD}\t${K}\t${EMPIRICAL_P}" \
    >> promoter_overlap_summary.tsv


    echo ""
    echo "Observed overlap : ${OBS_PERCENT}%"
    echo "Random mean      : ${RANDOM_MEAN}%"
    echo "Random SD        : ${RANDOM_SD}%"
    echo "Random >= obs    : ${K}/${N_PERM}"
    echo "Empirical P      : ${EMPIRICAL_P}"

    rm -rf "${TMPDIR}"
}


# ==============================================================================
# Initialize summary table
# ==============================================================================

echo -e \
"dataset\ttotal_peaks\tobserved_count\tobserved_percent\trandom_mean\trandom_sd\tn_random_ge_observed\tempirical_p" \
> promoter_overlap_summary.tsv


# ==============================================================================
# GlyDIP-seq
# ==============================================================================

run_permutation \
    "${GLYDIP}" \
    "GlyDIP"


# ==============================================================================
# OxiDIP-seq
# ==============================================================================

run_permutation \
    "${OXIDIP}" \
    "OxiDIP"


echo ""
echo "========================================"
echo "All analyses completed"
echo "========================================"

cat promoter_overlap_summary.tsv
