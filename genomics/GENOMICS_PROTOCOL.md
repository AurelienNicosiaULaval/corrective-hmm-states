# Genomic application: exploratory protocol

Recorded before inspecting processed profiles or fitting HMMs. This is a local analysis protocol, not an externally registered study. Earlier animal and economic pilots are retained in ETAT_SCIENTIFIQUE.md.

## Scientific question and data

Can a model distinguish sustained chromosomal allelic imbalance from changes in the allele measured at successive SNP markers? The target is regional allelic balance or imbalance, not clinical diagnosis or an estimate of the true number of cancer cell populations. A component corresponding to a low B-allele fraction and its high-fraction counterpart can describe the same imbalance. The identity of the allele designated B need not remain on the same parental chromosome along the sequence.

Use the first two accession-ordered HCC1143 tumor-normal pairs from GSE13372: GSM337641/GSM337662 and GSM337642/GSM337663. Four original CEL files are preserved. AS-CRMAv2 performs single-array preprocessing; TumorBoost will correct tumor allele fractions using the matched normal. Retain raw and corrected fractions. Coordinates for comparisons to Chiang et al. (2009) must be hg18. No anonymous, scrambled, simulated or merely genotype-called demonstration file is substituted for these measurements.

## Exploratory and validation separation

Chromosome 1 of pair 1 is the initial development series, chosen by chromosome number before profile inspection. Pair 1 chromosomes 2–11 are development extensions; 12–17 are validation and 18–22 remain confirmation until the specification is frozen. Pair 2 is a technical replication, not a second independent biological sample. No claim of patient-level generalization is supported by these two pairs.

Retain autosomal SNPs with finite measurements, unique positive genomic positions, and matched-normal B-allele fraction between 0.3 and 0.7. This filter identifies candidate heterozygous loci without inspecting the tumor profile. Thresholds 0.25–0.75 and 0.35–0.65 are sensitivity analyses. Do not clip corrected tumor values at zero or one. Reset the hidden chain at chromosome boundaries and gaps larger than one megabase. Record exclusions and actual distances. Check subsampling for effects of local linkage and marker spacing.

## First pilot

Fit univariate Gaussian-emission HMMs to the corrected tumor B-allele fractions on development chromosome 1. Compare G2, G3, G5 and mixtures M21 (two states with two and one components) and M221 (three states with two, two and one components). The mixture labels represent balance and degrees of imbalance only if their fitted components support that reading. Do not assign a copy number solely from an allele fraction. Use multiple starts, refine the best solutions, enforce the exact M21-to-G3 and M221-to-G5 inclusions, and report any likelihood nesting failure. The primary standard-deviation floor is 0.01 on the fraction scale; use 0.005 and 0.02 for sensitivity if the case progresses.

Quantify the fitted emissions, transition rows, occupancy, decoded run counts, convergence, competing maxima and boundary hits. The first pilot is for feasibility, not confirmatory performance. No model is chosen because a lower-likelihood local solution has appealing state labels.

## Evidence required before manuscript adoption

1. Biological interpretation supported by component locations and by the spatial profiles of raw measurements.
2. Comparisons at the same total number of components, with both density fit and scientifically relevant segmentation assessed.
3. Validation on chromosomes excluded from parameter estimation; technical replication reported separately.
4. Comparisons with established symmetry-aware processing of B-allele fractions and with a standard genomic segmentation method. A naive Gaussian baseline alone is insufficient.
5. A measured stability analysis under allele-label exchange, variance bounds, marker thinning and normal heterozygosity thresholds.
6. A comparison to independent sequencing-derived copy-ratio segments after coordinate verification. These are estimates from another assay, not error-free truth and not direct labels of allelic balance.
7. Explicit reporting of any persistent differences between the low- and high-fraction component transition rows. Extra states cannot simply be called spurious when the data support distinct dynamics.

The application may be rejected. Core manuscript theory, notation, title and EM exposition are not altered to accommodate a favorable empirical outcome.
