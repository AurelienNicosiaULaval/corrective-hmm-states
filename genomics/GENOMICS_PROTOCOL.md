# Genomic analysis protocol

The analysis examines whether complementary low and high B-allele fractions
can share regional transition dynamics. The fitted states describe allelic
balance and degrees of imbalance. Allele fractions alone do not identify
absolute copy number or the number of cell populations.

This guide describes the study design and candidate analyses. The
[original protocol](audit/protocol_before_confirmation.md), recorded before
profile inspection and HMM fitting, is preserved without changes. Its SHA-256
matches the protocol entry in the [confirmation freeze](audit/confirmation_freeze.json).
The [provenance record](audit/protocol_provenance.json) documents this guide's
editorial update. The study was exploratory and was not externally registered.

## Data and preprocessing

The data are HCC1143 tumor and HCC1143BL matched-normal SNP arrays from
[GEO GSE13372](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13372),
associated with [Chiang et al. (2009)](https://doi.org/10.1038/nmeth.1276).
The two pairs, selected in accession order, are GSM337641/GSM337662 and
GSM337642/GSM337663. The second pair is a technical repeat of the same cell lines.

AS-CRMAv2 processes the original CEL files, and TumorBoost corrects tumor
fractions using the matched normal. Analysis coordinates use hg18. Retained
markers have finite measurements, unique positive genomic positions and
matched-normal fractions in [0.3, 0.7]. Alternative filters use [0.25, 0.75]
and [0.35, 0.65]. Corrected tumor fractions are not clipped to [0, 1]. Each
chromosome boundary and gap exceeding one megabase starts a new sequence.
Marker thinning assesses sensitivity to spacing and local dependence.

## Development and evaluation

Chromosome 1 of the first pair was chosen for development by chromosome
number before profile inspection. Chromosomes 2–11 were reserved for possible
extensions and were not used. Chromosomes 12–17 provide validation data.
Chromosomes 18–22 were evaluated after the fitted models and specification
were frozen. They form a later held-out block within the exploratory design.
Chromosome 1 of the second pair assesses technical repeatability. None of
these comparisons supplies independent patient-level validation.

## Candidate models and estimation

Gaussian HMMs G2, G3 and G5 have two, three and five states. Mixture models
M21 and M221 have component counts (2, 1) and (2, 2, 1). The comparisons
M21–G3 and M221–G5 hold the total component count fixed. Multiple starts,
refinement of the best solutions and starts at exact nested representations
check optimization. The primary standard-deviation floor is 0.01; sensitivity
analyses use 0.005 and 0.02.

The recorded outputs include emission parameters, transition rows, posterior
occupation, decoded runs, convergence diagnostics, competing solutions and
boundary estimates. Selection within each candidate family uses likelihood.
The detailed fitting settings and symmetry-aware comparators are described
in Supporting Information S3.2.

## Interpretation and sensitivity

The assessment combines observed fraction profiles, emission components,
comparisons at equal component count, fixed-parameter prediction on excluded
chromosomes and technical repeatability. It also examines allelic symmetry,
allele-label exchange, variance constraints, marker thinning and normal-fraction
filters. Differences between paired transition rows are reported.

PSCBS supplies an established comparator using the same arrays and additional
total-intensity information. Published sequencing segments supply a separate
assay of total copy ratio after coordinate alignment. Neither provides
error-free labels for the HMM allelic profiles. The results, including the
advantages of competing models and the instability of short segments, are
reported in the article and Supporting Information S3.
