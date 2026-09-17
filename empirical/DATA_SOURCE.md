# Elk GPS data used in the empirical application

The analysis reads the `elk_data` object distributed with R package `moveHMM`.
The original coordinate table is not redistributed in this repository. At run
time, `numerics/scripts/run_elk_application.R` loads the object directly into
memory and verifies its structure and SHA-256 digest before fitting any model.
No CSV containing the original coordinates or derived coordinate endpoints is
created.
The archived results were generated from `moveHMM` version 1.10. The object
contains 735 ordered GPS locations for four elk (`elk-115`, `elk-163`,
`elk-287`, and `elk-363`) released in east-central Ontario, Canada.

The data originate from:

Morales, J. M., Haydon, D. T., Frair, J., Holsinger, K. E., and Fryxell,
J. M. (2004). Extracting more out of relocation data: building movement
models as mixtures of random walks. *Ecology*, 85, 2436-2445.
https://doi.org/10.1890/03-0269

They are distributed with:

Michelot, T., Langrock, R., and Patterson, T. A. (2016). moveHMM: an R package
for the statistical modelling of animal movement data using hidden Markov
models. *Methods in Ecology and Evolution*, 7, 1308-1315.
https://doi.org/10.1111/2041-210X.12578

The `moveHMM` package is distributed under GPL-3. The analysis input is loaded
with:

```r
library(moveHMM)
data(elk_data)
```

Expected SHA-256 of the serialized R object:

```text
b18bb2b3a6f2b71bd7e95efcb59d9f5c907bd88a7664989e17ecea5ccbb1ea49
```

The primary fitted response uses the animal identifier and ordered UTM coordinates.
Daily step length is the Euclidean distance between consecutive locations
within an animal. Following the public analyses of these data, locations are
treated as approximately 24 hours apart because exact observation times are
not included in the open data.


## Matched supplementary measurements for the second revision

Morales et al. (2004), Supplement 1, archived with the original publication:
https://doi.org/10.6084/m9.figshare.3523667.v1
File: https://ndownloader.figshare.com/files/5594180
SHA-256: `36d423b2832aecd929f7d539da8779965b93a48fa6aa1321fdb7b194e68db2a5`.
The archive metadata identify the supplementary deposit as CC0.

After removing 370 blank rows, all 735 IDs and coordinate pairs match the
package data exactly. The second revision uses published turning angles and
distance to open forest as auxiliary descriptive measurements. These are not
independently observed behavioural labels. A targeted check for elk 163 uses
the published interval-adjusted movement rates. No raw coordinate file is
redistributed in the reproducibility package. The source file is downloaded
temporarily at run time and verified before use.

Parton, Alison, and Blackwell, Paul G. (2017). Bayesian inference for multistate
step and turn animal movement in continuous time. JABES, 22, 373-392.
https://doi.org/10.1007/s13253-017-0286-5
Their treatment of these data explains the approximately 24-hour assumption
and the absence of exact observation times.
