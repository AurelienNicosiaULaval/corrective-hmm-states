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

The application uses only the animal identifier and ordered UTM coordinates.
Daily step length is the Euclidean distance between consecutive locations
within an animal. Following the public analyses of these data, locations are
treated as approximately 24 hours apart because exact observation times are
not included in the open data.
