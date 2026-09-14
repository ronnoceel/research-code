# nsqip-freeflap-frailty-replication

Analysis code accompanying our replication of Jubbal et al. (*Microsurgery* 2017; PMID
28573680) on age and microsurgical free-flap outcomes, using ACS-NSQIP 2022–2024. This
repository is a companion to the paper and assumes the reader has it.

## Contents

- `analysis.R` — the analysis pipeline (cohort construction, models, sensitivity analyses).
- `nsqip_load.R` — reads the required columns across the three PUF years.

## Running

1. Obtain the ACS-NSQIP PUF 2022–2024 tab-delimited files through a participating institution.
2. Point `nsqip_load.R`'s `nsqip_path()` at your local PUF directory.
3. `Rscript analysis.R` (writes result tables to a git-ignored `outputs/`).

Requires R (≥ 4.2) with `tidyverse`, `broom`, `mice`, and `logistf`.

## Data

Case-level ACS-NSQIP data cannot be redistributed under the ACS-NSQIP Data Use Agreement and
is not included; obtain it directly from the American College of Surgeons.

## Citation

See `CITATION.cff`. Original study: Jubbal KT, Zavlin D, Suliman A. *Microsurgery.*
2017;37(8):858–864. PMID 28573680.

Released under the MIT License (repository root `LICENSE`).
