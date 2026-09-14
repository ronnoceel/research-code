# nsqip_load.R -- pooled reader for the local ACS-NSQIP Participant Use Files.
# ---------------------------------------------------------------------------
# The PUFs live git-ignored at studies/_master/data/nsqip/puf{22,23,24}/ as
# tab-delimited, latin1-encoded flat files (~1M cases/yr). This is the NSQIP
# analogue of common/seer_extract.R: read ONLY the columns a study needs across
# the requested years and row-bind them. Study-specific recodes live in each
# analysis.R (per CLAUDE.md), not here.
#
# Gotchas baked in from the data:
#   * encoding is latin1 (occurrence labels carry non-breaking spaces);
#   * numeric labs/times use -99 as a MISSING sentinel (handle per-variable);
#   * occurrence variables use "No Complication" as their null level;
#   * column sets drift by year -- every column you request must exist in every
#     requested year or the read errors (that is deliberate: catch drift early).
#
# Usage (from inside a study folder):
#   source(file.path("..","..","common","nsqip_load.R"))
#   raw <- load_nsqip(cols = c("Age","SEX","PRHCT", ...))          # all 3 yrs
#   raw <- load_nsqip(cols = c(...), years = c(2023, 2024))        # subset
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

nsqip_path <- function(year) {
  yy <- substr(as.character(year), 3, 4)
  file.path("..", "_master", "data", "nsqip", paste0("puf", yy),
            paste0("acs_nsqip_puf", yy, ".txt"))
}

# Read `cols` (PUFYEAR is always added) from each year and row-bind. Columns are
# read as-is; coerce/recode downstream. Forces character read then re-parses to
# sidestep cross-year type-guess mismatches on row-bind.
load_nsqip <- function(cols, years = c(2022, 2023, 2024)) {
  cols <- unique(c("PUFYEAR", cols))
  parts <- lapply(years, function(y) {
    p <- nsqip_path(y)
    if (!file.exists(p)) stop("missing PUF for ", y, ": ", p)
    message("reading ", p)
    read_tsv(
      p,
      col_select = all_of(cols),
      col_types  = cols(.default = col_character()),
      locale     = locale(encoding = "latin1"),
      na         = c("", "NA"),
      progress   = FALSE
    )
  })
  out <- bind_rows(parts)
  # re-parse numeric-looking columns (readr's type_convert, but silent)
  suppressMessages(readr::type_convert(out, guess_integer = TRUE))
}
