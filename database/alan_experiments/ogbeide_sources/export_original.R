# Export a data object only; do not execute the deposited analysis code.
source_path <- "../gbif_annual/alan_experiment_sources/ogbeide2026/selected/emergence.rds"
target_path <- "../gbif_annual/alan_experiment_sources/ogbeide2026/selected/emergence_export.csv"
x <- readRDS(source_path)
stopifnot(is.data.frame(x), nrow(x) == 192L)
write.csv(x, target_path, row.names = FALSE)
