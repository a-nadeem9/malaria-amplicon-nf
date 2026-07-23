# Shared runtime helpers for the core SIMPLseq R scripts.

workflow_threads <- function(default = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv("SIMPLSEQ_THREADS", unset = as.character(default))))
  if (is.na(value) || value < 1L) default else value
}

register_bounded_backend <- function() {
  detected <- suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)))
  if (is.na(detected) || detected < 1L) detected <- workflow_threads()
  workers <- max(1L, min(workflow_threads(), detected))
  doMC::registerDoMC(cores = workers)
  workers
}

require_input_file <- function(path, label) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop(label, " not found: ", ifelse(is.null(path), "<not provided>", path))
  }
  invisible(path)
}
