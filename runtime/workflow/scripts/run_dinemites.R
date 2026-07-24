#!/usr/bin/env Rscript
# run_dinemites.R
# ---------------------------------------------------------------------------
# Runs the DINEMITES package on prepared longitudinal input data.
# Supports three model types: simple, clustering, bayesian.
# Produces allele probabilities, molFOI, new infections, and per-subject plots.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", script_arg[[1]]), fixed = TRUE)
script_dir <- dirname(normalizePath(script_path))
source(file.path(script_dir, "lib", "analysis_common.R"))

# ── CLI arguments ──────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input", type = "character",
              help = "Path to DINEMITES input TSV (from simplseq_to_dinemites.R)"),
  make_option("--qpcr_times", type = "character", default = "",
              help = "Optional TSV of PCR-positive visits with missing genotypes"),
  make_option("--n_imputations", type = "integer", default = 10,
              help = "Complete genotype imputations for PCR-positive missing visits [default %default]"),
  make_option("--model", type = "character", default = "simple",
              help = "Model type: simple, clustering, or bayesian [default %default]"),
  make_option("--outdir", type = "character",
              help = "Output directory for DINEMITES results"),
  make_option("--n_lags", type = "integer", default = 3,
              help = "Simple model: number of recent samples to check [default %default]"),
  make_option("--t_lag", type = "character", default = "Inf",
              help = "Simple model: time window in days, or Inf for no day cutoff [default %default]"),
  make_option("--seed", type = "integer", default = 1,
              help = "Random seed for stochastic DINEMITES models [default %default]"),
  make_option("--refresh", type = "integer", default = 100,
              help = "Stan/clustering progress refresh interval [default %default]"),
  make_option("--bayesian_lag_days", type = "character", default = "30",
              help = "Comma-separated lag windows in days for Bayesian covariates [default %default]"),
  make_option("--bayesian_chains", type = "integer", default = 1,
              help = "Bayesian model Stan chains [default %default]"),
  make_option("--bayesian_parallel_chains", type = "integer", default = 1,
              help = "Bayesian model parallel Stan chains [default %default]"),
  make_option("--bayesian_iter_warmup", type = "integer", default = 500,
              help = "Bayesian model warmup iterations [default %default]"),
  make_option("--bayesian_iter_sampling", type = "integer", default = 500,
              help = "Bayesian model sampling iterations [default %default]"),
  make_option("--bayesian_adapt_delta", type = "double", default = 0.99,
              help = "Bayesian model Stan adapt_delta [default %default]"),
  make_option("--bayesian_drop_out", type = "character", default = "false",
              help = "Use the Bayesian drop-out model: true/false [default %default]"),
  make_option("--infection_general_covariates", type = "character", default = "none",
              help = "Bayesian infection-level covariates: none, auto, or comma-separated column names [default %default]"),
  make_option("--plot_width", type = "double", default = 14,
              help = "DINEMITES plot width in inches [default %default]"),
  make_option("--plot_height", type = "double", default = 0,
              help = "DINEMITES plot height in inches; 0 auto-scales [default %default]"),
  make_option("--skip_plots", action = "store_true", default = FALSE,
              help = "Skip per-subject plot generation")
)
args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args$input) || is.null(args$outdir)) {
  stop("Required arguments: --input, --outdir")
}

model_type <- tolower(args$model)
if (!model_type %in% c("simple", "clustering", "bayesian")) {
  stop("Model type must be one of: simple, clustering, bayesian. Got: ", model_type)
}

parse_t_lag <- function(value) {
  raw_value <- trimws(as.character(value))
  if (tolower(raw_value) %in% c("", "inf", "infinity", "none")) {
    return(Inf)
  }
  parsed <- suppressWarnings(as.numeric(raw_value))
  if (is.na(parsed) || parsed < 0) {
    stop("--t_lag must be a non-negative number of days or Inf. Got: ", raw_value)
  }
  parsed
}

parse_bool <- function(value) {
  raw_value <- tolower(trimws(as.character(value)))
  if (raw_value %in% c("1", "true", "yes", "y", "on")) {
    return(TRUE)
  }
  if (raw_value %in% c("0", "false", "no", "n", "off", "")) {
    return(FALSE)
  }
  stop("Boolean value must be true or false. Got: ", value)
}

parse_lag_days <- function(value) {
  raw_values <- trimws(unlist(strsplit(as.character(value), ",")))
  parsed <- suppressWarnings(as.integer(raw_values))
  if (length(parsed) == 0 || any(is.na(parsed)) || any(parsed < 1)) {
    stop("--bayesian_lag_days must contain positive whole days. Got: ", value)
  }
  sort(unique(parsed))
}

parse_covariates <- function(value, dataset) {
  raw_value <- trimws(as.character(value))
  if (tolower(raw_value) %in% c("", "none", "false", "off", "0")) {
    return(character())
  }
  if (tolower(raw_value) == "auto") {
    candidates <- intersect(
      c("covariate_season", "covariate_season_missing"),
      colnames(dataset)
    )
  } else {
    candidates <- trimws(unlist(strsplit(raw_value, ",")))
    candidates <- candidates[nchar(candidates) > 0]
  }
  if (length(candidates) == 0) {
    return(character())
  }
  missing_cols <- setdiff(candidates, colnames(dataset))
  if (length(missing_cols) > 0) {
    stop("DINEMITES covariate column(s) not found: ", paste(missing_cols, collapse = ", "))
  }
  non_numeric <- candidates[!vapply(dataset[candidates], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("DINEMITES covariate column(s) must be numeric: ", paste(non_numeric, collapse = ", "))
  }
  candidates
}

if (is.na(args$n_lags) || args$n_lags < 1) {
  stop("--n_lags must be at least 1. Got: ", args$n_lags)
}
if (is.na(args$n_imputations) || args$n_imputations < 1 || args$n_imputations > 100) {
  stop("--n_imputations must be between 1 and 100. Got: ", args$n_imputations)
}
args$t_lag <- parse_t_lag(args$t_lag)
args$bayesian_drop_out <- parse_bool(args$bayesian_drop_out)
args$bayesian_lag_days <- parse_lag_days(args$bayesian_lag_days)

if (is.na(args$seed) || args$seed < 1) {
  stop("--seed must be at least 1. Got: ", args$seed)
}
if (is.na(args$refresh) || args$refresh < 0) {
  stop("--refresh must be at least 0. Got: ", args$refresh)
}
if (is.na(args$bayesian_chains) || args$bayesian_chains < 1) {
  stop("--bayesian_chains must be at least 1. Got: ", args$bayesian_chains)
}
if (is.na(args$bayesian_parallel_chains) || args$bayesian_parallel_chains < 1) {
  stop("--bayesian_parallel_chains must be at least 1. Got: ",
       args$bayesian_parallel_chains)
}
if (is.na(args$bayesian_iter_warmup) || args$bayesian_iter_warmup < 1) {
  stop("--bayesian_iter_warmup must be at least 1. Got: ",
       args$bayesian_iter_warmup)
}
if (is.na(args$bayesian_iter_sampling) || args$bayesian_iter_sampling < 1) {
  stop("--bayesian_iter_sampling must be at least 1. Got: ",
       args$bayesian_iter_sampling)
}
if (is.na(args$bayesian_adapt_delta) ||
    args$bayesian_adapt_delta <= 0 ||
    args$bayesian_adapt_delta >= 1) {
  stop("--bayesian_adapt_delta must be greater than 0 and less than 1. Got: ",
       args$bayesian_adapt_delta)
}

cat("[DINEMITES/run] Input:", args$input, "\n")
cat("[DINEMITES/run] Model:", model_type, "\n")
cat("[DINEMITES/run] Output dir:", args$outdir, "\n")
cat("[DINEMITES/run] Simple windows: n_lags=", args$n_lags,
    ", t_lag=", ifelse(is.infinite(args$t_lag), "Inf", args$t_lag), "\n", sep = "")
cat("[DINEMITES/run] Stochastic settings: seed=", args$seed,
    ", refresh=", args$refresh, "\n", sep = "")
cat("[DINEMITES/run] Bayesian settings: lag_days=", paste(args$bayesian_lag_days, collapse = ","),
    ", chains=", args$bayesian_chains,
    ", parallel_chains=", args$bayesian_parallel_chains,
    ", warmup=", args$bayesian_iter_warmup,
    ", sampling=", args$bayesian_iter_sampling,
    ", adapt_delta=", args$bayesian_adapt_delta,
    ", drop_out=", args$bayesian_drop_out, "\n", sep = "")

build_time_axis_labels <- function(dataset) {
  time_map <- dataset %>%
    dplyr::distinct(.data$time, .keep_all = TRUE) %>%
    dplyr::arrange(.data$time)

  labels <- as.character(time_map$time)

  if ("date_full" %in% colnames(time_map)) {
    dates <- tryCatch(as.Date(time_map$date_full),
                      error = function(e) rep(as.Date(NA), nrow(time_map)))
    labels <- ifelse(!is.na(dates),
                     paste0(time_map$time, "\n", format(dates, "%d %b %Y")),
                     labels)
  } else if ("date_label" %in% colnames(time_map)) {
    date_labels <- trimws(as.character(time_map$date_label))
    parsed_dates <- suppressWarnings(as.Date(date_labels, format = "%d %b %Y"))
    labels <- ifelse(!is.na(parsed_dates),
                     paste0(time_map$time, "\n", format(parsed_dates, "%d %b %Y")),
                     ifelse(!is.na(date_labels) & nchar(date_labels) > 0,
                            paste0(time_map$time, "\n", date_labels),
                            labels))
  }

  stats::setNames(labels, as.character(time_map$time))
}

build_allele_key <- function(dataset) {
  if (!"allele" %in% colnames(dataset)) {
    return(data.frame(short_allele_id = character(),
                      locus = character(),
                      allele = character(),
                      stringsAsFactors = FALSE))
  }

  locus_values <- if ("locus" %in% colnames(dataset)) {
    dataset$locus
  } else {
    rep("Allele", nrow(dataset))
  }

  key_source <- data.frame(
    locus = as.character(locus_values),
    allele = as.character(dataset$allele),
    stringsAsFactors = FALSE
  ) %>%
    mutate(locus = ifelse(is.na(.data$locus) | trimws(.data$locus) == "",
                          "Allele", .data$locus)) %>%
    filter(!is.na(.data$allele), trimws(.data$allele) != "") %>%
    distinct(.data$locus, .data$allele) %>%
    arrange(.data$locus, .data$allele) %>%
    group_by(.data$locus) %>%
    mutate(short_allele_id = paste0(.data$locus, "-", sprintf("%02d", row_number()))) %>%
    ungroup()

  key_source %>%
    select(short_allele_id, locus, allele)
}

build_allele_axis_labeler <- function(allele_key) {
  if (is.null(allele_key) || nrow(allele_key) == 0) {
    return(function(values) as.character(values))
  }
  label_map <- stats::setNames(allele_key$short_allele_id, allele_key$allele)
  function(values) {
    full_values <- as.character(values)
    short_values <- unname(label_map[full_values])
    ifelse(is.na(short_values) | short_values == "", full_values, short_values)
  }
}

safe_subject_dir_name <- function(value) {
  safe_value <- gsub("[^A-Za-z0-9._-]+", "_", as.character(value))
  safe_value <- gsub("^_+|_+$", "", safe_value)
  ifelse(nchar(safe_value) == 0, "subject", safe_value)
}

ensure_cmdstan_path <- function() {
  current_path <- tryCatch(cmdstanr::cmdstan_path(),
                           error = function(e) "")
  if (nzchar(current_path) && dir.exists(current_path)) {
    return(current_path)
  }

  rscript_path <- Sys.which("Rscript")
  rscript_prefix <- if (nzchar(rscript_path)) {
    normalizePath(file.path(dirname(rscript_path), ".."),
                  winslash = "/", mustWork = FALSE)
  } else {
    ""
  }
  r_home_prefix <- normalizePath(file.path(R.home(), "..", ".."),
                                 winslash = "/", mustWork = FALSE)

  candidates <- unique(c(
    Sys.getenv("CMDSTAN", unset = ""),
    file.path(Sys.getenv("CONDA_PREFIX", unset = ""), "bin", "cmdstan"),
    file.path(r_home_prefix, "bin", "cmdstan"),
    file.path(rscript_prefix, "bin", "cmdstan")
  ))
  candidates <- candidates[nzchar(candidates) & dir.exists(candidates)]
  if (length(candidates) == 0) {
    stop("[DINEMITES/run] ERROR: CmdStanR is installed, but no CmdStan directory was found.")
  }

  cmdstanr::set_cmdstan_path(candidates[1])
  cat("[DINEMITES/run] CmdStan path:", candidates[1], "\n")
  invisible(candidates[1])
}

compile_dinemites_stan_model <- function(model_type, bayesian_drop_out = FALSE) {
  package_path <- system.file(package = "dinemites")
  model_files <- instantiate::stan_package_model_files(package_path)
  model_pattern <- if (model_type == "clustering") {
    "model_infection_probabilities_clusters\\.stan$"
  } else if (bayesian_drop_out) {
    "model_infection_probabilities_bayesian_drop_out\\.stan$"
  } else {
    "model_infection_probabilities_bayesian\\.stan$"
  }
  selected_models <- model_files[grepl(model_pattern, model_files)]
  if (length(selected_models) != 1) {
    stop("[DINEMITES/run] ERROR: Could not find packaged Stan model for ",
         model_type, ".")
  }
  cat("[DINEMITES/run] Preparing Stan model:", basename(selected_models), "\n")
  instantiate::stan_package_compile(models = selected_models, quiet = TRUE)
  invisible(selected_models)
}

PREVALENCE_STRIP_LABEL <- "Allele prevalence (% visits)"
READABLE_PLOT_THEME <- theme(
  text = element_text(size = 12, color = "grey10"),
  plot.title = element_text(size = 18, color = "grey10"),
  axis.title = element_text(size = 13, color = "grey10"),
  axis.text = element_text(size = 11, color = "grey10"),
  strip.text = element_text(size = 12, color = "grey10"),
  legend.text = element_text(size = 11, color = "grey10"),
  legend.title = element_text(size = 12, color = "grey10")
)

prevalence_strip <- function() {
  display_label <- sub(" \\(", "\n(", PREVALENCE_STRIP_LABEL)
  ggplot(data.frame(l = display_label, x = 1, y = 1)) +
    geom_text(aes(.data$x, .data$y, label = .data$l), angle = 270, size = 5,
              lineheight = 0.95, color = "grey10") +
    theme_void() +
    coord_cartesian(clip = "off")
}

apply_time_axis_to_plot <- function(plot_out, x_breaks, x_labels, x_limit_max,
                                    allele_labeler = function(values) as.character(values)) {
  patchwork_plots <- list()
  if (inherits(plot_out, "patchwork")) {
    patchwork_plots <- plot_out$plots
    if (length(patchwork_plots) == 0 && !is.null(plot_out$patches$plots)) {
      patchwork_plots <- plot_out$patches$plots
    }
  }

  x_scale <- scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    limits = c(min(x_breaks), x_limit_max),
    guide = guide_axis(n.dodge = 2)
  )
  x_label <- labs(x = "Study day\nCollection date")
  x_theme <- theme(axis.text.x = element_text(size = 11, color = "grey10", lineheight = 0.95))
  allele_axis_theme <- theme(
    axis.text.y = element_text(size = 10, color = "grey10"),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 12, color = "grey10")
  )
  top_axis_theme <- theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )

  if (inherits(plot_out, "patchwork") && length(patchwork_plots) >= 3) {
    top_plot <- suppressMessages(patchwork_plots[[1]] + x_scale + labs(x = NULL) +
                                   READABLE_PLOT_THEME + top_axis_theme)
    allele_plot <- suppressMessages(
      patchwork_plots[[2]] + x_scale + x_label + labs(y = "Allele ID") +
        scale_y_discrete(labels = allele_labeler) + READABLE_PLOT_THEME +
        x_theme + allele_axis_theme
    )
    prevalence_label <- prevalence_strip()

    return(patchwork::wrap_plots(
      list(top_plot, allele_plot, prevalence_label),
      design = c(patchwork::area(1, 1),
                 patchwork::area(2, 1),
                 patchwork::area(2, 2, 2, 2)),
      heights = c(1, 5),
      widths = c(18, 1.5)
    ))
  }

  if (inherits(plot_out, "patchwork") && length(patchwork_plots) == 2) {
    top_plot <- suppressMessages(patchwork_plots[[1]] + x_scale + labs(x = NULL) +
                                   READABLE_PLOT_THEME + top_axis_theme)
    allele_plot <- suppressMessages(
      patchwork_plots[[2]] + x_scale + x_label + labs(y = "Allele ID") +
        scale_y_discrete(labels = allele_labeler) + READABLE_PLOT_THEME +
        x_theme + allele_axis_theme
    )
    return(patchwork::wrap_plots(
      list(top_plot, allele_plot, prevalence_strip()),
      design = c(patchwork::area(1, 1),
                 patchwork::area(2, 1),
                 patchwork::area(2, 2, 2, 2)),
      heights = c(1, 5),
      widths = c(18, 1.5)
    ))
  }

  suppressMessages(
    plot_out + x_scale + x_label + labs(y = "Allele ID") +
      scale_y_discrete(labels = allele_labeler) + READABLE_PLOT_THEME +
      x_theme + allele_axis_theme
  )
}

# ── Check DINEMITES package availability ───────────────────────────────────
if (!requireNamespace("dinemites", quietly = TRUE)) {
  stop("[DINEMITES/run] ERROR: the managed DINEMITES package is missing. ",
       "Re-run malaria-amplicon-nf runtime setup.")
}

library(dinemites)

# Check Stan package availability for Bayesian/clustering models.
if (model_type %in% c("bayesian", "clustering")) {
  stan_package_names <- c("instantiate", "cmdstanr", "rstan", "posterior")
  if (model_type == "clustering") {
    stan_package_names <- c(stan_package_names, "linkcomm")
  }
  missing_stan_packages <- stan_package_names[
    !vapply(stan_package_names, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_stan_packages) > 0) {
    stop("[DINEMITES/run] ERROR: ", model_type,
         " model requires missing R package(s): ",
         paste(missing_stan_packages, collapse = ", "))
  }
  invisible(ensure_cmdstan_path())
  invisible(compile_dinemites_stan_model(model_type, args$bayesian_drop_out))
}

# ── 1. Read input data ────────────────────────────────────────────────────
dataset <- read.csv(
  args$input,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
required_columns <- c("allele", "time", "subject")
missing_columns <- setdiff(required_columns, colnames(dataset))
if (length(missing_columns) > 0) {
  stop("[DINEMITES/run] Missing required input columns: ",
       paste(missing_columns, collapse = ", "))
}

# DINEMITES accepts numeric allele/locus identifiers, but our output joins and
# filenames need stable character keys across every model and input source.
dataset <- dataset %>%
  mutate(
    allele = na_if(trimws(as.character(.data$allele)), ""),
    subject = as.character(.data$subject),
    time = suppressWarnings(as.numeric(.data$time))
  )
if ("locus" %in% colnames(dataset)) {
  dataset$locus <- as.character(dataset$locus)
}
if (any(is.na(dataset$time))) {
  stop("[DINEMITES/run] time must contain numeric day values.")
}
infection_general_covariates <- parse_covariates(args$infection_general_covariates, dataset)
if (length(infection_general_covariates) == 0) {
  # DINEMITES distinguishes NULL (no covariates) from character(0), which its
  # numeric validation interprets as a non-numeric empty data frame.
  infection_general_covariates <- NULL
}

cat("[DINEMITES/run] Loaded", nrow(dataset), "rows (",
    length(unique(dataset$subject)), "subjects,",
    length(unique(stats::na.omit(dataset$allele))), "alleles )\n")
cat("[DINEMITES/run] Visit placeholders:", sum(is.na(dataset$allele)),
    "(preserved as visit metadata; never treated as alleles)\n")
cat("[DINEMITES/run] Infection general covariates:",
    ifelse(length(infection_general_covariates) > 0,
           paste(infection_general_covariates, collapse = ", "),
           "none"),
    "\n")

# ── 2. Fill in dataset (complete allele × subject × time grid) ─────────────
cat("[DINEMITES/run] Filling in dataset (complete grid)...\n")
dataset_filled <- fill_in_dataset(dataset)
if (any(is.na(dataset_filled$allele) |
        trimws(as.character(dataset_filled$allele)) == "")) {
  stop("[DINEMITES/run] Internal error: the filled dataset contains a blank allele.")
}
cat("[DINEMITES/run] Filled dataset:", nrow(dataset_filled), "rows\n")

qpcr_times <- data.frame(subject = character(), time = numeric())
has_qpcr_imputation <- FALSE
imputed_datasets <- NULL
probability_matrix <- NULL
if (nzchar(trimws(args$qpcr_times)) && file.exists(args$qpcr_times)) {
  qpcr_times <- read.delim(
    args$qpcr_times,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
  missing_qpcr_columns <- setdiff(c("subject", "time"), colnames(qpcr_times))
  if (length(missing_qpcr_columns) > 0) {
    stop("[DINEMITES/run] qPCR-only input is missing column(s): ",
         paste(missing_qpcr_columns, collapse = ", "))
  }
  qpcr_times <- qpcr_times %>%
    transmute(
      subject = as.character(.data$subject),
      time = suppressWarnings(as.numeric(.data$time))
    ) %>%
    filter(!is.na(.data$time), nzchar(trimws(.data$subject))) %>%
    distinct()
}

if (nrow(qpcr_times) > 0) {
  cat("[DINEMITES/run] Marking", nrow(qpcr_times),
      "PCR-positive visits with missing genotypes...\n")
  dataset_filled <- add_qpcr_times(dataset_filled, qpcr_times = qpcr_times)
  detected_cores <- suppressWarnings(parallel::detectCores())
  if (is.na(detected_cores) || detected_cores < 2) {
    detected_cores <- 2L
  }
  n_cores <- max(1L, min(4L, detected_cores - 1L))
  set.seed(args$seed)
  cat("[DINEMITES/run] Creating", args$n_imputations,
      "complete genotype imputations on", n_cores, "cores...\n")
  imputed_datasets <- impute_dataset(
    dataset_filled,
    n_imputations = args$n_imputations,
    n_cores = n_cores,
    verbose = TRUE
  )
  dataset_filled <- add_probability_present(dataset_filled, imputed_datasets)
  has_qpcr_imputation <- TRUE
}

prepare_bayesian_dataset <- function(model_dataset) {
  lag_columns <- paste0("lag_", args$bayesian_lag_days)
  lag_infection_columns <- paste0("lag_infection_", args$bayesian_lag_days)
  prepared <- model_dataset %>%
    add_present_infection() %>%
    add_persistent_column() %>%
    add_persistent_infection()
  for (lag_days in args$bayesian_lag_days) {
    prepared <- prepared %>%
      add_lag_column(lag_time = lag_days) %>%
      add_lag_infection(lag_time = lag_days)
  }
  list(
    dataset = prepared,
    lag_columns = lag_columns,
    lag_infection_columns = lag_infection_columns
  )
}

run_selected_model <- function(model_dataset, run_seed) {
  if (model_type == "simple") {
    return(list(
      dataset = model_dataset,
      result = determine_probabilities_simple(
        model_dataset,
        n_lags = args$n_lags,
        t_lag = args$t_lag
      )
    ))
  }
  if (model_type == "clustering") {
    return(list(
      dataset = model_dataset,
      result = determine_probabilities_clustering(
        model_dataset,
        refresh = args$refresh,
        seed = run_seed
      )
    ))
  }

  bayesian <- prepare_bayesian_dataset(model_dataset)
  list(
    dataset = bayesian$dataset,
    result = determine_probabilities_bayesian(
      bayesian$dataset,
      infection_persistence_covariates = c(
        "persistent_infection", bayesian$lag_infection_columns
      ),
      infection_general_covariates = infection_general_covariates,
      alleles_persistence_covariates = c("persistent", bayesian$lag_columns),
      chains = args$bayesian_chains,
      parallel_chains = args$bayesian_parallel_chains,
      iter_warmup = args$bayesian_iter_warmup,
      iter_sampling = args$bayesian_iter_sampling,
      refresh = args$refresh,
      adapt_delta = args$bayesian_adapt_delta,
      seed = run_seed,
      drop_out = args$bayesian_drop_out
    )
  )
}

# ── 3. Run selected model ─────────────────────────────────────────────────
model_runs <- if (has_qpcr_imputation) args$n_imputations else 1L
cat("[DINEMITES/run] Running", model_type, "model across", model_runs,
    ifelse(model_runs == 1, "dataset...\n", "imputed datasets...\n"))
t_start <- Sys.time()
analysis_dataset <- dataset_filled
model_outputs <- vector("list", model_runs)
probability_columns <- vector("list", model_runs)
for (index in seq_len(model_runs)) {
  model_dataset <- dataset_filled
  if (has_qpcr_imputation) {
    model_dataset$present <- imputed_datasets[, index]
  }
  cat("[DINEMITES/run] Model dataset", index, "of", model_runs, "\n")
  model_outputs[[index]] <- run_selected_model(model_dataset, args$seed + index - 1L)
  probability_columns[[index]] <- model_outputs[[index]]$result$probability_new
  if (length(probability_columns[[index]]) != nrow(model_outputs[[index]]$dataset)) {
    stop("[DINEMITES/run] Model returned ", length(probability_columns[[index]]),
         " probabilities for ", nrow(model_outputs[[index]]$dataset), " dataset rows.")
  }
}

dataset_filled <- if (model_type == "bayesian") {
  model_outputs[[1]]$dataset
} else {
  analysis_dataset
}
dataset_filled$present <- analysis_dataset$present
if ("probability_present" %in% colnames(analysis_dataset)) {
  dataset_filled$probability_present <- analysis_dataset$probability_present
}
results <- model_outputs[[1]]$result
if (has_qpcr_imputation) {
  probability_matrix <- do.call(cbind, probability_columns)
  colnames(probability_matrix) <- paste0("imputation_", seq_len(model_runs))
  dataset_filled <- add_probability_new(
    dataset_filled,
    probability_matrix,
    "probability_new"
  )
  results$probability_new <- dataset_filled$probability_new
} else {
  dataset_filled$probability_new <- results$probability_new
}

t_elapsed <- difftime(Sys.time(), t_start, units = "secs")
cat("[DINEMITES/run] Model completed in", round(as.numeric(t_elapsed), 1), "seconds\n")

# ── 4. Add probability column to dataset ───────────────────────────────────
# determine_probabilities_* returns a list with $probability_new (vector)
# and $fit (model object or NULL). We attach the probabilities to the dataset.
cat("[DINEMITES/run] Attaching probabilities to dataset...\n")
if (length(results$probability_new) != nrow(dataset_filled)) {
  stop("[DINEMITES/run] Model returned ", length(results$probability_new),
       " probabilities for ", nrow(dataset_filled), " dataset rows.")
}
dataset_filled$probability_new <- results$probability_new
allele_key <- build_allele_key(dataset_filled)
if (nrow(allele_key) > 0 && "locus" %in% colnames(dataset_filled)) {
  dataset_filled <- dataset_filled %>%
    left_join(allele_key, by = c("locus", "allele"))
} else {
  dataset_filled$short_allele_id <- NA_character_
}

# ── 5. Calculate new infections and molFOI ─────────────────────────────────
cat("[DINEMITES/run] Estimating new infections...\n")
if (has_qpcr_imputation) {
  estimated_new_infections_for_plot <- estimate_new_infections(
    dataset_filled,
    imputation_mat = imputed_datasets,
    probability_mat = probability_matrix
  )
  new_infection_matrix <- as.matrix(estimated_new_infections_for_plot)
  new_infection_subjects <- rownames(new_infection_matrix)
  if (is.null(new_infection_subjects)) {
    new_infection_subjects <- sort(unique(dataset_filled$subject))
  }
  new_infections <- data.frame(
    subject = new_infection_subjects,
    new_infections_mean = rowMeans(new_infection_matrix),
    new_infections_sd = if (ncol(new_infection_matrix) > 1) {
      apply(new_infection_matrix, 1, stats::sd)
    } else {
      rep(0, nrow(new_infection_matrix))
    },
    new_infections_lower = apply(
      new_infection_matrix, 1, stats::quantile, probs = 0.025
    ),
    new_infections_upper = apply(
      new_infection_matrix, 1, stats::quantile, probs = 0.975
    ),
    stringsAsFactors = FALSE
  )
} else {
  estimated_new_infections_for_plot <- estimate_new_infections(dataset_filled)
  new_infections <- estimated_new_infections_for_plot
  if (!"subject" %in% colnames(new_infections)) {
    new_infections$subject <- rownames(new_infections)
    new_infections <- new_infections %>%
      select(subject, everything())
  }
  rownames(new_infections) <- NULL
}

cat("[DINEMITES/run] Calculating molecular FOI...\n")
molfoi <- compute_molFOI(dataset_filled, method = "sum_then_max")
if (!"subject" %in% colnames(molfoi)) {
  molfoi$subject <- rownames(molfoi)
  molfoi <- molfoi %>%
    select(subject, everything())
}
rownames(molfoi) <- NULL
subjects <- unique(dataset_filled$subject)

# ── 6. Write outputs ──────────────────────────────────────────────────────
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
plots_dir <- file.path(args$outdir, "dinemites_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
subjects_dir <- file.path(args$outdir, "dinemites_subjects")
dir.create(subjects_dir, showWarnings = FALSE, recursive = TRUE)

if (has_qpcr_imputation) {
  imputation_key <- dataset_filled %>%
    transmute(
      row_index = seq_len(nrow(dataset_filled)),
      subject = as.character(.data$subject),
      time = as.numeric(.data$time),
      locus = as.character(.data$locus),
      allele = as.character(.data$allele)
    )
  write.table(
    imputation_key,
    file = file.path(args$outdir, "dinemites_imputation_row_key.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  write.table(
    imputed_datasets,
    file = file.path(args$outdir, "dinemites_imputed_presence_matrix.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  write.table(
    probability_matrix,
    file = file.path(args$outdir, "dinemites_probability_matrix.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  saveRDS(
    imputed_datasets,
    file = file.path(args$outdir, "dinemites_imputed_presence_matrix.rds")
  )
  cat("[DINEMITES/run] Wrote reproducible imputation matrices.\n")
}

# Allele probabilities
allele_probs_path <- file.path(args$outdir, "dinemites_allele_probabilities.tsv")
write.table(dataset_filled, file = allele_probs_path, sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("[DINEMITES/run] Wrote allele probabilities:", allele_probs_path, "\n")

# Allele key for static plot row IDs
allele_key_path <- file.path(args$outdir, "dinemites_allele_key.tsv")
write.table(allele_key, file = allele_key_path, sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("[DINEMITES/run] Wrote allele key:", allele_key_path, "\n")

# New infections per subject
new_inf_path <- file.path(args$outdir, "dinemites_new_infections.tsv")
write.table(new_infections, file = new_inf_path, sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("[DINEMITES/run] Wrote new infections:", new_inf_path, "\n")

# molFOI per subject
molfoi_path <- file.path(args$outdir, "dinemites_molfoi.tsv")
write.table(molfoi, file = molfoi_path, sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("[DINEMITES/run] Wrote molFOI:", molfoi_path, "\n")

# Per-subject outputs for patient-facing review.
for (subj in subjects) {
  subject_dir <- file.path(subjects_dir, safe_subject_dir_name(subj))
  dir.create(subject_dir, showWarnings = FALSE, recursive = TRUE)

  subject_dataset_all <- dataset_filled %>%
    filter(.data$subject == subj)
  subject_observed_alleles <- subject_dataset_all %>%
    filter(!is.na(.data$allele), trimws(as.character(.data$allele)) != "",
           .data$present > 0) %>%
    distinct(.data$locus, .data$allele)

  subject_dataset <- subject_dataset_all[0, ]
  if (nrow(subject_observed_alleles) > 0 && "locus" %in% colnames(subject_dataset_all)) {
    subject_dataset <- subject_dataset_all %>%
      semi_join(subject_observed_alleles, by = c("locus", "allele"))
  }

  subject_key <- allele_key
  if (nrow(subject_observed_alleles) > 0 && nrow(allele_key) > 0) {
    subject_key <- allele_key %>%
      semi_join(subject_observed_alleles, by = c("locus", "allele"))
  } else if (nrow(allele_key) > 0) {
    subject_key <- allele_key[0, ]
  }

  write.table(subject_dataset,
              file = file.path(subject_dir, "dinemites_allele_probabilities.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(subject_key,
              file = file.path(subject_dir, "dinemites_allele_key.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(new_infections %>% filter(.data$subject == subj),
              file = file.path(subject_dir, "dinemites_new_infections.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(molfoi %>% filter(.data$subject == subj),
              file = file.path(subject_dir, "dinemites_molfoi.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  subject_summary_lines <- c(
    "{",
    paste0('  "model_type": "', model_type, '",'),
    paste0('  "subject": "', subj, '",'),
    paste0('  "n_alleles": ', nrow(subject_observed_alleles), ','),
    paste0('  "n_timepoints": ', length(unique(subject_dataset_all$time)), ','),
    paste0('  "plot": "subject_', subj, '.png"'),
    "}"
  )
  writeLines(subject_summary_lines, file.path(subject_dir, "dinemites_summary.json"))
  cat("[DINEMITES/run] Wrote subject outputs:", subject_dir, "\n")
}

# ── 7. Generate per-subject plots ──────────────────────────────────────────
cat("[DINEMITES/run] Generating per-subject plots...\n")
empty_treatments <- data.frame(subject = character(), time = numeric())
allele_labeler <- build_allele_axis_labeler(allele_key)
plot_width <- max(args$plot_width, 11)

if (!args$skip_plots) tryCatch({
  plot_list <- plot_dataset(dataset_filled,
                            treatments = empty_treatments,
                            estimated_new_infections = estimated_new_infections_for_plot,
                            output = NULL,
                            height = 10,
                            width = plot_width)
  for (subj in subjects) {
    plot_path <- file.path(plots_dir, paste0("subject_", subj, ".png"))
    plot_out <- plot_list[[as.character(subj)]]
    if (!is.null(plot_out)) {
      subject_data <- dataset_filled %>% filter(.data$subject == subj)
      subject_breaks <- sort(unique(subject_data$time))
      subject_axis_labels <- build_time_axis_labels(subject_data)
      subject_labels <- unname(subject_axis_labels[as.character(subject_breaks)])
      subject_limit_max <- ifelse(
        max(subject_breaks) == min(subject_breaks),
        max(subject_breaks) + 1,
        max(subject_breaks) + max(7, 0.08 * diff(range(subject_breaks)))
      )
      subject_present_alleles <- length(unique(
        subject_data$allele[subject_data$present > 0 & !is.na(subject_data$allele)]
      ))
      subject_loci <- if ("locus" %in% colnames(subject_data)) {
        max(1, length(unique(na.omit(subject_data$locus))))
      } else {
        1
      }
      subject_plot_height <- if (args$plot_height > 0) {
        args$plot_height
      } else {
        min(18, max(5.5, 4.2 + (0.38 * subject_present_alleles) + (0.35 * subject_loci)))
      }
      plot_out <- apply_time_axis_to_plot(
        plot_out, subject_breaks, subject_labels, subject_limit_max, allele_labeler
      )
      ggsave(plot_path, plot = plot_out, height = subject_plot_height,
             width = plot_width, dpi = 300, limitsize = FALSE)
      cat("[DINEMITES/run] Plot:", plot_path, "\n")
      subject_plot_path <- file.path(subjects_dir, safe_subject_dir_name(subj),
                                     paste0("subject_", subj, ".png"))
      invisible(file.copy(plot_path, subject_plot_path, overwrite = TRUE))
    }
  }
}, error = function(e) {
  cat("[DINEMITES/run] WARNING: Could not generate DINEMITES plots:",
      conditionMessage(e), "\n")
})

# ── 8. Summary report ─────────────────────────────────────────────────────
summary_path <- file.path(args$outdir, "dinemites_summary.json")
model_diagnostics <- list(fits_evaluated = model_runs)
if (model_type == "bayesian") {
  diagnostic_rows <- lapply(seq_along(model_outputs), function(index) {
    fit <- model_outputs[[index]]$result$fit
    fit_summary <- tryCatch(fit$summary(), error = function(e) NULL)
    sampler_summary <- tryCatch(fit$diagnostic_summary(), error = function(e) NULL)
    finite_rhat <- if (!is.null(fit_summary)) {
      fit_summary$rhat[is.finite(fit_summary$rhat)]
    } else {
      numeric()
    }
    finite_ess <- if (!is.null(fit_summary)) {
      fit_summary$ess_bulk[is.finite(fit_summary$ess_bulk)]
    } else {
      numeric()
    }
    list(
      imputation = index,
      max_rhat = if (length(finite_rhat)) max(finite_rhat) else NA_real_,
      min_ess_bulk = if (length(finite_ess)) min(finite_ess) else NA_real_,
      divergent_transitions = if (!is.null(sampler_summary)) {
        sum(sampler_summary$num_divergent, na.rm = TRUE)
      } else {
        NA_real_
      },
      max_treedepth_transitions = if (!is.null(sampler_summary)) {
        sum(sampler_summary$num_max_treedepth, na.rm = TRUE)
      } else {
        NA_real_
      }
    )
  })
  rhat_values <- vapply(diagnostic_rows, function(row) row$max_rhat, numeric(1))
  ess_values <- vapply(diagnostic_rows, function(row) row$min_ess_bulk, numeric(1))
  divergence_values <- vapply(
    diagnostic_rows, function(row) row$divergent_transitions, numeric(1)
  )
  treedepth_values <- vapply(
    diagnostic_rows, function(row) row$max_treedepth_transitions, numeric(1)
  )
  model_diagnostics$max_rhat <- if (any(is.finite(rhat_values))) {
    max(rhat_values, na.rm = TRUE)
  } else {
    NA_real_
  }
  model_diagnostics$min_ess_bulk <- if (any(is.finite(ess_values))) {
    min(ess_values, na.rm = TRUE)
  } else {
    NA_real_
  }
  model_diagnostics$divergent_transitions <- sum(divergence_values, na.rm = TRUE)
  model_diagnostics$max_treedepth_transitions <- sum(treedepth_values, na.rm = TRUE)
  model_diagnostics$per_imputation <- diagnostic_rows
}
summary_data <- list(
  model_type     = model_type,
  seed           = args$seed,
  full_visit_calendar = TRUE,
  qpcr_positive_missing_genotypes = nrow(qpcr_times),
  used_genotype_imputation = has_qpcr_imputation,
  n_imputations = if (has_qpcr_imputation) args$n_imputations else 0,
  refresh        = args$refresh,
  bayesian_lag_days = paste(args$bayesian_lag_days, collapse = ","),
  bayesian_chains = args$bayesian_chains,
  bayesian_parallel_chains = args$bayesian_parallel_chains,
  bayesian_iter_warmup = args$bayesian_iter_warmup,
  bayesian_iter_sampling = args$bayesian_iter_sampling,
  bayesian_adapt_delta = args$bayesian_adapt_delta,
  bayesian_drop_out = args$bayesian_drop_out,
  infection_general_covariates = paste(infection_general_covariates, collapse = ","),
  n_subjects     = length(subjects),
  n_alleles      = length(unique(stats::na.omit(dataset_filled$allele))),
  n_timepoints   = length(unique(dataset_filled$time)),
  runtime_secs   = round(as.numeric(t_elapsed), 1),
  subjects       = subjects,
  plots          = paste0("subject_", subjects, ".png"),
  diagnostics    = model_diagnostics
)

write_json_file(summary_data, summary_path)

cat("[DINEMITES/run] Wrote summary:", summary_path, "\n")
cat("[DINEMITES/run] Done.\n")
