# Shared helpers for optional downstream analyses.

analysis_palette <- list(
  teal = "#20988c",
  teal_dark = "#005c68",
  teal_pale = "#eef8f7",
  pink = "#ec1850",
  pink_pale = "#fff1f5",
  neutral = "#e8eeee"
)

require_columns <- function(dataset, required, context) {
  missing <- setdiff(required, colnames(dataset))
  if (length(missing) > 0) {
    stop("[", context, "] Missing required columns: ", paste(missing, collapse = ", "))
  }
}

write_json_file <- function(value, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to write analysis metadata.")
  }
  jsonlite::write_json(value, path, pretty = TRUE, auto_unbox = TRUE, na = "null")
}

analysis_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#edf2f2", linewidth = 0.3),
      axis.text = ggplot2::element_text(color = analysis_palette$teal_dark),
      axis.title = ggplot2::element_text(color = analysis_palette$teal_dark),
      plot.title = ggplot2::element_text(color = analysis_palette$teal_dark, face = "bold")
    )
}
