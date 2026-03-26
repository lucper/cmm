library(tidyverse)

# ---------- helpers ----------

# Read "A B" edge list (whitespace separated) into a tibble(from, to)
read_edges <- function(path) {
  readr::read_table2(
    file = path,
    col_names = c("from", "to"),
    col_types = cols(.default = col_character())
  )
}

# Read FASTA where header is >ID and sequence may span multiple lines.
# Returns tibble(id, label) where label is concatenated sequence (no spaces).
read_fasta_labels <- function(path) {
  lines <- readr::read_lines(path)

  tibble(line = lines) %>%
    mutate(is_header = str_starts(line, ">")) %>%
    mutate(id = if_else(is_header, str_remove(line, "^>"), NA_character_)) %>%
    tidyr::fill(id) %>%
    filter(!is_header) %>%
    mutate(line = str_replace_all(line, "\\s+", "")) %>%
    group_by(id) %>%
    summarise(label = paste0(line, collapse = ""), .groups = "drop")
}

# Build mapping id -> integer code (1..n) using order of first appearance.
# If you prefer alphabetical: replace distinct() with distinct(id) %>% arrange(id)
build_id_map <- function(edges, labels) {
  ids_in_order <- c(edges$from, edges$to, labels$id)

  tibble(id = ids_in_order) %>%
    filter(!is.na(id), id != "") %>%
    distinct(id) %>%                      # keeps first occurrence
    mutate(code = row_number())
}

# ---------- main conversion ----------

convert_graph_formats <- function(edge_path,
                                  label_path,
                                  out_edge_path,
                                  out_label_path) {
  edges  <- read_edges(edge_path)
  labels <- read_fasta_labels(label_path)

  id_map <- build_id_map(edges, labels)

  # Re-encode edges
  edges_out <- edges %>%
    left_join(id_map, by = c("from" = "id")) %>%
    rename(left = code) %>%
    left_join(id_map, by = c("to" = "id")) %>%
    rename(right = code) %>%
    select(left, right)

  # Re-encode labels
  labels_out <- labels %>%
    left_join(id_map, by = c("id" = "id")) %>%
    transmute(left = code, right = label) %>%
    arrange(left)

  # Optional: sanity checks
  if (anyNA(edges_out$left) || anyNA(edges_out$right)) {
    stop("Some edge endpoints were not mapped (unexpected).")
  }
  if (anyNA(labels_out$left)) {
    stop("Some labels were not mapped (unexpected).")
  }

  # Write with semicolon separator and header "left;right"
  readr::write_delim(edges_out,  out_edge_path,  delim = ";")
  readr::write_delim(labels_out, out_label_path, delim = ";")

  invisible(list(
    id_map = id_map,
    edges_out = edges_out,
    labels_out = labels_out
  ))
}

res <- convert_graph_formats(
  edge_path     = "./data/aj-yeast/PEabovecutoff.lis",
  label_path    = "./data/aj-yeast/PEabovecutoff_filt2.f",
  out_edge_path = "./data/aj-yeast/edge_list.csv",
  out_label_path= "./data/aj-yeast/node_labels.csv"
)

head(res$edges_out)
head(res$labels_out)