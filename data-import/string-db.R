if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("STRINGdb", quietly = TRUE)) { BiocManager::install("STRINGdb") }
if (!requireNamespace("R.untils", quietly = TRUE)) { BiocManager::install("R.utils") }

library(data.table)
library(STRINGdb)
library(R.utils)

# --- parameters ---
taxon_id <- 195
string_version <- "12.0"
out_dir <- "data/generated/string_195_v12"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- URLs / paths ---
links_url <- sprintf(
  "https://stringdb-downloads.org/download/protein.links.v%s/%d.protein.links.v%s.txt.gz",
  string_version, taxon_id, string_version
)
links_gz <- file.path(out_dir, sprintf("%d.protein.links.v%s.txt.gz", taxon_id, string_version))

seq_url <- sprintf(
  "https://stringdb-downloads.org/download/protein.sequences.v%s/%d.protein.sequences.v%s.fa.gz",
  string_version, taxon_id, string_version
)
seq_gz <- file.path(out_dir, sprintf("%d.protein.sequences.v%s.fa.gz", taxon_id, string_version))

# --- download ---
download.file(links_url, destfile = links_gz, mode = "wb", quiet = TRUE)
download.file(seq_url,   destfile = seq_gz,   mode = "wb", quiet = TRUE)

# --- read links ---
# Typical columns: protein1 protein2 combined_score (whitespace-separated)
dt <- fread(links_gz)

# Optional: filter by score if desired (uncomment and set threshold)
# score_threshold <- 700
# dt <- dt[combined_score >= score_threshold]

edges <- dt[, .(protein1, protein2)]

# --- build integer node IDs from proteins in the edge list ---
proteins <- sort(unique(c(edges$protein1, edges$protein2)))
node_id <- seq_along(proteins)
names(node_id) <- proteins

edge_list <- data.table(
  left  = unname(node_id[edges$protein1]),
  right = unname(node_id[edges$protein2])
)

# --- read FASTA and map protein -> sequence ---
read_fasta_gz <- function(path_gz) {
  # returns named character vector: names = IDs, values = sequences (no whitespace)
  con <- gzfile(path_gz, open = "rt")
  on.exit(close(con), add = TRUE)

  x <- readLines(con, warn = FALSE)
  if (length(x) == 0) stop("FASTA file is empty: ", path_gz)

  header_idx <- which(startsWith(x, ">"))
  if (length(header_idx) == 0) stop("No FASTA headers found in: ", path_gz)

  header_end <- c(header_idx[-1] - 1, length(x))

  ids <- character(length(header_idx))
  seqs <- character(length(header_idx))

  for (i in seq_along(header_idx)) {
    h <- x[header_idx[i]]
    # Take the first token after '>' as the sequence ID
    ids[i] <- sub("^>(\\S+).*", "\\1", h)

    if (header_idx[i] == header_end[i]) {
      seqs[i] <- ""  # header with empty sequence (shouldn't happen, but handle)
    } else {
      s <- x[(header_idx[i] + 1):header_end[i]]
      s <- gsub("\\s+", "", s)  # remove whitespace/newlines just in case
      seqs[i] <- paste0(s, collapse = "")
    }
  }

  # If duplicates exist, keep the first (rare; depends on source)
  keep <- !duplicated(ids)
  seqs <- seqs[keep]
  ids  <- ids[keep]

  setNames(seqs, ids)
}

seq_map <- read_fasta_gz(seq_gz)

# --- create node_labels using sequences as labels ---
# Match by STRING protein ID. Some IDs in edges might be missing in FASTA (or vice versa).
seq_for_nodes <- unname(seq_map[proteins])

missing <- which(is.na(seq_for_nodes) | seq_for_nodes == "")
if (length(missing) > 0) {
  warning(sprintf(
    "Missing sequences for %d/%d nodes. Example missing IDs: %s",
    length(missing), length(proteins),
    paste(head(proteins[missing], 10), collapse = ", ")
  ))
  # Decide what you want: stop, or keep NA/empty. Here we stop to enforce correctness.
  stop("Cannot write node_labels.csv because some node sequences are missing.")
}

node_labels <- data.table(
  left  = node_id,
  right = seq_for_nodes
)

# --- write outputs in your required simple format (semicolon-separated) ---
fwrite(edge_list,
       file = file.path(out_dir, "edge_list.csv"),
       sep = ";", col.names = TRUE)

# Important: sequences can be long; write quoted to be safe CSV-wise
fwrite(node_labels,
       file = file.path(out_dir, "node_labels.csv"),
       sep = ";", col.names = TRUE, quote = TRUE)

# --- sanity checks ---
stopifnot(nrow(node_labels) == length(proteins))