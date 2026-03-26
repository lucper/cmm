library(tidyverse)
library(scales)
library(tools)
library(glue)


# Define the file path (replace with your actual file)
file_path <- glue("./logs/head10000_8_3.log")
log_id <- file_path_sans_ext(basename(file_path))
figure_dir <- "./data/figures/"

dir.create(figure_dir, recursive =  TRUE)

# Ingest and parse the data
parsed_data <- read_lines(file_path) %>%
  tibble(raw_text = .) %>%
  # Extract components matching the pattern: (Seq1, Seq2, Score)
  extract(
    col = raw_text,
    into = c("Sequence_1", "Sequence_2", "Score"),
    regex = "^\\(([^,]+),\\s*([^,]+),\\s*(\\d+)\\)$",
    convert = TRUE # Automatically coerces 'Score' to an integer
  )

head(parsed_data)

# Visualize the distribution of the 'Score' column
distribution_plot <- parsed_data %>%
  ggplot(aes(x = Score)) +
  geom_histogram(binwidth = 1) +
  geom_density(
    aes(y = after_stat(count)),
    adjust = 20,
    color = "red",
    linewidth = 1,
    alpha = 0.3,
  ) +
  labs(
    title = "Distribution of Motif Mining Scores",
    x = "Score",
    y = "Frequency"
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10),
    breaks = c(0, 1, 10, 100, 1000) # Explicitly define breaks to maintain readability
  )

print(distribution_plot)

ggsave(glue("{figure_dir}/distribution_{log_id}.pdf"))

