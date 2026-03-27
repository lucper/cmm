library(tidyverse)

input_id <- "aj-yeast"

# ── 1. Read data ──────────────────────────────────────────────
edges  <- read_delim(glue("./data/{input_id}/edge_list.csv"),  delim = ";", col_types = "ii")
labels <- read_delim(glue("./data/{input_id}/node_labels.csv"), delim = ";", col_types = "ic")

# ── 2. Node degree distribution ──────────────────────────────
# Each edge (left, right) contributes one degree to each endpoint.
degrees <- bind_rows(
  edges %>% select(node = left),
  edges %>% select(node = right)
) %>%
  count(node, name = "degree")

# Summary statistics
summary(degrees$degree)

# Histogram
ggplot(degrees, aes(x = degree)) +
  geom_histogram(binwidth = 1, fill = "steelblue", colour = "white") +
  labs(
    title = "Distribution of Node Degrees",
    x     = "Degree",
    y     = "Count"
  ) +
  theme_minimal()

# Log-log ECDF (useful if the degree distribution is heavy-tailed)
ggplot(degrees, aes(x = degree)) +
  stat_ecdf(geom = "step") +
  scale_x_log10() +
  scale_y_continuous(trans = "reverse", labels = scales::percent) +
  labs(
    title = "Complementary CDF of Node Degrees (log-log)",
    x     = "Degree (log scale)",
    y     = "P(Degree ≥ x)"
  ) +
  theme_minimal()

# ── 3. Node label length distribution ────────────────────────
label_lengths <- labels %>%
  mutate(label_length = nchar(right))

# Summary statistics
summary(label_lengths$label_length)

# Histogram
ggplot(label_lengths, aes(x = label_length)) +
  geom_histogram(bins = 50, fill = "darkorange", colour = "white") +
  labs(
    title = "Distribution of Node Label Lengths",
    x     = "Label Length (characters)",
    y     = "Count"
  ) +
  theme_minimal()