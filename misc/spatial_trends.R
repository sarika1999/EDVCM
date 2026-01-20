# Load required libraries
library(tidyverse)
library(sf)
#devtools::install_github("UrbanInstitute/urbnmapr")
library(urbnmapr)
library(viridis)

# Example: if your data is in a CSV
# It should have columns like: fips, flood_count
floods <- read_csv("~/Desktop/panic/misc/flood_counts.csv")
floods <- floods %>%
  mutate(fips = str_pad(fips, width = 5, pad = "0"))

# Example structure if already in memory:
# floods <- tibble(
#   fips = c("06037", "06059", "48113"),  # LA, Orange (CA), Dallas (TX)
#   flood_count = c(5, 3, 10)
# )

# Load county geometries (from Urban Institute)
counties <- get_urbn_map("counties", sf = TRUE)
states <- get_urbn_map("states", sf = TRUE)

# Inspect to confirm FIPS column name
# names(counties)
# -> includes "county_fips"

# Join using FIPS code
county_floods_map <- counties %>%
  left_join(floods, by = c("county_fips" = "fips"))

# Filter out Alaska and Hawaii
county_floods_map <- county_floods_map %>%
  filter(!state_name %in% c("Alaska", "Hawaii"))

states <- states %>%
  filter(!state_name %in% c("Alaska", "Hawaii"))

# Plot with borders visible (even for NA)
ggplot() +
  geom_sf(data = county_floods_map, aes(fill = frequency), color = "white", size = 0.1) +  # county borders
  geom_sf(data = states, fill = NA, color = "black", size = 0.4) +  # state outlines
  scale_fill_viridis(
    option = "plasma",
    name = "Frequency",
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

county_floods_map %>%
  group_by(state_name) %>%
  summarise(has_data = any(!is.na(frequency))) %>%
  filter(!has_data) %>%
  pull(state_name)
