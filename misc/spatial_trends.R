# Load required libraries
library(tidyverse)
library(sf)
#devtools::install_github("UrbanInstitute/urbnmapr")
library(urbnmapr)
library(viridis)

setwd('./multinomial_GP')

dat <- readRDS('data.rds') #read in data 

dat <- dat %>% select(floodcty_id, duration) %>% unique() %>% mutate(fips = substr(floodcty_id, nchar(floodcty_id) - 4, nchar(floodcty_id)))

summary_stats <- dat %>% group_by(fips) %>% summarise(avg = mean(duration), counts = n())

summary_stats <- summary_stats %>%
  mutate(fips = str_pad(fips, width = 5, pad = "0"))

# Load county geometries (from Urban Institute)
counties <- get_urbn_map("counties", sf = TRUE)
states <- get_urbn_map("states", sf = TRUE)

# Join using FIPS code
county_floods_map <- counties %>%
  left_join(summary_stats, by = c("county_fips" = "fips"))

# Remove Alaska and Hawaii
county_floods_map <- county_floods_map %>%
  filter(!state_name %in% c("Alaska", "Hawaii"))

states <- states %>%
  filter(!state_name %in% c("Alaska", "Hawaii"))

# Frequency 
freq <- ggplot() +
  geom_sf(data = county_floods_map, aes(fill = counts), color = "white", size = 0.1) +  # county borders
  geom_sf(data = states, fill = NA, color = "black", size = 0.5) +  # state outlines
  scale_fill_viridis(
    option = "magma",
    name = "Frequency",
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

county_floods_map %>%
  group_by(state_name) %>%
  summarise(has_data = any(!is.na(counts))) %>%
  filter(!has_data) %>%
  pull(state_name)

# Average duration
dur <- ggplot() +
  geom_sf(data = county_floods_map, aes(fill = avg), color = "white", size = 0.1) +  # county borders
  geom_sf(data = states, fill = NA, color = "black", size = 0.5) +  # state outlines
  scale_fill_viridis(
    option = "magma",
    name = "Average duration",
    breaks = seq(0, 10, by = 2),
    limits = c(0, 10),
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

library(ggpubr)
ggarrange(freq, dur, labels = c("A", "B"), ncol = 2)
