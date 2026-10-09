
# Load libraries 

library(rgbif)
library(purrr)
library(stringr)
library(tidyr)
library(data.table)
library(sf)
library(maps)
library(dplyr)
library(data.table)
library(lubridate)
library(slider)
library(here)
library(ggplot2)

#1. ##########

# Read data #
setwd("/Users/juliag.dealedo/ONE/Postdoc/ANTENNA/WP3/DATA/EBMS/antenna_eBM_20250328")
DATA = fread("ebms_count.csv")
ebms_visit = fread("ebms_visit.csv")
ebms_coord = fread("ebms_transect_coord.csv")

#2. ##########
# Taxonomic standarization #
sp <- DATA %>% distinct(species_name) %>% filter(!is.na(species_name))
taxon_df <- name_backbone_checklist(name_data = sp$species_name)
taxon_df$species_name = taxon_df$verbatim_name

taxon_df_check = taxon_df %>%
  #filter(confidence < 90 | matchType != "EXACT") %>%
  mutate(
    species = if_else(is.na(species), species_name, species),
    taxonomic_source = if_else(is.na(species), "manual", "GBIF"))

DATA = DATA %>% left_join (taxon_df_check, by="species_name")


# Save species taxonomic information
fwrite(taxon_df_check, "ebms_taxon.csv")
rm (sp, taxon_df, taxon_df_check)

#3. ##########

# Join visit and count data to fill gaps in both sides #

DATA_simple = DATA %>% select(c("transect_id", "visit_date", "year", "species_name", "count"))
ebms_visit = ebms_visit %>% select(c("transect_id", "visit_date", "year"))

ebms_visit <- ebms_visit %>%
  mutate(transect_id = as.character(transect_id), visit_date = as.Date(visit_date))

DATA_simple <- DATA_simple %>%
  mutate(transect_id = as.character(transect_id), visit_date = as.Date(visit_date))

# Full dataset with the all visits: Those not found in counts.csv and those not found in visits.csv
data_full = full_join(ebms_visit, DATA_simple, by=c("transect_id", "visit_date", "year"))
# From now onwards we work on data_full
rm (DATA_simple, ebms_visit)


#4. ##########

# Filter transects have at least 10 consecutive years with a moving window. 
# Now ALL transects that have more than 10 consecutive years are retained. Even if there is a gap in between, but the 10 years still happen.
# We still need to decide on which timeseries to cut (e.g. 10 consecutive years from 2005-2015, gap in 2016, and other 10 years from 2017-2026)

transect_year = data_full %>% select(transect_id, year) %>% distinct()

  

transect_year_filter <- transect_year %>%
  distinct(transect_id, year) %>%
  arrange(transect_id, year) %>%
  group_by(transect_id) %>%
  summarise(
    consecutive_10 = any(slide_lgl(year, ~ length(.x) == 10 && all(diff(.x) == 1), .before = 9, .complete = TRUE)),
    n_years=n_distinct(year),
    .groups = "drop") %>% filter(consecutive_10)

valid_transects = (transect_year_filter$transect_id)

data_full_filtered = data_full %>% filter (transect_id %in% valid_transects) 
# filtered dataset without the short and not consecutive timeseries.

rm(data_full, transect_year, transect_year_filter)

#5. ##########

# Transects that have species with more than 50% of absences along the years need to be filtered out. 

visit_date_transect = data_full_filtered %>% select (visit_date, year, transect_id) %>% distinct()
transect_species = data_full_filtered %>% select (transect_id, species_name) %>% distinct()
visit_date_transect_species_com = transect_species %>% left_join(visit_date_transect, by = "transect_id")
visit_date_transect_species_count <- visit_date_transect_species_com %>% 
  left_join(data_full_filtered %>% 
              select(c("transect_id", "species_name", "year", "visit_date", "count")), 
            by = c("transect_id", "species_name", "year", "visit_date"))

visit_date_transect_species_count_grouped = visit_date_transect_species_count %>% 
  group_by(year,transect_id, species_name)  %>% 
  summarize(observed=!all(is.na(count)))

transect_species_valid_by_50 = visit_date_transect_species_count_grouped %>% 
  group_by(transect_id, species_name)  %>% 
  summarize(pct_true = mean(observed) * 100) %>%
  filter(pct_true >= 50.0) #filter by the 50% of presences
 
# Filter by the transects of interest
useful_transect_species= unique (paste0 (transect_species_valid_by_50$transect_id, transect_species_valid_by_50$species_name))
visit_date_transect_species_count$code_id = paste0(visit_date_transect_species_count$transect_id, visit_date_transect_species_count$species_name)
visit_date_transect_species_count_filtered = visit_date_transect_species_count %>% filter (code_id %in% useful_transect_species) 

rm (visit_date_transect, visit_date_transect_species_com, visit_date_transect_species_count, visit_date_transect_species_count_grouped)

#6. ##########

# Now, the transect-species-date combination needs to be filled with 0s. 
visit_date_transect_species_count_filtered_0 = visit_date_transect_species_count_filtered %>% mutate (count = ifelse (is.na(count), 0, count ))

# Date operations: preparing to average by month
visit_date_transect_species_count_filtered_0 = visit_date_transect_species_count_filtered_0 %>%  
  mutate(year_month = format(floor_date(visit_date, "month"), "%Y-%m"), month=month(visit_date))

# final number of models to run
length(unique(visit_date_transect_species_count_filtered_0$code_id))

# Mean number of counts per month
mean_all_months_count <- visit_date_transect_species_count_filtered_0 %>% 
  group_by(transect_id, species_name, year, month) %>% 
  summarize (mean_count = mean(count),
             sd_count = sd(count)) 

#7. ##########

# Infer 0 across non-fieldwork months (e.g. september-december)
# using data.frame to accelerate... 

setDT(mean_all_months_count)

months <- unique(mean_all_months_count[, .(transect_id, species_name, year)])[
  , .(month = 1:12), by = .(transect_id, species_name, year)]

mean_all_months_count <- merge(
  months,
  mean_all_months_count,
  by = c("transect_id", "species_name", "year", "month"),
  all.x = TRUE
)

mean_all_months_count[is.na(mean_count), mean_count := 0]
mean_all_months_count[is.na(sd_count), sd_count := 0]

#8. ##########

# Save df
fwrite(mean_all_months_count, "mean_all_months_count.csv")

mean_all_months_count %>%
  filter(transect_id == "ESBMS.223859") 

# Plot to check a time series
plot_example = mean_all_months_count %>%
  filter(species_name == "Pararge aegeria", transect_id == "ES-CTBMS.9") 

plot_example %>%  
  mutate(year_month = make_date(year, month, 1)) %>%
  ggplot(aes(x = year_month, y = mean_count)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  geom_point(size=3, alpha=0.5) +
  geom_area(alpha=0.2)+
  theme_minimal()+labs(
    title = plot_example$species_name,
    subtitle = plot_example$transect_id,
    x = "Year",
    y = "Mean count"
  ) 

here()
ggsave("example.png", width=18, height=9)

# 2. Definir las coordenadas originales (EPSG:3035)

mean_all_months_count %>% summarize(n_distinct(species_name), n_distinct(transect_id)) 
mean_all_months_count %>% select(c(species_name, transect_id)) %>% n_distinct()

datos = ebms_coord %>% filter (transect_id %in% mean_all_months_count$transect_id, !is.na(transect_lon)) 

puntos <- st_as_sf(datos, coords = c("transect_lon", "transect_lat"), crs = 3035)

# 3. Convertir a longitud y latitud en grados
xy <- st_coordinates(st_transform(puntos, crs = 4326))
library(ggplot2)
library(maps)
library(dplyr)

world <- map_data("world")
xy_df <- data.frame(lon = xy[, 1], lat = xy[, 2], puntos$transect_id)

ggplot() +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = "#fad5b2", color = "black", linewidth = 0.2) +
  geom_point(data = xy_df,
             aes(x = lon, y = lat),
             color = "steelblue", alpha = 0.3, size = 2) +
  coord_quickmap(xlim = c(-25, 45), ylim = c(27, 72), expand = FALSE) +
  labs(
    title = "Filtered Transect Locations",
    subtitle = "213 species · 2,538 transects · 47,074 time series",
    caption = "Species monitored for ≥10 consecutive years and detected in ≥50% of years",
    x = NULL, y = NULL
  ) + theme_void()

ggsave("filtered_transect_locations.png")


library(leaflet)

leaflet(xy_df) %>%
  addProviderTiles("Esri.WorldTopoMap") %>%
  addCircleMarkers(
    lng = ~lon,
    lat = ~lat,
    radius = 4,
    color = "steelblue",
    fillOpacity = 0.6,
    stroke = FALSE,
    label = ~puntos.transect_id
  ) %>%
  fitBounds(lng1 = -25, lat1 = 27, lng2 = 45, lat2 = 72)

