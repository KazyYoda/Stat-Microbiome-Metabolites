############################################################
# Geographic Distribution of Sampling Sites
#
# Purpose:
#   Visualize the geographic distribution of sampling sites
#   across Chiang Mai and Chiang Rai, Thailand.
#
# Outputs:
#   1. Thailand overview map with sampling sites
#   2. Zoomed map of Chiang Mai and Chiang Rai
#   3. Sampling-site map with terrain/hillshade
#   4. Sampling-site map with elevation
#   5. Combined elevation and hillshade map
############################################################


# ==========================================================
# 1. Install and Load Required Packages
# ==========================================================
# Install packages if they are not already available.
# This step only needs to be run once per R installation.

install.packages(c("sf", "ggplot2", "rnaturalearth", "rnaturalearthdata"))


library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)


# ==========================================================
# 2. Define Sampling Sites
# ==========================================================
# Create a table containing the sampling-site information.
#
# Lon and Lat:
#   Geographic coordinates in decimal degrees (WGS84).
#
# Site:
#   Unique identifier combining ethnicity and sampling area.
#
# Ethnicity:
#   Ethnic group represented at each sampling site.

sites <- data.frame(
  Site = c("Akha-CM",
           "Akha-CR",
           "Lahu-CM",
           "Lahu-CR",
           "Lisu-CR",
           "Khuen-CM"),
  Ethnicity = c("Akha","Akha",
                "Lahu","Lahu",
                "Lisu","Khuen"),
  Lon = c(99.3058,
          99.6800,
          99.2600,
          99.6500,
          99.6800,
          98.8794),
  Lat = c(19.4474,
          20.1500,
          19.4600,
          20.1200,
          20.1200,
          18.5533)
)


# Convert the sampling-site coordinates into an sf spatial object.
# CRS 4326 corresponds to WGS84 geographic coordinates.

sites_sf <- st_as_sf(
  sites,
  coords = c("Lon","Lat"),
  crs = 4326
)


# ==========================================================
# 3. Obtain Thailand Country Boundary
# ==========================================================
# Download the Thailand country boundary from Natural Earth.
# The resulting object is returned as an sf spatial object.

thailand <- ne_countries(
  country = "Thailand",
  scale = "medium",
  returnclass = "sf"
)


# ==========================================================
# 4. Thailand Overview Map
# ==========================================================
# Display all sampling sites on a map of Thailand.
#
# Sampling sites are represented by different shapes according
# to ethnic group.

ggplot() +
  geom_sf(data = thailand,
          fill = "grey95",
          color = "black") +
  
  geom_sf(data = sites_sf,
          aes(shape = Ethnicity),
          size = 3) +
  
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  coord_sf(
    xlim = c(97, 106),
    ylim = c(5, 21)
  )


# ==========================================================
# 5. Prepare Sampling-Site Labels
# ==========================================================
# Extract the coordinates from the sf object so that site names
# can be positioned using geom_text_repel().

library(ggrepel)

coords <- st_coordinates(sites_sf)

sites_label <- cbind(
  sites,
  X = coords[,1],
  Y = coords[,2]
)


# ==========================================================
# 6. Load and Prepare Province Boundaries
# ==========================================================
# Select Chiang Mai and Chiang Rai provinces for the detailed map.
#
# Note:
#   'tha_prov' must already exist in the R environment at this
#   point. It is generated later in this script using GADM data.

library(dplyr)

north_prov <- subset(
  tha_prov,
  Geographic_location %in% c("Chiang Mai", "Chiang Rai")
)


# ==========================================================
# 7. Zoomed Map of Sampling Sites
# ==========================================================
# Display the sampling sites within Chiang Mai and Chiang Rai.
# Site labels are automatically positioned using ggrepel to
# minimize overlap.

geom_sf_text(
  data = sites_sf,
  aes(label = Site),
  nudge_y = 0.08,
  size = 3
)


ggplot() +
  geom_sf(
    data = north_prov,
    aes(fill = Geographic_location),
    color = "black",
    linewidth = 0.4
  ) +
  
  geom_sf(
    data = sites_sf,
    aes(shape = Ethnicity),
    size = 3
  ) +
  
  geom_text_repel(
    data = sites_label,
    aes(
      x = X,
      y = Y,
      label = Site
    ),
    size = 4,
    max.overlaps = Inf
  ) +
  
  scale_fill_manual(
    values = c(
      "Chiang Mai" = "darkgoldenrod2",
      "Chiang Rai" = "deepskyblue3"
    )
  ) +
  
  coord_sf(
    xlim = c(98.5, 100.0),
    ylim = c(18.3, 20.4)
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude"
  )


# ==========================================================
# 8. Download and Prepare Elevation Data
# ==========================================================
# Install packages required to obtain and process elevation data.
#
# terra:
#   Raster and spatial data processing.
#
# elevatr:
#   Download elevation data for specified geographic locations.

install.packages(c("terra", "elevatr"))
library(terra)
library(elevatr)
library(sf)


# Select Chiang Mai and Chiang Rai from the province shapefile.

north_prov <- subset(
  tha_prov,
  NAME_1 %in% c("Chiang Mai", "Chiang Rai")
)


# Download a digital elevation model (DEM) for the study region.
#
# z = 8 controls the spatial resolution of the downloaded
# elevation data.

dem <- get_elev_raster(
  locations = north_prov,
  z = 8,
  clip = "locations"
)

dem <- rast(dem)


# ==========================================================
# 9. Calculate Terrain Variables
# ==========================================================
# Derive slope and aspect from the digital elevation model.
# These variables are used to generate a hillshade layer.

slope <- terrain(dem, "slope", unit = "radians")
aspect <- terrain(dem, "aspect", unit = "radians")

hill <- shade(slope, aspect)


# Convert the hillshade raster into a data frame for ggplot2.

hill_df <- as.data.frame(hill, xy = TRUE)
names(hill_df)[3] <- "hillshade"


# Use province names as geographic-location labels.

north_prov$Geographic_location <- north_prov$NAME_1


# ==========================================================
# 10. Sampling Sites with Terrain Hillshade
# ==========================================================
# Overlay the sampling sites and province boundaries on a
# hillshade background to provide geographic terrain context.

ggplot() +
  
  geom_raster(
    data = hill_df,
    aes(x = x, y = y, fill = hillshade)
  ) +
  scale_fill_gradient(
    low = "white",
    high = "grey40",
    guide = "none"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  geom_sf(
    data = north_prov,
    aes(fill = Geographic_location),
    color = "black",
    linewidth = 0.5,
    alpha = 0.35
  ) +
  
  scale_fill_manual(
    name = "Geographic location",
    values = c(
      "Chiang Mai" = "darkgoldenrod2",
      "Chiang Rai" = "deepskyblue3"
    )
  ) +
  
  geom_sf(
    data = sites_sf,
    aes(shape = Ethnicity),
    size = 3
  )   + 
  
theme_bw() +
  theme(
    panel.grid = element_blank()
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude"
    ) 


# ==========================================================
# 11. Obtain Thailand Province Boundaries
# ==========================================================
# Download first-level administrative boundaries (provinces)
# for Thailand using the GADM database.
#
# The resulting object is converted to an sf object for use
# with ggplot2 and other spatial operations.

install.packages("geodata")
library(geodata)
library(sf)


tha_prov <- geodata::gadm(
  country = "THA",
  level = 1,
  path = tempdir()
)


tha_prov <- st_as_sf(tha_prov)


# ==========================================================
# 12. Thailand Province Map
# ==========================================================
# Display all Thai provinces together with the sampling sites.

ggplot() +
  geom_sf(
    data = tha_prov,
    fill = "grey95",
    color = "grey40",
    linewidth = 0.2
  ) +
  geom_sf(
    data = sites_sf,
    aes(shape = Ethnicity),
    size = 3
  ) +
  coord_sf(
    xlim = c(97, 106),
    ylim = c(5, 21)
  ) +
  theme_bw()


# ==========================================================
# 13. Highlight Chiang Mai and Chiang Rai
# ==========================================================
# Create a categorical variable identifying the two study
# provinces while assigning all other provinces to "Other".

names(tha_prov)
tha_prov$Geographic_location <- "Other"

tha_prov$Geographic_location[tha_prov$NAME_1 == "Chiang Mai"] <- "Chiang Mai"
tha_prov$Geographic_location[tha_prov$NAME_1 == "Chiang Rai"] <- "Chiang Rai"


# Plot Thailand with Chiang Mai and Chiang Rai highlighted.

ggplot() +
  geom_sf(
    data = tha_prov,
    aes(fill = Geographic_location),
    color = "NA",
    linewidth = 0.1
  ) +
  geom_sf(
    data = sites_sf,
    aes(shape = Ethnicity),
    size = 3
  ) +
  scale_fill_manual(
    name = "Geographic location",
    values = c(
      "Chiang Mai" = "darkgoldenrod2",
      "Chiang Rai" = "deepskyblue3",
      "Other" = "grey85"
    )
  ) +
  coord_sf(
    xlim = c(97, 106),
    ylim = c(5, 21)
  ) +
  theme_bw() +
  
  labs(
    x = "Longitude",
    y = "Latitude"
  ) 


# ==========================================================
# 14. Sampling Sites with Elevation
# ==========================================================
# Visualize elevation across the study region and overlay
# province boundaries and sampling sites.

library(ggplot2)
library(ggnewscale)
library(viridis)


# Convert the DEM raster into a data frame for ggplot2.

elev_df <- as.data.frame(dem, xy = TRUE)
names(elev_df)[3] <- "Elevation"


north_prov$Geographic_location <- north_prov$NAME_1


ggplot() +
  
  # Elevation background
  geom_raster(
    data = elev_df,
    aes(x = x, y = y, fill = Elevation)
  ) +
  
  scale_fill_gradientn(
    colours = terrain.colors(10),
    name = "Elevation (m)"
  )+
  
  ggnewscale::new_scale_fill() +
  
  # Province overlay
  geom_sf(
    data = north_prov,
    aes(fill = Geographic_location),
    color = "black",
    linewidth = 0.5,
    alpha = 0.30
  ) +
  
  scale_fill_manual(
    name = "Geographic location",
    values = c(
      "Chiang Mai" = "darkgoldenrod2",
      "Chiang Rai" = "deepskyblue3"
    )
  ) +
  
  # Sampling sites
  geom_sf(
    data = sites_sf,
    aes(shape = Ethnicity),
    size = 3
  ) +
  
  coord_sf(
    xlim = c(98.5, 100.0),
    ylim = c(18.3, 20.4)
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude"
  )


# ==========================================================
# 15. Combined Elevation and Hillshade Map
# ==========================================================
# Combine the digital elevation model with hillshade to provide
# both elevation information and visual representation of terrain.
#
# The two raster layers use separate fill scales. ggnewscale is
# used to allow independent legends/scales for these layers.

# Convert the DEM raster into a data frame and remove missing
# elevation values.

elev_df <- as.data.frame(dem, xy = TRUE, na.rm = TRUE)
names(elev_df)[3] <- "Elevation"


ggplot() +
  
  # Hillshade background
  geom_raster(
    data = hill_df,
    aes(x = x, y = y, fill = hillshade)
  ) +
  
  scale_fill_gradientn(
    colours = terrain.colors(10),
    guide = "none"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  # Semi-transparent elevation layer
  geom_raster(
    data = elev_df,
    aes(x = x, y = y, fill = Elevation),
    alpha = 0.5
  ) +
  
  scale_fill_gradientn(
    colours = terrain.colors(10),
    name = "Elevation (m)"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  # Province overlay
  geom_sf(
    data = north_prov,
    aes(fill = Geographic_location),
    color = "black",
    linewidth = 0.5,
    alpha = 0.25
  ) +
  
  scale_fill_manual(
    name = "Geographic location",
    values = c(
      "Chiang Mai" = "darkgoldenrod2",
      "Chiang Rai" = "deepskyblue3"
    )
  ) +
  
  # Sampling sites
  geom_sf(
    data = sites_sf,
    aes(shape = Ethnicity),
    size = 3
  ) +
  
  coord_sf(
    xlim = c(98.5, 100.0),
    ylim = c(18.3, 20.4)
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank()
  ) +
  
  labs(
    x = "Longitude",
    y = "Latitude"
  )
