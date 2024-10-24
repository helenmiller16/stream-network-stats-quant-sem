library(sf)
library(ggplot2)
library(data.table)

use_id <- 44230000861
source("C:/Users/helen/OneDrive/Documents/mekong/tsl-nutrients/03_figures/themes.R")

file_dir <- here::here("data/output/daily_medians")
files <- list.files(file_dir)
files <- files[grepl("daily_medians_sentinel_100mbuffer_",  files)]

r_ids <- sub('daily_medians_sentinel_100mbuffer_', '', files) 
r_ids <- sub('\\.csv', '', r_ids) |> 
  as.numeric()

data_list <- lapply(files, \(f) {
  id <- sub('daily_medians_sentinel_100mbuffer_', '', f) 
  id <- sub('\\.csv', '', id)
  df <- fread(file.path(file_dir, f))
  if (nrow(df) == 0) return(NULL)
  df$reach_id <- id
  setnames(df, names(df)[grepl("stackstac", names(df))], "value")
  df
})
data <- rbindlist(data_list, use.names = TRUE)
data[band != "count", value := ifelse( (value < 1) , value + 0.1, value / 10000)] 

data <- data[, .(
  time, 
  variable = band, 
  value, 
  reach_id
)]

# very basic QC
data <- data[ !(variable != "count" & value > 0.3) ] 


# make wide
wide = dcast(data, time + reach_id ~ variable, value.var = 'value')
wide <- wide[complete.cases(wide)]

wide[, reach_id := as.numeric(reach_id)]

wide[, ndvi := (nir - red)/(nir + red)]
wide[, ndci := (rededge1 - red)/(rededge1 + red)]
wide[, ndti := (red - green)/(red + green)]

# After
data_after <- fread(('modeling_output_spacetime_ndti.csv'))
data_spatial <- fread("data_spatial_only_ndti.csv")

# Make a map of data at first three days

## spatial data ----- 
river_geom <- st_read("river_geom_modeled.geojson") 
pred_space_sf <- merge(river_geom, data_spatial,by = "reach_id", all.x = TRUE, all.y = TRUE)
data_sf <- merge(river_geom, wide, by = "reach_id", all.x = TRUE, all.y = TRUE)
pred_sf <- merge(river_geom, data_after, by = "reach_id", all.x = TRUE, all.y = TRUE)

## Make maps ----- 
library(basemaps)

# convert everything to web mercator to work with basemaps
data_sf <- st_transform(data_sf, "EPSG:3857")
pred_sf <- st_transform(pred_sf, "EPSG:3857")
pred_space_sf <- st_transform(pred_space_sf, "EPSG:3857")
river_geom <- st_transform(river_geom, "EPSG:3857")

basemap_ggplot(map_service = "esri", 
               map_type = "world_hillshade",
               ext = st_buffer(data_sf, 10000), 
) + 
  geom_sf(data = river_geom, 
          color = "gray30") + 
  tsl_theme + 
  theme(axis.title.x = element_blank(),  # Remove x-axis label
        axis.title.y = element_blank(),  # Remove y-axis label
        axis.text.x = element_blank(),   # Remove x-axis tick labels
        axis.text.y = element_blank(),   # Remove y-axis tick labels
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.ticks.y = element_blank(),   # Remove y-axis ticks
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))
  

png("ndti_data_map.png", width = 8, height = 3, units = "in", res = 240)
basemap_ggplot(map_service = "esri", 
               map_type = "world_hillshade",
               ext = st_buffer(data_sf, 10000), 
               ) + 
  geom_sf(data = river_geom, 
          color = "gray30")+
  geom_sf(data = data_sf[data_sf$time %in% as.Date(c("2023-05-20", 
                                                     "2023-05-25", 
                                                     "2023-05-30", 
                                                     "2023-06-04")),], 
          aes(col = ndti), lwd = 2) + 
  scale_color_gradientn(colors = c(
    "#143123",
    "#7c985e", 
    "#f3bf5f", 
    "#fcf95a"
  ), 
  limits = c(-.3, 0.4), 
  name = "NDTI") + 
  #expand_limits(x = 0, y = 0) + 
  tsl_theme + 
  # scale_color_gradientn(colors = c(low = "#030149",
  #                                  mid = "#b2019a",
  #                                  high = "#fcdd5f"),
  #                       limits = c(-.3, 0.4)
  #                       ) +
  facet_grid(cols = vars(time)) +
  theme(axis.title.x = element_blank(),  # Remove x-axis label
        axis.title.y = element_blank(),  # Remove y-axis label
        axis.text.x = element_blank(),   # Remove x-axis tick labels
        axis.text.y = element_blank(),   # Remove y-axis tick labels
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.ticks.y = element_blank(),   # Remove y-axis ticks
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))
dev.off()


png("ndti_preds_spatial_map.png", width = 8, height = 3, units = "in", res = 240)
basemap_ggplot(map_service = "esri", 
               map_type = "world_hillshade",
               ext = st_buffer(data_sf, 10000)) + 
  geom_sf(data = pred_space_sf, 
          aes(col = pred), lwd = 2) + 
  tsl_theme + 
  scale_color_gradientn(colors = c("#143123",
                                   "#7c985e", 
                                   "#f3bf5f", 
                                   "#fcf95a"
  ), 
  limits = c(-.3, 0.4),
  oob = scales::oob_squish) +
  theme(axis.title.x = element_blank(),  # Remove x-axis label
        axis.title.y = element_blank(),  # Remove y-axis label
        axis.text.x = element_blank(),   # Remove x-axis tick labels
        axis.text.y = element_blank(),   # Remove y-axis tick labels
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.ticks.y = element_blank(),   # Remove y-axis ticks
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank()
  ) +
  facet_grid(cols = vars(date)) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  guides(color = "none")
dev.off()



png("ndti_preds_spacetime_map.png", width = 8, height = 3, units = "in", res = 240)
basemap_ggplot(map_service = "esri", 
               map_type = "world_hillshade",
               ext = st_buffer(data_sf, 10000), 
) + 
  geom_sf(data = pred_sf[pred_sf$date %in% as.Date(c("2023-05-20", 
                                              "2023-05-25", 
                                              "2023-05-30", 
                                              "2023-06-04")),], 
          aes(col = pred), lwd = 2) + 
  tsl_theme + 
  scale_color_gradientn(colors = c(
    "#143123",
    "#7c985e", 
    "#f3bf5f", 
    "#fcf95a"
  ), 
                        limits = c(-.3, 0.4), 
                        oob = scales::oob_squish) +
  facet_grid(cols = vars(date))+
  theme(axis.title.x = element_blank(),  # Remove x-axis label
        axis.title.y = element_blank(),  # Remove y-axis label
        axis.text.x = element_blank(),   # Remove x-axis tick labels
        axis.text.y = element_blank(),   # Remove y-axis tick labels
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.ticks.y = element_blank(),   # Remove y-axis ticks
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))
dev.off()

# make a png for each day to put together into a gif
days <- unique(pred_sf$date)[order(unique(pred_sf$date))]
#for (i in seq_along(days)) {
for (i in seq_along(days)) {
  png(paste0("animation_pngs/ndti_preds_", i, ".png"), width = 400, height = 400)
  # p <- basemap_ggplot(map_service = "esri",
  #                map_type = "world_hillshade",
  #                ext = st_buffer(data_sf, 10000),
  # ) +
  p <- ggplot() +
    geom_sf(data = pred_sf[pred_sf$date == days[i],], 
            aes(col = pred), lwd = 2) + 
    tsl_theme + 
    scale_color_gradientn(colors = c(
      "#143123",
      "#7c985e", 
      "#f3bf5f", 
      "#fcf95a"
    ), 
    limits = c(-.3, 0.4), 
    oob = scales::oob_squish) +
    theme(axis.title.x = element_blank(),  # Remove x-axis label
          axis.title.y = element_blank(),  # Remove y-axis label
          axis.text.x = element_blank(),   # Remove x-axis tick labels
          axis.text.y = element_blank(),   # Remove y-axis tick labels
          axis.ticks.x = element_blank(),  # Remove x-axis ticks
          axis.ticks.y = element_blank(),   # Remove y-axis ticks
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank()
    ) +
    ggtitle(strftime(days[i], "%B")) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  print(p)
  dev.off()
}

# Get reach id closest to a lat/lon
#13.251633, 107.360310
point <- st_point(c(107.360310, 13.251633)) |>
  st_sfc(crs = 4326) |> 
  st_transform(crs = st_crs(pred_sf))
river_geom[which.min(st_distance(point, pred_sf)), ]$reach_id
