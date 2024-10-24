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

nm_map = c(
  "coastal" = 443,
  "red"= 665, 
  "green" = 560, 
  "blue" = 490
)

wide[, avw := (coastal + red + green + blue) / 
    (coastal/nm_map['coastal'] + 
      red/nm_map['red'] + 
       green/nm_map['green'] + 
       blue/nm_map['blue'])
]

wide[, avw := (-7.4719643E-10* avw^5 + 
                 1.8794584E-06*avw^4 + 
                 -1.8924228E-03*avw^3 +
                 9.5069314E-01*avw^2 + 
                 -2.3623942E+02*avw + 
                 2.3384674E+04)
       ]
wide[, ndvi := (nir - red)/(nir + red)]
wide[, ndci := (rededge1 - red)/(rededge1 + red)]
wide[, ndti := (red - green)/(red + green)]
wide[, brightness := (red + green + blue + coastal)]

wide[, qc_percent := count / max(.SD$count), reach_id]
wide[, rgb := rgb(red, green, blue, maxColorValue = 0.3)]
wide[, niir := rgb(rededge1, nir08, red, maxColorValue = 0.3)]
wide[, nrg := rgb(rededge1, red, green, maxColorValue = 0.3)]
hsv <- rgb2hsv(wide$red, wide$green, wide$blue, maxColorValue = 0.5) |> t()
colnames(hsv) <- c("hue", "sat", "val")
wide <- cbind(wide, hsv)
# Need SWIR for FAI
# colors normalized
wide[, green_norm := green / (coastal + blue + green + red + nir)]
wide[, red_norm := red / (coastal + blue + green + red + nir)]
wide[, blue_norm := blue / (coastal + blue + green + red + nir)]
wide[, coastal_norm := coastal / (coastal + blue + green + red + nir)]
wide[, nir_norm := nir / (coastal + blue + green + red + nir)]
wide[, rededge1_norm := rededge1 / (coastal + blue + green + red + nir)]
wide[, turb := (289.29*pi*red)/(1-pi*red/.1641)]


# After

use_id <- 44240200011
data_after <- fread(('modeling_output_spacetime_ndti.csv'))
data_after_time <- fread(paste0('modeling_output_time_', use_id, '_ndti.csv'))
data_after[, reach_id := as.numeric(reach_id)]

 # id = 44230000861
# #id = 44240500111
# data_after[, reach_id := as.numeric(reach_id)]
# reach_after <- data_after[reach_id == id]
# reach <- wide[reach_id == id]
# 
# # use_id = 44230000861
# use_id = c(44240500111, 44240600011, 44240700011)
# use_id <- 44240200011
# 
# use_ids = c(#44240100011,
#             44240300011,
#             44240200011)
# 
# use_ids = c(44250000011, 44240200011)
# 
# 
# 
# use_id <- 44240200011
# Before
ggplot() + 
  geom_point(data = wide[reach_id %in% use_id & qc_percent > 0.25], 
             mapping = aes(
               x = time, 
               y = ndti
             ), 
             color = "#383024",
             #color = wide[reach_id %in% use_id & qc_percent > 0.25][['rgb']], 
             size = 2) + 
  # geom_line(data = wide[reach_id %in% use_id & qc_percent > 0.25], 
  #           mapping = aes(x = time, y = red, group = reach_id), color = "#8C372E")+
  tsl_theme +
  labs(x = "",
       y = "NDCI",
       title = "Turbidity index: water year 2023") +
  lims(
    x = c(min(wide$time), max(wide$time)), 
    y = c(-0.3, 0.3)
  ) +
  data_colors_dark

ggsave("timeseries_data.png")

# After
# time only
ggplot() + 
  geom_ribbon(data = data_after_time, 
              mapping = aes(x = time, 
                            ymin = (x - se), 
                            ymax = (x + se)), 
              alpha = 0.3, 
              fill = "#383024") + 
  geom_point(data = wide[reach_id %in% use_id & qc_percent > 0.25], 
             mapping = aes(
               x = time, 
               y = ndti
             ), 
             color = "#383024",
             #color = wide[reach_id %in% use_ids & qc_percent > 0.25][['rgb']], 
             size = 2) + 
  geom_line(data = data_after_time, 
            mapping = aes(x = time, y = x), color = "#383024")+
  
  tsl_theme +
  lims(y = c(-0.3, 0.3), 
       x = c(min(wide$time), max(wide$time)), )+
  labs(x = "", 
       y = "NDTI", 
       title = "Predicted Turdibity index (time only model)")
# ggsave("timeseries_temporal_interpolation.png")

data_after$reach_id <- as.character(data_after$reach_id)
# spacetime
ggplot() + 
  geom_ribbon(data=data_after[reach_id %in% use_id],
              mapping = aes(x = date,
                            ymin = response - se_response,
                            ymax = response + se_response,
                            group = reach_id, 
                            fill = reach_id),
              alpha = 0.3, 
              fill = "#383024"
              ) +
  # geom_line(data = wide[reach_id %in% use_id & qc_percent > 0.25],
  #           mapping = aes(x = time, y = red, group = reach_id), color = "#8C372E")+
  geom_point(data = wide[reach_id %in% use_id & qc_percent > 0.25], 
             mapping = aes(
               x = time, 
               y = ndti, 
               color = reach_id
             ), 
             color = "#383024",
             #color = wide[reach_id %in% use_ids & qc_percent > 0.25][['rgb']], 
             size = 2) + 
  geom_line(data =  data_after[reach_id %in% use_id], 
            mapping = aes(x = date, y = response, color = reach_id, group = reach_id), 
            color = "#383024"
            )+
  tsl_theme +
  lims(y = c(-0.3, 0.3), 
       x = c(min(wide$time), max(wide$time)))+
  labs(x = "", 
       y = "NDTI",
       title = "Predicted Turbidity index (spatio-temporal model)") +
  data_colors_dark
ggsave("timeseries_spatial_interpolation.png")

ggplot() + 
  geom_ribbon(data=data_after[reach_id == use_id], 
              mapping = aes(x = date, 
                  ymin = response - 2*se_response, 
                  ymax = response + 2*se_response),
              fill = 'lightblue', 
              alpha = 0.4) + 
  
  geom_point(data = wide[reach_id == use_id & qc_percent > 0.25], 
             mapping = aes(x = time, y = avw), 
             color = 'darkblue') + 
  geom_line(data = data_after[reach_id == use_id], 
            mapping = aes(x = date, y = response), 
            color = 'darkblue') + 
  geom_point(data = data_after[reach_id == use_id], 
            mapping = aes(x = date, y = response), 
            color = 'darkblue', alpha = 0.2) + 

# 
#   geom_point(data = wide[reach_id == id-20 & qc_percent > 0.25],
#              mapping = aes(x = time, y = avw),
#              color = "yellow") +
#   geom_line(data = data_after[reach_id == id-20],
#             mapping = aes(x = date, y = response),
#             color = "yellow") +
# 
#   geom_point(data = wide[reach_id == id-30 & qc_percent > 0.25],
#              mapping = aes(x = time, y = avw),
#              color = "orange") +
#   geom_line(data = data_after[reach_id == id-30],
#             mapping = aes(x = date, y = response),
#             color = "orange") +
# 
#   geom_point(data = wide[reach_id == 44240500021 & qc_percent > 0.25],
#              mapping = aes(x = time, y = avw),
#              color = "red") +
#   geom_line(data = data_after[reach_id == 44240500021],
#             mapping = aes(x = date, y = response),
#             color = "red") +
#   # 
#   geom_point(data = wide[reach_id == 44240200021 & qc_percent > 0.25],
#              mapping = aes(x = time, y = avw),
#              color = "darkred") +
#   geom_line(data = data_after[reach_id == 44240200021],
#             mapping = aes(x = date, y = response),
#             color = "darkred") +
  # 
  
  
  tsl_theme +
  labs(x = "", 
       y = "AVW",
       title = "Water color water year 2023") 


 # some other fun plots

ggplot(wide[qc_percent > 0.25]) + 
  geom_point(aes(x = sat,  y = ndci), 
             color = wide[qc_percent > 0.25 ][['rgb']], 
             size = 5) + 
  #geom_path(aes(x = ndti,  y = ndci, group = reach_id)) +  
  facet_wrap(~(ndti > -.025)) +
  tsl_theme  

ggplot(wide[qc_percent > 0.25]) + 
  geom_point(aes(x = green_norm,  y = ndci), 
             color = wide[qc_percent > 0.25 ][['rgb']], 
             size = 5) + 
  #geom_path(aes(x = ndti,  y = ndci, group = reach_id)) + 
  tsl_theme  

ggplot(wide[ qc_percent > 0.2 & reach_id == use_id]) + 
  geom_line(data = wide[reach_id == use_id], 
            mapping = aes(x = time,  y = qc_percent/2), alpha = 0.2) + 
  geom_line(aes(x = time,  y = red_norm), 
             color = 'orange') + 
  geom_line(aes(x = time,  y = green_norm), 
             color = 'darkgreen') + 
  geom_line(aes(x = time,  y = rededge1_norm), 
             color = 'red') + 
  geom_line(aes(x = time,  y = nir_norm), 
             color = 'darkred') + 
  geom_line(aes(x = time,  y = blue_norm), 
             color = 'blue') + 
  geom_line(aes(x = time,  y = coastal_norm), 
             color = 'purple') + 
 # geom_line(aes(x = time,  y = rededge1)) +
  #geom_path(aes(x = time,  y = red_norm, group = reach_id)) + 
  tsl_theme

ggplot(wide[reach_id == use_id]) + 
  
  geom_line(aes(x = time,  y = qc_percent)) + 
  geom_line(aes(x = time, y = ndci)) + 
  # geom_line(aes(x = time,  y = rededge1)) +
  #geom_path(aes(x = time,  y = red_norm, group = reach_id)) + 
  tsl_theme

ggplot(wide[qc_percent > 0.9 & reach_id == use_id]) + 
  
  geom_line(aes(x = time,  y = red_norm), 
            color = 'orange') + 
  # geom_line(aes(x = time,  y = rededge1)) +
  #geom_path(aes(x = time,  y = red_norm, group = reach_id)) + 
  tsl_theme


