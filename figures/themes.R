library(ggplot2)
# Theme and colors -----
location_colors <- c("#47824c", "#c99246", "#a77eb2")
light_colors <- c( "#ea8579", "#828da5", "#a8c49b", "#eabd79", "#9e6d7c" )
dark_colors <- c("#8c372e", "#458a93","#2c5619",  "#a3640d", "#7e5389" )
data_colors_light <- scale_color_discrete(type = light_colors) 
data_colors_dark <- scale_color_discrete(type = dark_colors) 
tsl_theme <- theme_minimal() + 
  theme(panel.border = element_rect(fill = NA, color = "gray80"), 
        text = element_text(size=11))
inundation_time_labels <- c("300" = "Edge and open water \n (inundation days >= 300)", 
                             "100" = "Floodplain \n (inundation days < 300)")
theme_colors = c("lightgreen" = "#d0e0d2", 
                 'medgreen' = "#a8c49b", 
                 'darkgreen' = "#47824c", 
                 'darkdarkgreen' = "#2c5619")

scale_color_gradient_tsl <- function(low = '#194212', 
                                     high = '#cff756',
                                     ...) {
  scale_colour_gradient(
    low = low,
    high = high,
    ...
  )
}

