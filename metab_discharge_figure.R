library(data.table)
library(ggplot2)
skg <- fread("C:/minidots/SKG2-END.csv")
skg_metab <- fread("C:/github/metabolism-3s/data/generated/preds_k600_eq_1_skg.csv")
skg_metab[, GPP := ifelse(GPP < 0, 0, GPP)]

# Now load MRC data
water_level <- fread("C:/github/mekong-remote-sensing/data/mrc20240117042527/Water Level.Telemetry@KH_014003_[Koh Key].csv")

# And predicted NDCI
ndci <- fread("C:/github/lmb-metabolism-sensing/predict_reflectance/quant_seminar/modeling_output_spacetime_ndci.csv")
ndci[, reach_id := as.character(reach_id)]
ndci <- ndci[reach_id == "44240200011", ]

ndti <- fread("C:/github/lmb-metabolism-sensing/predict_reflectance/quant_seminar/modeling_output_spacetime_ndti.csv")
ndti[, reach_id := as.character(reach_id)]
ndti <- ndti[reach_id == "44240200011", ]

# Convert timestamp to datetime
skg[,datetime := as.POSIXct(time, tz="Asia/Phnom_Penh")]
water_level[, datetime := as.POSIXct(`Timestamp (UTC+07:00)`, 
                                     format = "%Y-%m-%dT%H:%MZ", 
                                     tz = "Asia/Phnom_Penh")]
skg_metab[, date := as.POSIXct(date, tz = "Asia/Phnom_Penh")]
ndci[, date := as.POSIXct(date, tz = "Asia/Phnom_Penh")]
ndti[, date := as.POSIXct(date, tz = "Asia/Phnom_Penh")]


# Trim down dates
o_trimmed <- skg[datetime > as.POSIXct("2023-05-17 00:00:00", tz = "Asia/Phnom_Penh") & 
                    datetime < as.POSIXct("2023-11-19 00:00:00", tz = "Asia/Phnom_Penh"),]
water_trimmed <- water_level[datetime > as.POSIXct("2023-05-17 00:00:00", tz = "Asia/Phnom_Penh") & 
                               datetime < as.POSIXct("2023-11-19 00:00:00", tz = "Asia/Phnom_Penh"),]
gpp_trimmed <- skg_metab[date> as.POSIXct("2023-05-17 00:00:00", tz = "Asia/Phnom_Penh") &
                           date < as.POSIXct("2023-11-19 00:00:00", tz = "Asia/Phnom_Penh"), ]
ndci_trimmed <- ndci[date> as.POSIXct("2023-05-17 00:00:00", tz = "Asia/Phnom_Penh") &
                       date < as.POSIXct("2023-11-19 00:00:00", tz = "Asia/Phnom_Penh"), ]
ndti_trimmed <- ndti[date> as.POSIXct("2023-05-17 00:00:00", tz = "Asia/Phnom_Penh") &
                       date < as.POSIXct("2023-11-19 00:00:00", tz = "Asia/Phnom_Penh"), ]

# constants to translate water data to be around the same range as O data
w_mult <- 1.8
w_add <- 4

png("skg_DO_vs_ndti_preds.png", width = 1000, height = 500, type = "windows", res = 120, bg = "transparent")
ggplot() + 
  # geom_area(data = water_trimmed, aes(x = datetime, y= w_mult*Value+w_add), 
  #           fill = "#208dcc60", 
  #           color = "#208dcc") + 
  # geom_ribbon(data = water_trimmed,
  #             aes(x = water_trimmed$datetime,
  #                 ymin = w_add,
  #                 ymax = w_mult*Value+w_add),
  #             fill = "#208dcc20",
  #             color = "#208dcc40") +
  geom_line(data = o_trimmed, aes(x = datetime, y = DO), alpha = 0.6)+
  # geom_point(data = o_trimmed[datetime %in% as.POSIXct(c("2023-05-20 10:07:00", 
  #                                                        "2023-08-13 10:07:00"),
  #                                                      tz = "Asia/Phnom_Penh")],
  #            aes(x = datetime, y = DO),
  #            color = "#ffab40ff",
  #            size = 3)+
  # geom_line(data = gpp_trimmed,
  #           aes(x = date, y = -ER/3 + 4)) +
  geom_ribbon(data=ndti_trimmed,
              mapping = aes(x = date,
                            ymin = (response+.15)*30+5 - se_response*50,
                            ymax = (response+.15)*30+5 + se_response*50,),
              alpha = 0.2, 
              fill = "#383024"
  ) +
  geom_line(data = ndti_trimmed,
            aes(x = date, y = (pred+.15)*30+5),
            color = "#383024",
            linewidth = 1.2) +
  
  # geom_line(data = ndti_trimmed,
  #            aes(x = date, y = (pred + .1) * 100 + 5),
  #            color = "darkgreen") +
  geom_point(data = ndti_trimmed,
             aes(x = date, y = (ndti + .15) * 30 + 5),
             color = "#383024", 
             size = 2) +
  theme_minimal() + 
  labs(x = "", y = "DO (mg/L)") +
  theme(
    panel.background = element_rect(fill='transparent', color = NA), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    #panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'), #transparent legend panel
    panel.border = element_rect(fill = "transparent", colour = "gray50"), 
    text = element_text(size = 12)
  ) + 
  scale_x_datetime(expand = c(0,0)) + 
  scale_y_continuous(name = "DO (mg/L))",
                     limits =  c(2, 18),
                     sec.axis = sec_axis(transform = ~(. - 5)/50 - .15, name = "NDTI")) 

dev.off()
