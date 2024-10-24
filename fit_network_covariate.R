suppressPackageStartupMessages({
  library(sf)
  library(data.table)
  library(sfnetworks)
  library(tidygraph)
  library(igraph)
  library(Matrix)
})

## Helpers
logit <- function(p) {
  log(p/(1-p))
}
invlogit <- function(x) {
  1/(1+exp(-x))
}


use_id <- 44230000861
source("C:/Users/helen/OneDrive/Documents/mekong/tsl-nutrients/03_figures/themes.R")

# load data ----- 
## reflectance data -----
file_dir <- here::here("data/external/outputs")
files <- list.files(file_dir)
files <- files[grepl("daily_medians_sentinel_100mbuffer_",  files)]

r_ids <- sub('daily_medians_sentinel_100mbuffer_', '', files) 
r_ids <- sub('\\.csv', '', r_ids) |> 
  as.numeric()

data_list <- lapply(files, \(f) {
  id <- sub('daily_medians_sentinel_100mbuffer_', '', f) 
  id <- sub('\\.csv', '', id)
  df <- fread(file.path(file_dir, f))
  df$reach_id <- id
  df
})
data <- rbindlist(data_list, use.names = TRUE)

data <- data[, .(
  time, 
  variable = band, 
  value = data[[names(data)[grepl("stackstac", names(data))]]] + 0.1, 
  reach_id
)]

# make wide
wide = dcast(data, time + reach_id ~ variable, value.var = 'value')

wide[, qc_percent := count / max(.SD$count), reach_id]
setnames(wide, "time", "date")

## spatial data ----- 
river_geom <- st_read(here::here("data/rivers.geojson")) |> 
  st_transform("EPSG:32648")
river_geom <- river_geom[river_geom$reach_id %in% wide$reach_id, ]
network <- as_sfnetwork(river_geom) 

reflectance_data <- merge(river_geom['reach_id'], wide, by = "reach_id", all.x=FALSE, sort = FALSE)

## Prepare to fit ----- 
lmb <- reflectance_data |>
  as_sfnetwork()
lmb <- st_transform(lmb, "EPSG:32648")


edges = st_as_sf(activate(network,"edges"))
nodes = st_as_sf(activate(network,"nodes"))

# get hydro variable
# 44240500011, 44240400011 are in the reservoir
# 44240300021 is downstream. 
edges$hydro <- 0
edges[edges$reach_id %in% c(44240300021), ]$hydro <- 1

nodes$hydro <- 0
nodes[edges$from, ]$hydro <- edges$hydro


data = st_as_sf(activate(lmb,"edges"))

st_geometry(data) <- st_geometry(nodes[data$from, ])

# remove rows with qa_percent less than 25%
data <- data[data$qc_percent > 0.25, ]

# then add rows to make sure we plot every date

days <- seq(min(data$date), max(data$date), by = 5)
x <- matrix(NA, nrow = length(days), ncol = ncol(data))
y <- data.frame(x)
colnames(y) <- names(data)
y$date <- days
y$geometry <- data$geometry[1:nrow(y)]
y <- st_as_sf(y)
data_to_plot <- rbind(data, y)

table = data.frame( from = edges$from,
                    to = edges$to,
                    # convert to km
                    length_km = units::drop_units(st_length(edges))/1000 )
N <- nrow(nodes)
timesteps <- length(days)


# make data table 
y_nt <- matrix(NA, nrow = N, ncol = timesteps)
for(row in 1:nrow(data)) {
  x <- data[row, ]
  y_nt[x$from, which(days == x$date)] <- x[['red']]
}



# get flow and sources
sources <- which(!1:nrow(nodes) %in% table$to)
# use width for weights
# get mean reach width for each source node
source_flow <- edges$width[sapply( sources, \(x) st_nearest_feature( x=nodes[x,], edges ))]

flow <- rep(0, N)
flow[sources] <- source_flow
for (source in sources) {
  # get list of all downstream points
  connected <- sapply(1:nrow(nodes), 
                      \(x) ifelse(source == x, 0, edge_connectivity(network, source, x)))
  # add source flow to all downstream points
  flow <- flow + unlist(connected)*flow[source]
}



library(TMB)
compile( "predict_reflectance/quant_seminar/network_spacetime_covariate.cpp" )
dyn.load( dynlib("predict_reflectance/quant_seminar/network_spacetime_covariate") )

theta <- .004 # autocorrelation parameter
alpha = -.2 # intercept
beta <- .01 # effect of x on ndti
sigma <- 1 # variance parameter for y
rho_t <- 0.9
beta_t <- 0.1

Params <- list(
  logtheta = log(theta), # autocorrelation parameter
  logsigma = log(sigma),
  alpha = alpha, # intercept
  logbeta1 = log(beta), # effect of x on red
  logbeta2 = log(beta),
  logit_rhoW = 1,
  h = 0, 
  psi_n = matrix(rnorm(nrow(y_nt))), 
  omega_nt = matrix(rnorm(length(y_nt)), 
                    nrow = N, ncol = timesteps)# spatio-temporal effect
)

Data <- list(
  #red_n = (nodes$red +1 )/2,
  n_t = timesteps,
  y_nt = y_nt, 
  from_e = edges$from -1, # index from 0
  to_e = edges$to -1, # index from 0
  dist_e = table$length_km,
  flow_n = flow, 
  hydro_n = nodes$hydro, # index from 0
  source_s = sources -1 # index from 0
)

Random <- "omega_nt"

Obj <- MakeADFun(data = Data, 
                 parameters = Params,
                 random = Random, 
                 DLL = "network_spacetime_covariate")


invisible(capture.output( {
  
  opt <- nlminb(start = Obj$par, 
                obj = Obj$fn,
                grad = Obj$gr)
  SD <- sdreport(Obj, bias.correct = TRUE)
  
}))

summary(SD, "fixed")
