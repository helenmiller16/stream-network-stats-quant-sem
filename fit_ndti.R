fit_spacetime = FALSE
fit_spatial = FALSE
fit_temporal = TRUE

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


use_id <- 44240200011
source("C:/Users/helen/OneDrive/Documents/mekong/tsl-nutrients/03_figures/themes.R")

# load data ----- 
## reflectance data -----
file_dir <- ("C:/github/lmb-metabolism-sensing/data/output/daily_medians")
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
data <- data[ !(variable != "count" & value > 1) ] 

# make wide
wide = dcast(data, time + reach_id ~ variable, value.var = 'value')
wide <- wide[complete.cases(wide)]

wide[, qc_percent := count / max(.SD$count), reach_id]
wide[, ndti := (red - green)/(red + green)]
setnames(wide, "time", "date")

## spatial data ----- 

if (!file.exists("river_geom_modeled.geojson")) {
  river_geom <- st_read("C:/github/lmb-metabolism-sensing/data/rivers.geojson") |> 
    st_transform("EPSG:32648")
  river_geom_ <- river_geom[river_geom$reach_id %in% wide$reach_id, ]
  
  # add back in any sections of river which do not have any data 
  # to make sure they're part of the network
  fix_missing <- function(r_g) {
    missing <- r_g$rch_id_dn[!(r_g$rch_id_dn %in% r_g$reach_id) & 
                               (r_g$rch_id_dn %in% river_geom$reach_id)]
    if (length(missing) > 0) {
      r_g <- fix_missing(
        rbind(r_g, river_geom[river_geom$reach_id %in% missing, ]) 
      )
    } 
    return(r_g)
  }
  river_geom <- fix_missing(river_geom_)
  st_write(river_geom, "river_geom_modeled.geojson")
} else {
  river_geom <- st_read("river_geom_modeled.geojson")
}

network <- as_sfnetwork(river_geom) 

reflectance_data <- merge(river_geom['reach_id'], wide, by = "reach_id", all.x=TRUE, sort = FALSE)

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
# only model every 5 days... 
# and snap other days to nearest
days <- seq(as.IDate("2023-04-30"), as.IDate("2024-04-29"), by = 5)

data[!(data$date %in% days), ]$date <- 
  (data[!(data$date %in% days), ]$date - as.difftime(2, unit="days"))

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
  y_nt[x$from, which(days == x$date)] <- x[['ndti']]
}
# transform y to 0-1
# AVW can range from coastal to red wavelength
# y_nt <- (y_nt - min(nm_map)) / (max(nm_map) - min(nm_map))

library(TMB)

if (fit_temporal) {
  # Fit timeseries ----- 
  compile("time_only.cpp")
  dyn.load(dynlib("time_only"))
  
  y_t <- y_nt[edges[edges$reach_id == use_id, ]$from, ]
  
  sigma_x <- 1
  sigma_y <- 1
  alpha <- -.2
  beta <- 0.1
  
  x_t <- rnorm(length(y_t))
  
  Params <- list(
    logsigma_y = log(sigma_y), 
    logsigma_x = log(sigma_x),
    #alpha = alpha, 
    #logbeta = log(beta), 
    logit_rho = 0, 
    x_t = x_t
  )
  
  Data <- list(
    y_t = y_t
  )
  
  Obj <- MakeADFun(data = Data, 
                   parameters = Params,
                   random = c(
                     "x_t") ,
                   DLL = "time_only")
  Obj$env$beSilent()
  system.time({
    Opt <- nlminb( start = Obj$par, 
                   obj = Obj$fn, 
                   gr = Obj$gr, 
                   control=list(trace=1, eval.max=1e4, iter.max=1e4))
  })
  r <- Obj$report()
  Opt$SD <- sdreport(Obj, bias.correct = TRUE)
  
  rand <- summary(Opt$SD, "random")
  
  # plot estimates
  z_t <-  data.table(x = rand[, "Estimate"])
  z_t$time <- days
  z_t$se <- rand[, 'Std. Error']
  
  
  fwrite(z_t, paste0("modeling_output_time_", use_id, "_ndti.csv"))
  
}


# Prep Spatial data -----
theta <- .004 # autocorrelation parameter
alpha = -.2 # intercept
beta <- .01 # effect of x on ndti
sigma <- 1 # variance parameter for y
rho_t <- 0.9
beta_t <- 0.1

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

Params <- list(
  logsigma_y = log(sigma), 
  logtheta = log(theta), 
  logbeta = log(beta), 
  x_n = rep(0, N)
)

node_to_reachid <- data.table(node = edges$from, 
                              reach_id = edges$reach_id)
setorder(node_to_reachid, node)
space_data <- NULL


if (fit_spatial) {
  # Fit spatial -----
  compile( "network_tailup.cpp")
  dyn.load( dynlib("network_tailup") )
  
  
  for (date in as.IDate(c("2023-05-20", 
                          "2023-05-25", 
                          "2023-05-30", 
                          "2023-06-04")) ) {
    print("----------------------------------")
    print(date)
    
    d <- as.IDate(date)
    Data <- list(
      y_n = y_nt[, which(days == d)], 
      from_e = table$from -1, # index from 0
      to_e = table$to -1, # index from 0
      dist_e = table$length_km,
      flow_n = flow, 
      #hydro_h = which(nodes$hydro == 1)-1, # index from 0
      source_s = sources -1 # index from 0
    )
    Params <- list(
      logsigma_y = log(sigma), 
      alpha = 0, 
      logtheta = log(theta), 
      logbeta = log(beta), 
      x_n = rep(0, N)
    )
    
    Obj <- MakeADFun(data = Data, 
                     parameters = Params,
                     random = c("x_n") ,
                     DLL = "network_tailup")
    Obj$env$beSilent()
    system.time({
      Opt <- nlminb( start = Obj$par, 
                     obj = Obj$fn, 
                     gr = Obj$gr, 
                     control=list(trace=1, eval.max=1e4, iter.max=1e4))
    })
    r <- Obj$report()
    Opt$SD <- sdreport(Obj, bias.correct = TRUE)
    rand <- summary(Opt$SD, "random")
    dt <- data.table(
      date = d,
      node = 1:(N),
      reach_id = c(node_to_reachid$reach_id, NA),
      pred = rand[, "Estimate"], 
      se = rand[, "Std. Error"]
    )
    if (is.null(space_data)) {
      space_data <- dt
    } else {
      space_data <- rbind(space_data, dt)
    }
  }
  
  
  fwrite(space_data, "data_spatial_only_ndti.csv")
}



# Fit spacetime -----
if (fit_spacetime) {
  compile( "network_spacetime.cpp" )
  dyn.load( dynlib("network_spacetime") )
  
  
  
  omega_nt = matrix(rnorm(length(y_nt)), 
                    nrow = N, ncol = timesteps) |> 
    rbind(rep(1, timesteps)) # add row of 1's (hydro effect)
  
  # TRY fixing psi (also don't fit beta1)
  Params <- list(
    logtheta = log(theta), # autocorrelation parameter
    logsigma_y = log(sigma), # variance parameter for y
    alpha = 0, # intercept
    logbeta1 = log(beta), # effect of psi on ndti
    logbeta2 = log(beta_t), # effect of omega on ndti
    #logit_rhoW = logit(rho_t), # temporal autocorrelation parameter
    
    psi_n = rep(0, N), # spatial effect
    omega_nt = matrix(rnorm(length(y_nt)), 
                      nrow = N, ncol = timesteps)# spatio-temporal effect
  )
  
  
  
  Data <- list(
    #ndti_n = (nodes$ndti +1 )/2,
    n_t = timesteps,
    y_nt = y_nt,
    from_e = table$from -1, # index from 0
    to_e = table$to -1, # index from 0
    dist_e = table$length_km,
    flow_n = flow, 
    #hydro_h = which(nodes$hydro == 1)-1, # index from 0
    source_s = sources -1 # index from 0
  )
  
  
  Obj <- MakeADFun(data = Data, 
                   parameters = Params,
                   random = c(
                     "psi_n", 
                     "omega_nt") ,
                   # map = list(
                   #   psi_n = factor(rep(NA, N)), 
                   #   beta1 = factor(NA)
                   # ),
                   DLL = "network_spacetime")
  Obj$env$beSilent()
  system.time({
    Opt <- nlminb( start = Obj$par, 
                   obj = Obj$fn, 
                   gr = Obj$gr, 
                   control=list(trace=1, eval.max=1e4, iter.max=1e4))
  })
  r <- Obj$report()
  Opt$SD <- sdreport(Obj, bias.correct = TRUE)
  
  # Make figures of spatial covariance and precision
  # Precision
  Matrix::image(r$Q, 
                xlab = "", ylab = "", sub = "", yaxt = 'n', axis = \(...){}, main = "Precision")
  
  # Covariance
  Matrix::image(solve(r$I-r$Gamma) %*% solve(r$V) %*% solve(t(r$I - r$Gamma)), 
                xlab = "", ylab = "", sub = "", yaxt = 'n', axis = \(...){}, main = "Covariance")
  
  summary(Opt$SD, "fixed")
  rand <- summary(Opt$SD, "random")
  
  # plot estimates
  z_nt <- matrix(summary(Opt$SD, "report")[, 'Est. (bias.correct)'], 
                 nrow = N, ncol = timesteps)
  z_nt <- r$z_nt
  colnames(z_nt) <- days
  z_nt <- as.data.table(z_nt)
  z_nt$id <- 1:nrow(z_nt)
  z_dt <- melt(z_nt, id.vars = "id")
  z_dt$variable <- as.Date(as.numeric(as.character(z_dt$variable)))
  setnames(z_dt, c("variable", "value"), c("date", "pred"))
  z_dt$se <- summary(Opt$SD, "report")[, 'Std. Error']
  z_dt$omega <- rand[rownames(rand) == "omega_nt", "Estimate"]
  # add reach id
  z_dt <- merge(st_set_geometry(edges[c('from', 'reach_id')], NULL), z_dt, 
                all.y = TRUE, by.x = c('from'), by.y = c("id"))
  # add geom
  z_sf <- st_sf(z_dt, st_as_sfc(nodes[z_dt$from, ]))
  
  
  
  res <- merge(data, z_dt, all.y = TRUE, by.x = c("reach_id", "from", "date"), by.y = c("reach_id", "from", "date"))
  res$response <-  res$pred 
  res$se_response <- res$se 
  
  # save result
  st_set_geometry(res, NULL) |>
    fwrite('modeling_output_spacetime_ndti.csv')
  
  
}
