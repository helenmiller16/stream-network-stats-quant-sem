suppressPackageStartupMessages({
  library(sf)
  library(data.table)
  library(sfnetworks)
  library(tidygraph)
  library(igraph)
  library(Matrix)
})

# function from textbook
rmvnorm_prec <-
  function( mu, # estimated fixed and random effects
            prec, # estimated joint precision
            n.sims = 1) {
    
    require(Matrix)
    # Simulate values
    z0 = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    # Q = t(P) * L * t(L) * P
    L = Cholesky(prec, super=TRUE)
    # Calcualte t(P) * solve(t(L)) * z0 in two steps
    z = solve(L, z0, system = "Lt") # z = Lt^-1 * z
    z = solve(L, z, system = "Pt") # z = Pt    * z
    return(mu + as.matrix(z))
  }
## Helpers
logit <- function(p) {
  log(p/(1-p))
}
invlogit <- function(x) {
  1/(1+exp(-x))
}


use_id <- 44230000861
source("C:/Users/helen/OneDrive/Documents/mekong/tsl-nutrients/03_figures/themes.R")

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

# Only look at Nov 1 - Jan 1
# data <- data[(as.numeric(time) > as.numeric(as.IDate("2023-11-01"))) & 
#        (as.numeric(time) < as.numeric(as.IDate("2024-01-01")))]

# make wide
wide = dcast(data, time + reach_id ~ variable, value.var = 'value')

# Derived indices
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
setnames(wide, "time", "date")

river_geom <- st_read(here::here("data/rivers.geojson")) |> 
  st_transform("EPSG:32648")
river_geom <- river_geom[river_geom$reach_id %in% wide$reach_id, ]
network <- as_sfnetwork(river_geom) 

reflectance_data <- merge(river_geom['reach_id'], wide, by = "reach_id", all.x=FALSE, sort = FALSE)


# Load data
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



# png(paste0("data_avw.png"), width = 10, height = 10, units = "in", res = 120)
ggplot(data_to_plot[
  (as.numeric(data_to_plot$date) > as.numeric(as.IDate("2023-11-01"))) & 
    (as.numeric(data_to_plot$date) < as.numeric(as.IDate("2024-01-01"))), 
]) +
  geom_sf(data = edges, col = 'lightblue') +
  geom_sf(aes(color = .data[['avw']] )) +
  theme_minimal() +
  scale_color_gradientn(colours = c("#053e93",
                                    "#008b81",
                                    "#9cc845",
                                    "#fcfc71"),
                        na.value = "transparent",
                        name = 'avw') +
  facet_wrap(~date, ncol = 4)
# dev.off()

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
# transform y to 0-1
# AVW can range from coastal to red wavelength
# y_nt <- (y_nt - min(nm_map)) / (max(nm_map) - min(nm_map))


library(TMB)
compile( "predict_reflectance/network_spacetime.cpp" )
dyn.load( dynlib("predict_reflectance/network_spacetime") )

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

omega_nt = matrix(rnorm(length(y_nt)), 
                  nrow = N, ncol = timesteps) |> 
  rbind(rep(1, timesteps)) # add row of 1's (hydro effect)

# TRY fixing psi (also don't fit beta1)
Params <- list(
  logtheta = log(theta), # autocorrelation parameter
  logphi = log(sigma), # variance parameter for y
  alpha = 0, # intercept
  logbeta1 = log(beta), # effect of psi on ndti
  logbeta2 = log(beta_t), # effect of omega on ndti
  logit_rhoW = logit(rho_t), # temporal autocorrelation parameter
  
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

summary(Opt$SD, "fixed")
rand <- summary(Opt$SD, "random")

# plot estimates
z_nt <- matrix(summary(Opt$SD, "report")[, 'Est. (bias.correct)'], 
               nrow = N, ncol = timesteps)
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

# plot psi
psi <- rand[rownames(rand) == "psi_n", ]
nodes$psi <- psi
plot(st_geometry(edges))
plot(nodes['psi'], pch = 16, add = TRUE)


res <- merge(data, z_dt, all.y = TRUE, by.x = c("reach_id", "from", "date"), by.y = c("reach_id", "from", "date"))
res$response <-  res$pred 
res$se_response <- res$se 

plot(res[['red']], res$pred, pch = 16, cex = .5, 
     xlab = "obs", ylab = "predicted")
abline(0, 1)
segments(res[['red']], res$response-res$se_response, y1=res$response+res$se_response, 
         col = adjustcolor(1, alpha.f = .3))


# save result
st_set_geometry(res, NULL) |>
  fwrite(here::here('modeling_output_red.csv'))

