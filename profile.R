# model profiling
# library(lineprof)
# devtools::install_github("ncsu-landscape-dynamics/rpops", ref = "spatial-index")
library(PoPS)
library(raster)
library(lubridate)

folderfun::setff("In", "C:/Users/cmjone25/Desktop/Projects/pops_profiling/data/")

parameter_means <- c(1, 675, 0.99, 10000, 0, 0)
parameter_cov_matrix <- matrix(ncol = 6, nrow = 6, 0)
params_to_estimate <- c(T, T, T, T, F, F)
infected_file <- ffIn("int_infected.tif")
host_file <- ffIn("host100.tif")
total_populations_file <- ffIn("total_population.tif")
temp <- FALSE
temperature_coefficient_file <- ""
precip <- FALSE
precipitation_coefficient_file <- ""
model_type = "SI"
latency_period = 3
time_step <- "week"
season_month_start <- 1
season_month_end <- 12
start_date <- "2019-01-01"
end_date <- "2019-12-31"
use_lethal_temperature <- FALSE
temperature_file <- ""
lethal_temperature <- -30
lethal_temperature_month <- 1
mortality_on <- FALSE
mortality_rate <- 0
mortality_time_lag <- 0
management <- FALSE
treatment_dates <- c('2019-11-24')
treatments_file <- ""
treatment_method <- "ratio"
natural_kernel_type <- "exponential"
anthropogenic_kernel_type <- "cauchy"
natural_dir <- "NONE"
natural_kappa <- 0
anthropogenic_dir <- "NONE"
anthropogenic_kappa <- 0
pesticide_duration <- c(0)
pesticide_efficacy <- 1.0
mask <- NULL
success_metric <-
  "number of locations and total distance"
output_frequency <- "week"
output_frequency_n <- 1
movements_file <- ffIn("")
use_movements <- FALSE
start_exposed <- FALSE
generate_stochasticity <- TRUE
establishment_stochasticity <- TRUE
movement_stochasticity <- TRUE
deterministic <- FALSE
establishment_probability <- 0.5
dispersal_percentage <- 0.99
quarantine_areas_file <- ""
use_quarantine <- FALSE
random_seed <- c(42, 49, 55, 88, 95, 102, 200, 252, 294, 355)
use_spreadrates <- TRUE
calibration_method <- "ABC"
number_of_iterations <- 100
number_of_cores <- 5


config <- c()
config$random_seed <- random_seed
config$infected_file <- infected_file
config$host_file <- host_file
config$total_populations_file <- total_populations_file
config$parameter_means <- parameter_means
config$parameter_cov_matrix <- parameter_cov_matrix
config$temp <- temp
config$temperature_coefficient_file <- temperature_coefficient_file
config$precip <- precip
config$precipitation_coefficient_file <- precipitation_coefficient_file
config$model_type <- model_type
config$latency_period <- latency_period
config$time_step <- time_step
config$season_month_start <- season_month_start
config$season_month_end <- season_month_end
config$start_date <- start_date
config$end_date <- end_date
config$use_lethal_temperature <- use_lethal_temperature
config$temperature_file <- temperature_file
config$lethal_temperature <- lethal_temperature
config$lethal_temperature_month <- lethal_temperature_month
config$mortality_on <- mortality_on
config$mortality_rate <- mortality_rate
config$mortality_time_lag <- mortality_time_lag
config$management <- management
config$treatment_dates <- treatment_dates
config$treatments_file <- treatments_file
config$treatment_method <- treatment_method
config$natural_kernel_type <- natural_kernel_type
config$anthropogenic_kernel_type <- anthropogenic_kernel_type
config$natural_dir <- natural_dir
config$anthropogenic_dir <- anthropogenic_dir
config$pesticide_duration <- pesticide_duration
config$pesticide_efficacy <- pesticide_efficacy
config$output_frequency <- output_frequency
config$output_frequency_n <- output_frequency_n
config$movements_file <- movements_file
config$use_movements <- use_movements
config$start_exposed <- start_exposed
config$generate_stochasticity <- generate_stochasticity
config$establishment_stochasticity <- establishment_stochasticity
config$movement_stochasticity <- movement_stochasticity
config$deterministic <- deterministic
config$establishment_probability <- establishment_probability
config$dispersal_percentage <- dispersal_percentage
config$quarantine_areas_file <- quarantine_areas_file
config$use_quarantine <- use_quarantine
config$use_spreadrates <- use_spreadrates
config$number_of_iterations <- number_of_iterations
config$number_of_cores <- number_of_cores
# add function name for use in configuration function to skip
# function specific specifc configurations namely for validation and
# calibration.
config$function_name <- "pops"
config$failure <- NULL

config <- configuration(config)

# profvis({
# })
if (file.exists(ffIn("profile2.csv"))) {
  timing2 <- read.csv(ffIn("profile2.csv"))
} else {
  timing2 <- data.frame(host_percentage = 0, model_optimization = NA, run_time = 0)
}

# make sure to change host percentage and optimization when running
host_percentage <- 100
optimization <- "none"

for (i in seq_len(length(random_seed))) {
  
  start_time <- Sys.time()
  data <- pops_model(random_seed = config$random_seed[i],
                     use_lethal_temperature = config$use_lethal_temperature,
                     lethal_temperature = config$lethal_temperature,
                     lethal_temperature_month = config$lethal_temperature_month,
                     infected = config$infected,
                     exposed = config$exposed,
                     susceptible = config$susceptible,
                     total_populations = config$total_populations,
                     mortality_on = config$mortality_on,
                     mortality_tracker = config$mortality_tracker,
                     mortality = config$mortality,
                     quarantine_areas = config$quarantine_areas,
                     treatment_maps = config$treatment_maps,
                     treatment_dates = config$treatment_dates,
                     pesticide_duration = config$pesticide_duration,
                     resistant = config$resistant,
                     use_movements = config$use_movements,
                     movements = config$movements,
                     movements_dates = config$movements_dates,
                     weather = config$weather,
                     temperature = config$temperature,
                     weather_coefficient = config$weather_coefficient,
                     ew_res = config$ew_res,
                     ns_res = config$ns_res,
                     num_rows = config$num_rows,
                     num_cols = config$num_cols,
                     time_step = config$time_step,
                     reproductive_rate = config$reproductive_rate[1],
                     spatial_indices = config$spatial_indices,
                     mortality_rate = config$mortality_rate,
                     mortality_time_lag = config$mortality_time_lag,
                     season_month_start = config$season_month_start,
                     season_month_end = config$season_month_end,
                     start_date = config$start_date,
                     end_date = config$end_date,
                     treatment_method = config$treatment_method,
                     natural_kernel_type = config$natural_kernel_type,
                     anthropogenic_kernel_type =
                       config$anthropogenic_kernel_type,
                     use_anthropogenic_kernel =
                       config$use_anthropogenic_kernel,
                     percent_natural_dispersal =
                       config$percent_natural_dispersal[1],
                     natural_distance_scale = config$natural_distance_scale[1],
                     anthropogenic_distance_scale =
                       config$anthropogenic_distance_scale[1],
                     natural_dir = config$natural_dir,
                     natural_kappa = config$natural_kappa[1],
                     anthropogenic_dir = config$anthropogenic_dir,
                     anthropogenic_kappa = config$anthropogenic_kappa[1],
                     output_frequency = config$output_frequency,
                     output_frequency_n = config$output_frequency_n,
                     quarantine_frequency = config$quarantine_frequency,
                     quarantine_frequency_n = config$quarantine_frequency_n,
                     use_quarantine = config$use_quarantine,
                     spreadrate_frequency = config$spreadrate_frequency,
                     spreadrate_frequency_n = config$spreadrate_frequency_n,
                     use_spreadrates = config$use_spreadrates,
                     model_type_ = config$model_type,
                     latency_period = config$latency_period,
                     generate_stochasticity =
                       config$generate_stochasticity,
                     establishment_stochasticity =
                       config$establishment_stochasticity,
                     movement_stochasticity = config$movement_stochasticity,
                     deterministic = config$deterministic,
                     establishment_probability =
                       config$establishment_probability,
                     dispersal_percentage = config$dispersal_percentage
  )
  end_time <- Sys.time()
  total_time <- end_time - start_time
  
  timing2[nrow(timing2) + 1,] <- c(host_percentage, optimization, as.numeric(total_time))
  
}

write.csv(y, ffIn("profile2.csv"), row.names = FALSE)


comparison <- read.csv(ffIn("profile_combined.csv"))

library(ggplot2)
ggplot(data = comparison) +
  geom_smooth(aes(x = host_percentage, y = run_time, color = model_optimization))
