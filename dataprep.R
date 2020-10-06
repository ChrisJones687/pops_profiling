library(PoPS)
library(folderfun)
library(raster)
library(sf)
library(fasterize)

cols <- 500
rows <- 500
cells <- cols * rows

size10 <- round(cells *.1)

host100 <- matrix(round(runif(cells, 10, 70), digits = 0), nrow = rows, ncol = cols)
samp90 <- 
  sample(seq(1,cells,1), size = size10 * 9, replace = FALSE, prob = NULL)
samp80 <- 
  sample(samp90, size = size10 * 8, replace = FALSE, prob = NULL)
samp70 <- 
  sample(samp80, size = size10 * 7, replace = FALSE, prob = NULL)
samp60 <- 
  sample(samp70, size = size10 * 6, replace = FALSE, prob = NULL)
samp50 <- 
  sample(samp60, size = size10 * 5, replace = FALSE, prob = NULL)
samp40 <- 
  sample(samp50, size = size10 * 4, replace = FALSE, prob = NULL)
samp30 <- 
  sample(samp40, size = size10 * 3, replace = FALSE, prob = NULL)
samp20 <- 
  sample(samp30, size = size10 * 2, replace = FALSE, prob = NULL)
samp10 <- 
  sample(samp20, size = size10 * 1, replace = FALSE, prob = NULL)

host10 <- host100
host10[samp90] <- 0

host20 <- host100
host20[samp80] <- 0

host30 <- host100
host30[samp70] <- 0

host40 <- host100
host40[samp60] <- 0

host50 <- host100
host50[samp50] <- 0

host60 <- host100
host60[samp40] <- 0

host70 <- host100
host70[samp30] <- 0

host80 <- host100
host80[samp20] <- 0

host90 <- host100
host90[samp10] <- 0

sampinf <- 
  sample(samp10, size = 30, replace = FALSE, prob = NULL)

int_infected <- host10
int_infected[] <- 0
int_infected[sampinf] <- 2

resolution <- 100
host100_rast <- raster(host100, xmn = 0, xmx = nrow(host100) * resolution, 
                      ymn = 0, ymx = ncol(host90) * resolution)
host90_rast <- raster(host90, xmn = 0, xmx = nrow(host90) * resolution, 
                      ymn = 0, ymx = ncol(host90) * resolution)
host80_rast <- raster(host80, xmn = 0, xmx = nrow(host80) * resolution, 
                      ymn = 0, ymx = ncol(host90) * resolution)
host70_rast <- raster(host70, xmn = 0, xmx = nrow(host70) * resolution, 
                      ymn = 0, ymx = ncol(host70) * resolution)
host60_rast <- raster(host60, xmn = 0, xmx = nrow(host60) * resolution, 
                      ymn = 0, ymx = ncol(host90) * resolution)
host50_rast <- raster(host50, xmn = 0, xmx = nrow(host50) * resolution, 
                      ymn = 0, ymx = ncol(host50) * resolution)
host40_rast <- raster(host40, xmn = 0, xmx = nrow(host40) * resolution, 
                      ymn = 0, ymx = ncol(host90) * resolution)
host30_rast <- raster(host30, xmn = 0, xmx = nrow(host30) * resolution, 
                      ymn = 0, ymx = ncol(host30) * resolution)
host20_rast <- raster(host20, xmn = 0, xmx = nrow(host20) * resolution, 
                      ymn = 0, ymx = ncol(host90) * resolution)
host10_rast <- raster(host10, xmn = 0, xmx = nrow(host10) * resolution, 
                      ymn = 0, ymx = ncol(host10) * resolution)

total_population <- host10_rast
total_population[] <- 100

int_infected_rast <- raster(int_infected, xmn = 0, xmx = nrow(int_infected) * resolution, 
                      ymn = 0, ymx = ncol(int_infected) * resolution)

write.csv(int_infected, "data/int_infected.csv", col.names = FALSE, row.names = FALSE)
write.csv(host100, "data/host100.csv")
write.csv(host90, "data/host90.csv")
write.csv(host80, "data/host80.csv")
write.csv(host70, "data/host70.csv")
write.csv(host60, "data/host60.csv")
write.csv(host50, "data/host50.csv")
write.csv(host40, "data/host40.csv")
write.csv(host30, "data/host30.csv")
write.csv(host20, "data/host20.csv")
write.csv(host10, "data/host10.csv")
writeRaster(total_population, "data/total_population.tif")
writeRaster(host100_rast, "data/host100.tif", overwrite = TRUE)
writeRaster(host90_rast, "data/host90.tif", overwrite = TRUE)
writeRaster(host80_rast, "data/host80.tif", overwrite = TRUE)
writeRaster(host70_rast, "data/host70.tif", overwrite = TRUE)
writeRaster(host60_rast, "data/host60.tif", overwrite = TRUE)
writeRaster(host50_rast, "data/host50.tif", overwrite = TRUE)
writeRaster(host40_rast, "data/host40.tif", overwrite = TRUE)
writeRaster(host30_rast, "data/host30.tif", overwrite = TRUE)
writeRaster(host20_rast, "data/host20.tif", overwrite = TRUE)
writeRaster(host10_rast, "data/host10.tif", overwrite = TRUE)
writeRaster(int_infected_rast, "data/int_infected.tif", overwrite = TRUE)

s <- as.matrix(read.csv("data/int_infected.csv"))
