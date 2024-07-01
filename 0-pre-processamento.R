
# carregar pacotes
library(sf)
library(raster)
library(rgdal)
library(fasterize)
library(tmap)
# library(ggspatial)
library(tidyverse)
library(tidyr)
# library(GSIF)
library(aqp)
library(zoo)
library(mpspline2)
# library(ithir)
library(raster)
library(sp)
library(dplyr)

# diretorios --------------------------------------------------------------
# conferir diretorio



setwd("H:\\Meu Drive\\Projetos\\ml_soil_caatinga")
# setwd("/mnt/Windows/Users/ramon/OneDrive/Documentos/01-Principal/Projetos/Projetos_R/dsm_caatinga/Pontos_amostragem_PronaSolos_2020/")
getwd()

caatinga <- st_read("./caatinga_wgs84.shp")

data_all <- read_sf("./soils_database.shp")

data_all <- st_intersection(data_all ,caatinga)

colnames(data_all)
# data_all = subset(data_all, select = -c(prof) )

# Criando um vetor com as variáveis
vars <- c('C', 'N', 'CN', 'P', 'CP', 'NP', "densidade_", "areia_tota", "silte", "argila", "ph_h2o", "ph_kcl")

# Loop for para obter os valores das variáveis
for (var in vars) {
  print(var)
  # Selecione as colunas dinamicamente
  data <- data_all %>%
    select(c(codigo_pon, prof_s, prof_i, all_of(var), Latitude, Longitude))
  
  # data = subset(data, select = c(codigo_pon, prof_s, prof_i, var, Latitude, Longitude) )
  print(data)
  data <- data %>% 
    na.omit() %>% drop_na()
  
  data$prof <- data$prof_i - data$prof_s
  
  prof_i <- data[data$prof < 0,]$prof_s
  
  prof_s <- data[data$prof < 0,]$prof_i
  
  data[data$prof < 0,]$prof_s <- prof_s
  
  data[data$prof < 0,]$prof_i <- prof_i
  
  data = subset(data, select = -c(prof))
  
  data <- transform(data, sites = as.numeric(interaction(Latitude, Longitude, drop=TRUE)))
  
  
  
  dat <- data.frame("SID" = data$codigo_pon,
                    "UD" = data$prof_s,
                    "LD" = data$prof_i,
                    "VAL" = data[[var]],
                    "LAT" = data$Latitude,
                    "LONG" = data$Longitude,
                    stringsAsFactors = FALSE)
  
  
  # dat <- dat[dat$LD > 0,]
  
  spl_dat <- mpspline_tidy(obj = dat, d = c(0, 20, 40, 60, 80, 100, 200), var_name = 'VAL')
  
  dat_u <- unique(dat[c("SID", "LAT", "LONG" )])
  
  spl_dat_est_dcm <- spl_dat$est_dcm

  data_merge <-inner_join(dat_u, spl_dat_est_dcm, by = "SID")
  
  names(data_merge)[names(data_merge) == "SPLINED_VALUE"] <- var

  data_sf = st_as_sf(data_merge, coords = c("LONG", "LAT"), crs = 4326)
  
  st_crs(data_sf)

  st_write(data_sf,paste0("soils_database_", var, ".shp"), layer_options = "ENCODING=UTF-8", delete_layer = TRUE)
  
}

# Exibindo os valores inseridos pelo usuário
cat("\nValores inseridos:\n")
for (var in names(valores)) {
  cat(paste(var, ": ", valores[[var]], "\n"))
}


data = subset(data, select = c(codigo_pon, prof_s, prof_i, ph_kcl, Latitude, Longitude) )

data <- data %>% 
  na.omit() %>% drop_na()

data$prof <- data$prof_i - data$prof_s

prof_i <- data[data$prof < 0,]$prof_s

prof_s <- data[data$prof < 0,]$prof_i

data[data$prof < 0,]$prof_s <- prof_s

data[data$prof < 0,]$prof_i <- prof_i

data = subset(data, select = -c(prof))

data <- transform(data, sites = as.numeric(interaction(Latitude, Longitude, drop=TRUE)))



dat <- data.frame("SID" = data$codigo_pon,
                  "UD" = data$prof_s,
                  "LD" = data$prof_i,
                  "VAL" = data$ph_kcl,
                  "LAT" = data$Latitude,
                  "LONG" = data$Longitude,
                  stringsAsFactors = FALSE)


dat <- dat[dat$LD > 0,]

# spl_dat <- mpspline_tidy(obj = dat, d = c(0, 5, 15, 30, 60, 100,  200), var_name = 'VAL')

# (0–20 cm, 20–40cm, 40–60 cm, 60–80 cm and 80–100cm)

# spl_dat <- mpspline_tidy(obj = dat, d = c(0, 10, 30, 100, 200), var_name = 'VAL')

spl_dat <- mpspline_tidy(obj = dat, d = c(0, 20, 40, 60, 80, 100, 200), var_name = 'VAL')

# spl_dat <- mpspline_tidy(obj = dat, d = c(0, 30, 100, 200), var_name = 'VAL')

dat_u <- unique(dat[c("SID", "LAT", "LONG" )])

# dat_u$SID <- factor(dat_u$SID)

spl_dat_est_dcm <- spl_dat$est_dcm

# spl_dat_est_dcm$SID <- factor(spl_dat_est_dcm$SID)

data_merge <-inner_join(dat_u, spl_dat_est_dcm, by = "SID")

# 
# summary(data_merge[data_merge$LD <=5,])
# 
# count(data_merge[data_merge$LD <=5,])


names(data_merge)[names(data_merge) == "SPLINED_VALUE"] <- 'ph_kcl'

#data_merge$ORC_log <- log(data_merge$ORC)


data_sf = st_as_sf(data_merge, coords = c("LONG", "LAT"), crs = 4326)

st_crs(data_sf)

#data_sf$ORC_log[ data_sf$ORC_log == "-Inf"] <- NA
plot(data_sf)

st_write(data_sf, "soils_database_ph_kcl.shp", layer_options = "ENCODING=UTF-8", delete_layer = TRUE)

