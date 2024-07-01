
####    Section 1 - Load data, AOI, and remove duplicates    ####

setwd("H:\\Meu Drive\\Projetos\\dsm_caatinga")

# load packages
library(sf) #spatial points to manipulate in r
library(tidyverse) #manipulate and tidy data
library(yardstick) #model evaluation
library(ggpubr) #extra for plotting
library(caret)
library(randomForest)
library(doParallel) 
library(gstat)
library(rgee) #GEE for R

reticulate::py_run_string("import ee; ee.Initialize()")
rgee::ee_Initialize()


# ee_Initialize()
# ee_Authenticate(auth_mode = 'gcloud')
# ee_Initialize(user = 'ramonssouza93@gmail.com')

cl <- makeCluster(detectCores()-2, type='PSOCK')
registerDoParallel(cl)

# Define the studyarea
studyarea <- ee$FeatureCollection("projects/ee-ramonssouza93/assets/bioma_250")
studyarea <- studyarea$filter(ee$Filter$eq('Bioma', 'Caatinga'))$geometry()$bounds()
print(studyarea$getInfo())


# Name Variable
name_var <- 'P' # 'C', 'N', 'P', retencao_u

UD <- 0

Data = ee$FeatureCollection(paste0("projects/gee-ufpb/assets/soils_database_",name_var))

Data = Data$filter(ee$Filter$bounds(studyarea))

Data = Data$filter(ee$Filter$neq(name_var, 'NA'))

Data = Data$filter(ee$Filter$neq(name_var, NULL))

#Data =  Data$filter(ee$Filter$gt(name_var, 0))

Data = Data$filter(ee$Filter$eq('UD', UD))

# ## Z-score method
# ## Define o limite inferior e superior
# mean =  ee$Number(Data$reduceColumns(ee$Reducer$mean(), list(name_var))$get('mean'))
# std =  ee$Number(Data$reduceColumns(ee$Reducer$stdDev(), list(name_var))$get('stdDev'))
# upper = mean$add(3)$multiply(std)
# lower = mean$subtract(3)$multiply(std)
# # Remove os outliers
# Data =  Data$filter(ee$Filter$And(ee$Filter$lt(name_var,upper), ee$Filter$gt(name_var, lower)))

computeLog <- function(feature) {
  n <- ee$Number(feature$get(name_var))$add(1)$log()
  new_var_name <- paste0(name_var, '_log') # Concatenação para criar o novo nome da variável
  feature$set(setNames(list(n), new_var_name)) # Correção aplicada aqui
}




Data <- ee$FeatureCollection(Data)$map(computeLog)
name_var <- paste0(name_var, '_log')

print(Data$first()$getInfo())


# Define uma função para adicionar as longitudes e latitudes como propriedades
add_coords <- function(feature) {
  coords <- feature$geometry()$coordinates()
  return(feature$set(list(longitude = coords$get(0), latitude = coords$get(1))))
}


# Aplica a função em cada feature da coleção
Data = Data$map(add_coords)


print(Data$first()$getInfo())



mapbiomas = ee$Image('projects/mapbiomas-workspace/public/collection8/mapbiomas_collection80_integration_v1')$clip(studyarea)

fabdem = ee$ImageCollection("projects/sat-io/open-datasets/FABDEM") # 30m
fac = ee$Image("projects/gee-ufpb/assets/flow_acu_fabdem")$clip(studyarea)$rename('FAC')
sca = ee$Image("projects/gee-ufpb/assets/SCA_fabdem")$clip(studyarea)$rename('SCA')

wc = ee$Image('projects/gee-ufpb/assets/Bioclimatic')
SBIO_0_5cm = ee$Image("projects/crowtherlab/soil_bioclim/SBIO_v2_0_5cm")
SBIO_5_15cm = ee$Image("projects/crowtherlab/soil_bioclim/SBIO_v2_5_15cm")
aridity_index_yearly = ee$Image("projects/sat-io/open-datasets/global_ai/global_ai_yearly")$rename('IA')$clip(studyarea)
et_yearly = ee$Image("projects/sat-io/open-datasets/global_et0/global_et0_yearly")$rename('ET0')$clip(studyarea)
et_yearly_sd = ee$Image("projects/sat-io/open-datasets/global_et0/global_et0_yearly_sd")$rename('ET0_stdDev')

productivity = ee$Image("projects/gee-ufpb/assets/modis_productivity")$clip(studyarea)
Ts = ee$Image("projects/gee-ufpb/assets/GSHTD")$clip(studyarea)

index_database = ee$Image("projects/dataset-landsat/assets/landsat_dataset")$clip(studyarea)

index_database = index_database$select(c(
  "B","G","R","NIR","SWIR1","SWIR2","NDVI","NDWI","GNDVI","CFLUX","EVI","MSAVI2","TCW"))

productivity = productivity$select(c('Gpp_mean',
                                     'Npp_mean',
                                     'Gpp_std',
                                     'Npp_std'))$rename(c('Gpp_mean',
                                                          'Npp_mean',
                                                          'Gpp_stdDev',
                                                          'Npp_stdDev'))

SBIO_vars <- c('SBIO1_Annual_Mean_Temperature', 'SBIO2_Mean_Diurnal_Range', 'SBIO3_Isothermality', 'SBIO4_Temperature_Seasonality', 'SBIO5_Max_Temperature_of_Warmest_Month', 'SBIO6_Min_Temperature_of_Coldest_Month', 'SBIO7_Temperature_Annual_Range', 'SBIO8_Mean_Temperature_of_Wettest_Quarter', 'SBIO9_Mean_Temperature_of_Driest_Quarter', 'SBIO10_Mean_Temperature_of_Warmest_Quarter', 'SBIO11_Mean_Temperature_of_Coldest_Quarter')
SBIO_0_5cm <- SBIO_0_5cm$select(SBIO_vars)$rename(c('SBIO1_0_5cm', 'SBIO2_0_5cm', 'SBIO3_0_5cm', 'SBIO4_0_5cm', 'SBIO5_0_5cm', 'SBIO6_0_5cm', 'SBIO7_0_5cm', 'SBIO8_0_5cm', 'SBIO9_0_5cm', 'SBIO10_0_5cm', 'SBIO11_0_5cm'))
SBIO_5_15cm <- SBIO_5_15cm$select(SBIO_vars)$rename(c('SBIO1_5_15cm', 'SBIO2_5_15cm', 'SBIO3_5_15cm', 'SBIO4_5_15cm', 'SBIO5_5_15cm', 'SBIO6_5_15cm', 'SBIO7_5_15cm', 'SBIO8_5_15cm', 'SBIO9_5_15cm', 'SBIO10_5_15cm', 'SBIO11_5_15cm'))


class_mapbiomas <- c(22 , 23 , 24 , 30 , 25 , 26 , 33 , 31 , 27 , 1 , 3 , 4 , 5 , 49 , 10 , 11 , 12 , 32 ,
                     29 , 50 , 13 , 14 , 15 , 18 , 19 , 39 , 20 , 40 , 62 , 41 , 36 , 46 , 47 , 48 , 9 , 21)

new_class <- c(0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  0 , 1 , 1 , 1 , 1 ,  1 ,  1 ,  1 ,  1 ,  1 ,  1 ,
               1 ,  1 ,  2 ,  2 ,  2 ,  2 ,  2 ,  2 ,  2 ,  2 ,  2 ,  2 ,  2 ,  2 ,  2 , 2 , 2)

mask = mapbiomas$select('classification_2022')$remap(class_mapbiomas, new_class)$rename('lulc')

# #####
# 
# | 1 | formação Florestal | 3 |
# | 2 | formação Savanica - Brasil | 4 |
# | 3 | formação Campestre | 12 |
# | 4 | pastagem | 15 |
# | 5 | agricultura | 39, 20, 40, 41, 46, 47, 48 |
# | 6 | silvicultura - Brasil | 9 |
# | 7 | mosaicoAgriculturaPastagem | 21 |
#   
# #####

lulc1985 = mapbiomas$select('classification_1985')$remap(
  c(3, 4, 12, 15, 39, 20, 40, 41, 46, 47, 48, 9, 21),
  c(1, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 6, 7))$rename('lulc_1985')
lulc2022 = mapbiomas$select('classification_2022')$remap(
  c(3, 4, 12, 15, 39, 20, 40, 41, 46, 47, 48, 9, 21),
  c(1, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 6, 7))$rename('lulc_2022')

# Apply remap function to each band
lulc = lulc1985$addBands(lulc2022)

print(lulc$bandNames()$getInfo())

lulc_1985_dummy = lulc$select(c('lulc_1985'))$eq(c(1,2,3,4,5,6,7))$rename(c('formacaoFlorestal_1985','formacaoSavanica_1985', 'formacaoCampestre_1985', 'pastagem_1985', 'agricultura_1985', 'silvicultura_1985', 'mosaicoAgriculturaPastagem_1985'))

lulc_2022_dummy = lulc$select(c('lulc_2022'))$eq(c(1,2,3,4,5,6,7))$rename(c('formacaoFlorestal_2022','formacaoSavanica_2022', 'formacaoCampestre_2022', 'pastagem_2022', 'agricultura_2022', 'silvicultura_2022', 'mosaicoAgriculturaPastagem_2022'))

print(lulc_1985_dummy$bandNames()$getInfo())

agricultura = ee$Image("projects/gee-ufpb/assets/agricultura")$select(c('agricultura_2022'))$rename(c('agricultura_2022_idade'))
formacaoCampestre = ee$Image("projects/gee-ufpb/assets/formacaoCampestre")$select(c('formacaoCampestre_2022'))$rename(c('formacaoCampestre_2022_idade'))
formacaoFlorestal = ee$Image("projects/gee-ufpb/assets/formacaoFlorestal")$select(c('formacaoFlorestal_2022'))$rename(c('formacaoFlorestal_2022_idade'))
formacaoSavanica = ee$Image("projects/gee-ufpb/assets/formacaoSavanica")$select(c('formacaoSavanica_2022'))$rename(c('formacaoSavanica_2022_idade'))
mosaicoAgriculturaPastagem = ee$Image("projects/gee-ufpb/assets/mosaicoAgriculturaPastagem")$select(c('mosaicoAgriculturaPastagem_2022'))$rename(c('mosaicoAgriculturaPastagem_2022_idade'))
pastagem = ee$Image("projects/gee-ufpb/assets/pastagem")$select(c('pastagem_2022'))$rename(c('pastagem_2022_idade'))
silvicultura = ee$Image("projects/gee-ufpb/assets/silvicultura")$select(c('silvicultura_2022'))$rename(c('silvicultura_2022_idade'))

lulc_idade = agricultura$addBands(formacaoCampestre)$addBands(formacaoFlorestal)$addBands(formacaoSavanica)$addBands(mosaicoAgriculturaPastagem)$addBands(pastagem)$addBands(silvicultura)


terrain_metrics = ee$Image("projects/ee-rss/assets/terrainMetrics")$select(c(
  'Elevation', 'Slope', 'Northness', 'Eastness', 'HorizontalCurvature', 'VerticalCurvature', 'ShapeIndex'))$rename(c(
    'Elev', 'Declividade', 'Northness', 'Eastness', 'CurvHor', 'CurvVert', 'ShapeIndex'))

proj = terrain_metrics$select(0)$projection()

scale_proj = proj$nominalScale()$getInfo()

slope_deg = terrain_metrics$select('Declividade')

slope_p = slope_deg$divide(180)$multiply(ee$Image(pi))$tan()$multiply(100)

slope_rad = slope_deg$multiply(ee$Image(pi)$divide(180))

slope_rad = slope_rad$add(slope_rad$lte(0)$multiply(ee$Image(0.001)))

spi = fac$multiply(slope_rad$tan())$rename("SPI")$clip(studyarea)

twi = ((fac$divide(slope_rad$tan()))$log())$rename("TWI")$clip(studyarea)

elev = terrain_metrics$select("Elev")$clip(studyarea)

tpi = elev$subtract(elev$focalMean(50, 'square', 'pixels'))$rename("TPI")

topographic_database = terrain_metrics$select(c('Elev', 'Northness', 'Eastness', 'CurvHor', 'CurvVert', 'ShapeIndex'))$addBands(slope_p)$addBands(twi)$addBands(spi)$addBands(tpi)

wc = wc$rename(c('BIO1', 'BIO2', 'BIO3', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9', 'BIO10','BIO11', 'BIO12', 'BIO13', 'BIO14', 'BIO15', 'BIO16', 'BIO17', 'BIO18', 'BIO19'))


climate_database = wc$addBands(SBIO_0_5cm)$addBands(SBIO_5_15cm)$addBands(aridity_index_yearly)$addBands(et_yearly)$addBands(et_yearly_sd)$addBands(productivity)$addBands(Ts)$clip(studyarea)

climate_database = climate_database$resample('bicubic')$reproject(index_database$projection())$clip(studyarea)

# print(climate_database$bandNames()$getInfo())

geology = ee$FeatureCollection ("projects/gee-ufpb/assets/geol_area")
geology = geology$filter(ee$Filter$bounds(studyarea))
geology = geology$filter(ee$Filter$neq('nm_provinc', "Corpo d'água continental"))
propertyValues = geology$aggregate_array('nm_provinc')$distinct()

# Converte a lista de valores em uma lista de strings
class_names <- propertyValues$getInfo()

# Função para criar uma imagem para cada classe
create_class_image <- function(class_name) {
  # Filtra a coleção para a classe específica
  class_feature_collection <- geology$filter(ee$Filter$eq('nm_provinc', class_name))
  
  # Cria uma imagem binária para a classe (pixels da classe têm valor 1, outros têm valor 0)
  class_image <- ee$Image(0)$byte()$paint(class_feature_collection, 1)
  
  # Define o nome da banda como o nome da classe
  class_image <- class_image$rename(class_name)
  
  return(class_image)
}

# Mapeia a função sobre a lista de nomes de classe para obter uma lista de imagens
class_images <- lapply(class_names, create_class_image)

# Converte a lista de imagens em uma coleção de imagens
geology_dummy <- ee$ImageCollection(class_images)$toBands()

geology_dummy = geology_dummy$select(c(
  "0_Borborema","1_Cobertura Cenozoica","2_Costeira e Margem Continental","3_Mantiqueira","4_Recôncavo-Tucano-Jatobá","5_São Francisco","6_Parnaíba", "7_São Luís"))$rename(c(
    "Borborema","Cobertura_Cenozoica","Costeira_e_Margem_Continental","Mantiqueira","Reconcavo_Tucano_Jatoba","Sao_Francisco","Parnaiba", "Sao_Luis"))

# Imprime a coleção de imagens
print(geology_dummy$bandNames()$getInfo())

# Map$centerObject(studyarea, 4)
# 
# Map$addLayer(geology_dummy$clip(studyarea), list(
#   bands = c('Borborema'),
#   min = 0,
#   max = 1,
#   palette=c('black', 'green')), 'Borborema')


soil = ee$FeatureCollection("projects/gee-ufpb/assets/pedo_area")
soil = soil$filter(ee$Filter$bounds(studyarea))
soil = soil$filter(ee$Filter$neq('leg_ordem', "CORPO D'ÁGUA CONTINENTAL"))
soil = soil$filter(ee$Filter$neq('leg_ordem', 'OUTROS'))
soil = soil$filter(ee$Filter$neq('leg_ordem', 'DUNAS'))
propertyValues = soil$aggregate_array('leg_ordem')$distinct()

# Converta a lista de valores em uma lista de strings
class_names <- propertyValues$getInfo()


# Função para criar uma imagem para cada classe
create_class_image <- function(class_name) {
  # Filtra a coleção para a classe específica
  class_feature_collection <- soil$filter(ee$Filter$eq('leg_ordem', class_name))
  
  # Cria uma imagem binária para a classe (pixels da classe têm valor 1, outros têm valor 0)
  class_image <- ee$Image(0)$byte()$paint(class_feature_collection, 1)
  
  # Define o nome da banda como o nome da classe
  class_image <- class_image$rename(class_name)
  
  return(class_image)
}

# Mapeia a função sobre a lista de nomes de classe para obter uma lista de imagens
class_images <- lapply(class_names, create_class_image)

# Converte a lista de imagens em uma coleção de imagens
soil_dummy <- ee$ImageCollection(class_images)$toBands()

soil_dummy = soil_dummy$select(c(
  "0_PLINTOSSOLO", "1_NITOSSOLO", "2_ORGANOSSOLO", "3_VERTISSOLO", "4_ARGISSOLO", "5_CAMBISSOLO",
  "6_CHERNOSSOLO", "7_ESPODOSSOLO", "8_GLEISSOLO", "9_LATOSSOLO", "10_LUVISSOLO", "11_NEOSSOLO",
  "12_PLANOSSOLO"))$rename(c(
    "Plintossolo", "Nitossolo", "Organossolo", "Vertissolo", "Argissolo", "Cambissolo",
    "Chernossolo", "Espodossolo", "Gleissolo", "Latossolo", "LuvisSolo", "Neossolo",
    "Planossolo"))

print(soil_dummy$bandNames()$getInfo())

geom = ee$FeatureCollection ("projects/gee-ufpb/assets/geom_area")
geom = geom$filter(ee$Filter$bounds(studyarea))
geom = geom$filter(ee$Filter$neq('compartime', "Corpo d'água continental"))
propertyValues = geom$aggregate_array('compartime')$distinct()
propertyValues$getInfo()

# Converta a lista de valores em uma lista de strings
class_names = propertyValues$getInfo()


# Função para criar uma imagem para cada classe
create_class_image <- function(class_name) {
  # Filtra a coleção para a classe específica
  class_feature_collection <- geom$filter(ee$Filter$eq('compartime', class_name))
  
  # Cria uma imagem binária para a classe (pixels da classe têm valor 1, outros têm valor 0)
  class_image <- ee$Image(0)$byte()$paint(class_feature_collection, 1)
  
  # Define o nome da banda como o nome da classe
  class_image <- class_image$rename(class_name)
  
  return(class_image)
}

# Mapeia a função sobre a lista de nomes de classe para obter uma lista de imagens
class_images <- lapply(class_names, create_class_image)


# Converte a lista de imagens em uma coleção de imagens
geom_dummy <- ee$ImageCollection(class_images)$toBands()

geom_dummy = geom_dummy$select(c(
  "0_Depressões", "1_Chapadas", "2_Patamares", "3_Planaltos", "4_Tabuleiros", "5_Serras", "6_Planícies"))$rename(c(
    "Depressoes", "Chapadas", "Patamares", "Planaltos", "Tabuleiros", "Serras", "Planicies"))

print(geom_dummy$bandNames()$getInfo())


# Cria um raster com as coordenadas x e y
xy_image = ee$Image$pixelLonLat()$clip(studyarea)

# Set of variables to be used in the model

stack = topographic_database$addBands(climate_database)$addBands(index_database)$addBands(xy_image)
#$addBands(lulc_1985_dummy)$addBands(lulc_2022_dummy)$addBands(geology_dummy)$addBands(geom_dummy)$addBands(soil_dummy)
# stack = stack$updateMask(mask)
print(stack$bandNames()$getInfo())

# Add a random column (named random) and specify seed value for repeatability.
Data = Data$randomColumn()

# Função para remover a chave 'random'
remove_random <- function(feature) {
  return(feature$set('random', NULL))
}

# Separe 75% para treinamento, 25% para validação
split <- 0.75
training <- Data$filter(ee$Filter$lt('random', split))
training = training$map(remove_random)

testing <- Data$filter(ee$Filter$gte('random', split))
testing = testing$map(remove_random)
samples_testing_df <- ee_as_sf(testing)

print(training$first()$getInfo())

samples_training <- stack$sampleRegions(collection = training, properties = list(name_var), scale = scale_proj, tileScale = 16, geometries = TRUE)
samples_training_df <- ee_as_sf(samples_training)

# Remova a variável específica ("name_var") da lista usando drop

samples_training_df <- samples_training_df %>% select(-name_var)

samples_training_df <- samples_training_df %>% st_drop_geometry()

cor_matrix <- cor(samples_training_df, method = 'spearman')

low_corr <- cor_matrix[,-findCorrelation(cor_matrix,cutoff=0.95)] # remove high cor variable

length(colnames(low_corr))

colnames(low_corr)

s <- topographic_database$addBands(climate_database)$addBands(index_database)$addBands(xy_image)$addBands(lulc_1985_dummy)$addBands(lulc_2022_dummy)$addBands(geology_dummy)$addBands(geom_dummy)$addBands(soil_dummy)$addBands(lulc_idade)
# names_rasters <- s$bandNames()$getInfo()
# write.table(names_rasters, './names_rasters.csv', row.names = FALSE, sep = ";", dec = ',')

samples_training <- s$sampleRegions(collection = training, properties = list(name_var), scale = scale_proj, tileScale = 16, geometries = TRUE)

samples_training_df <- ee_as_sf(samples_training)

samples_training_df <- samples_training_df %>% st_drop_geometry()

write.table(samples_training_df, paste0('./results/tables/Training/', name_var,'_samples_training.csv'), row.names = FALSE, sep = ";", dec = ',')

stack <- stack$select(colnames(low_corr))$addBands(lulc_1985_dummy)$addBands(lulc_2022_dummy)$addBands(geology_dummy)$addBands(geom_dummy)$addBands(soil_dummy)$addBands(lulc_idade)

# stack <- stack$updateMask(mask)

print(stack$bandNames()$getInfo())

contagem <- length(stack$bandNames()$getInfo())

print(contagem)

bands <- stack$bandNames()


samples_training <- stack$sampleRegions(collection = training, properties = list(name_var), scale = scale_proj, tileScale = 16, geometries = TRUE)

samples_training_df <- ee_as_sf(samples_training)

samples_training_df <- samples_training_df %>% st_drop_geometry()

# write.table(samples_training_df, paste0('./results/tables/Training/', name_var,'_samples_training.csv'), row.names = FALSE, sep = ";", dec = ',')

samples_training_df <- samples_training_df %>% select(-name_var)



samples_testing <- stack$sampleRegions(collection = testing, properties = list(name_var), scale = scale_proj, tileScale = 16, geometries = TRUE)

samples_testing_df <- ee_as_sf(samples_testing)

samples_testing_df <- samples_testing_df %>% select(-name_var)

samples_testing_df <- samples_testing_df %>% st_drop_geometry()

# write.table(samples_testing_df, paste0('./results/tables/Testing/', name_var,'_samples_testing.csv'), row.names = FALSE, sep = ";", dec = ',')

k <- 10

seedList <- seq(1, 50, 1)
cat("The seeds are:", seedList, "\n")

# for(seed in seedList){
# 
#   subampleWithCovariatesRaw <- stack$reduceRegions(collection = training, reducer = ee$Reducer$first(), tileScale = 16)
# 
#   subampleWithCovariates <- subampleWithCovariatesRaw$filter(ee$Filter$notNull(subampleWithCovariatesRaw$first()$propertyNames()))
# 
#   subampleWithCovariatesAndFold <- subampleWithCovariates$randomColumn('CV_Fold', seed)$map(function(f) {f$set('CV_Fold', ee$Number(f$get('CV_Fold'))$multiply(k)$toInt())})
# 
#   # cat(subampleWithCovariates$size()$getInfo(), "\n")
# 
#   trainTableWithCovarites_Export <- ee$batch$Export$table$toAsset(
#     collection = subampleWithCovariatesAndFold,
#     description = paste0('Train_Table_seed_', seed, '_Exportation', sep = ''),
#     assetId = paste0('projects/ee-ramonssouza93/assets/Data/Gridsampled_Train_Table_', name_var, '_seed_', seed, sep = '')
#   )
# 
#   trainTableWithCovarites_Export$start()
# }
# cat('Covariates extraction is done!\n')


# generate the classifier list based on fullParameterSpace
classifierListsGenerator <- function(parameterSets, randomDiscrete = TRUE, 
                                     randomNumber = 12, nTrees = 20, 
                                     modelType = 'REGRESSION', bagFraction = 0.632, 
                                     Seed = 0) {
  
  # define an empty list to load the defined models for grid search
  classifierList <- list()
  
  if (randomDiscrete) {
    # check the randomNumber
    if (is.null(randomNumber)) {
      cat('Warning! an integer number needs to be allocated to <randomNumber>!\n')
    } else {
      cat(paste('A randomDiscrete approach has been applied to do grid search the parameter space!\n',
                'The random model number is: ', randomNumber, '!\n'))
    }
    
    # subset the fullParameterSpace randomly with the randomNumber
    set.seed(Seed)
    randomParameterApplied <- sample(parameterSets, randomNumber, replace = TRUE)
    
  } else {
    cat('The full space of the parameter sets is being running for grid search\n')
    set.seed(Seed)
    randomParameterApplied <- sample(parameterSets, randomNumber, replace = TRUE)
  }
  
  cat(Seed, '\n')
  cat('function use 20 as the default nTrees,\n',
      'You can define your own nTree value in the function argument settings!\n')
  
  # loop through the randomParameterApplied
  for (ParaSet in randomParameterApplied) {
    model_name <- paste('GridSeach_Model_', ParaSet[1], '_', ParaSet[2], '_', ParaSet[3], sep = '')
    
    # define the parameter setting of each model in the grid search 
    # and allocate those parameters into the feature
    # define the paramter setting of each model in the grid seach and allocate those parameters into the feature
    perRF <- ee$Feature(ee$Geometry$Point(c(0,0)))$set('ModelName',model_name,'PerClassifier',ee$Classifier$smileRandomForest(
      numberOfTrees = nTrees,
      variablesPerSplit = ParaSet[1],
      minLeafPopulation = ParaSet[2],
      maxNodes = ParaSet[3],
      bagFraction = bagFraction)$setOutputMode(modelType))
    
    
    classifierList <- append(classifierList, perRF)
  }
  
  return(classifierList)
}


# Define a função R^2 para uso com modelos de valores contínuos (ou seja, modelos baseados em regressão)
coefficientOfDetermination <- ee_utils_pyfunc(function(anyVariableTable, propertyOfInterest, propertyOfInterest_Predicted) {
  
  # Calcula a média da propriedade de interesse
  propertyOfInterestMean <- ee$Number(ee$Dictionary(ee$FeatureCollection(anyVariableTable)$select(c(propertyOfInterest))$reduceColumns(ee$Reducer$mean(), list(propertyOfInterest)))$get('mean'))
  
  # Calcula a soma total dos quadrados
  totalSoSFunction <- ee_utils_pyfunc(function(f) {
    return(f$set('Difference_Squared', ee$Number(ee$Feature(f)$get(propertyOfInterest))$subtract(propertyOfInterestMean)$pow(ee$Number(2L))))
  })
  
  totalSumOfSquares <- ee$Number(ee$Dictionary(ee$FeatureCollection(anyVariableTable)$map(totalSoSFunction)$select(c('Difference_Squared'))$reduceColumns(ee$Reducer$sum(), list('Difference_Squared')))$get('sum'))
  
  # Calcula a soma residual dos quadrados
  residualSoSFunction <- ee_utils_pyfunc(function(f) {
    return(f$set('Residual_Squared', ee$Number(ee$Feature(f)$get(propertyOfInterest))$subtract(ee$Number(ee$Feature(f)$get(propertyOfInterest_Predicted)))$pow(ee$Number(2L))))
  })
  
  residualSumOfSquares <- ee$Number(ee$Dictionary(ee$FeatureCollection(anyVariableTable)$map(residualSoSFunction)$select(c('Residual_Squared'))$reduceColumns(ee$Reducer$sum(), list('Residual_Squared')))$get('sum'))
  
  # Finaliza o cálculo
  r2 <- ee$Number(1L)$subtract(residualSumOfSquares$divide(totalSumOfSquares))
  
  return(ee$Number(r2))
})

# RMSE function
RMSE <- ee_utils_pyfunc(function(anyVariableTable,propertyOfInterest,propertyOfInterest_Predicted){
  # Compute the squared difference between observed and predicted
  propDiff <- ee_utils_pyfunc(function(f){
    diff = ee$Number(f$get(propertyOfInterest))$subtract(ee$Number(f$get(propertyOfInterest_Predicted)))
    
    return(f$set('diff', diff$pow(2L)))
  })
  # calculate RMSE from squared difference
  rmse = ee$Number(anyVariableTable$map(propDiff)$reduceColumns(ee$Reducer$mean(), list('diff'))$get('mean'))$sqrt()
  
  return(rmse)
})

# MAE function
MAE <- ee_utils_pyfunc(function(anyVariableTable,propertyOfInterest,propertyOfInterest_Predicted){
  # Compute the absolute difference between observed and predicted
  propDiff <- ee_utils_pyfunc(function(f){
    diff = ee$Number(f$get(propertyOfInterest))$subtract(ee$Number(f$get(propertyOfInterest_Predicted)))
    
    return(f$set('diff', diff$abs()))
  })
  # calculate RMSE from squared difference
  mae = ee$Number(anyVariableTable$map(propDiff)$reduceColumns(ee$Reducer$mean(), list('diff'))$get('mean'))
  
  return(mae)
})

# Define the function to calculate Lin's Concordance Correlation Coefficient (LCCC)
LCCC <- ee_utils_pyfunc(function(anyVariableTable, propertyOfInterest, propertyOfInterest_Predicted){
  # Select the properties of interest
  selected_properties = anyVariableTable$select(c(propertyOfInterest, propertyOfInterest_Predicted))
  
  # Calculate the Pearson correlation coefficient
  correlation = selected_properties$reduceColumns(ee$Reducer$pearsonsCorrelation(), list(propertyOfInterest, propertyOfInterest_Predicted))$get('correlation')
  
  # Compute the variance of observed and predicted values
  varObserved = ee$Number(anyVariableTable$reduceColumns(ee$Reducer$variance(), list(propertyOfInterest))$get('variance'))
  varPredicted = ee$Number(anyVariableTable$reduceColumns(ee$Reducer$variance(), list(propertyOfInterest_Predicted))$get('variance'))
  
  # Compute the mean of observed and predicted values
  meanObserved = ee$Number(anyVariableTable$reduceColumns(ee$Reducer$mean(), list(propertyOfInterest))$get('mean'))
  meanPredicted = ee$Number(anyVariableTable$reduceColumns(ee$Reducer$mean(), list(propertyOfInterest_Predicted))$get('mean'))
  
  # Compute Lin's Concordance Correlation Coefficient (LCCC)
  lccc = ee$Number(2L)$multiply(correlation)$multiply(varObserved$sqrt())$multiply(varPredicted$sqrt())$divide(varObserved$add(varPredicted)$add((meanObserved$subtract(meanPredicted))$pow(2L)))
  
  return(lccc)
})

computeCVAccuracy <- function(featureWithClassifier,
                              propertyOfInterest,
                              modelType,
                              kFoldAssignmentFC,
                              cvFoldString,
                              classProperty,
                              extractedVariableTable) {
  
  # Extrair o classificador do recurso
  cOI <- ee$Classifier(featureWithClassifier$get('PerClassifier'))
  
  # Criar uma função para mapear através das atribuições de dobras e calcular a precisão geral
  # para todas as dobras de validação
  computeAccuracyForFold <- function(foldFeature) {
    # Organizar os dados de treinamento e validação
    foldNumber <- ee$Number(ee$Feature(foldFeature)$get('Fold'))
    trainingData <- extractedVariableTable$filterMetadata(cvFoldString, 'not_equals', foldNumber)
    validationData <- extractedVariableTable$filterMetadata(cvFoldString, 'equals', foldNumber)
    
    # Treinar o classificador e classificar o conjunto de dados de validação
    trainedClassifier <- cOI$train(trainingData, classProperty, propertyOfInterest)
    outputtedPropName <- paste0(classProperty, '_Predicted')
    classifiedValidationData <- validationData$classify(trainedClassifier, outputtedPropName)
    
    # Criar uma instrução central if/then que determina o tipo de valores de precisão que são retornados
    if (modelType == 'CLASSIFICATION') {
      # Calcular a precisão geral da classificação
      errorMatrix <- classifiedValidationData$errorMatrix(classProperty, outputtedPropName)
      overallAccuracy <- ee$Number(errorMatrix$accuracy())
      return(foldFeature$set(list(accuracyMetricString = overallAccuracy)))
    } else {
      # Calcular o R^2 da regressão
      r2ToSet <- coefficientOfDetermination(classifiedValidationData, classProperty, outputtedPropName)
      rmseToSet <- RMSE(classifiedValidationData, classProperty, outputtedPropName)
      maeToSet <- MAE(classifiedValidationData, classProperty, outputtedPropName)
      lccToSet <- LCCC(classifiedValidationData, classProperty, outputtedPropName)
      return(foldFeature$set(list(R2 = r2ToSet, RMSE = rmseToSet, MAE = maeToSet, LCCC = lccToSet)))
    }
  }
  
  # Calcular os valores de precisão do classificador em todas as dobras
  accuracyFC <- kFoldAssignmentFC$map(computeAccuracyForFold)
  meanAccuracy <- accuracyFC$aggregate_mean('R2')
  sdAccuracy <- accuracyFC$aggregate_total_sd('R2')
  
  # Calcular média e desvio padrão de RMSE
  meanRMSE <- accuracyFC$aggregate_mean('RMSE')
  sdRMSE <- accuracyFC$aggregate_total_sd('RMSE')
  
  # Calcular média e desvio padrão de MAE
  meanMAE <- accuracyFC$aggregate_mean('MAE')
  sdMAE <- accuracyFC$aggregate_total_sd('MAE')
  
  # Calcular média e desvio padrão de LCCC
  meanLCCC <- accuracyFC$aggregate_mean('LCCC')
  sdLCCC <- accuracyFC$aggregate_total_sd('LCCC')
  
  # Computar o recurso para retornar
  featureToReturn <- featureWithClassifier$select(list('cName'))$set(list(
    Mean_R2 = meanAccuracy,
    StDev_R2 = sdAccuracy,
    Mean_RMSE = meanRMSE,
    StDev_RMSE = sdRMSE,
    Mean_MAE = meanMAE,
    StDev_MAE = sdMAE,
    Mean_LCCC = meanLCCC,
    StDev_LCCC = sdLCCC
  ))
  
  return(featureToReturn)
}

propertyOfInterest <- stack$bandNames()$getInfo()

gridSearchEarthEngine <- function(inputTrainTable, # train data table in ee$FeatureCollection format
                                  propertyOfInterest, # list of predictors
                                  classProperty = 'C', # response variable name in Google earth engine
                                  nTrees = 20, # number of trees, default is 100
                                  variablesPerSplitList = seq(3, 21, by = 3), # list
                                  minLeafPopulationList = seq(2, 20, by = 2), # list
                                  maxNodesList = seq(10, 100, by = 10), # list
                                  bagFraction = 0.632,
                                  randomDiscrete = TRUE, # boolean
                                  randomNumber = 1, # if random discrete is True, you must set this value
                                  foldsValue = 10,
                                  modelType = 'REGRESSION',
                                  cvFoldString = 'CV_Fold',
                                  pyramidingPolicy = 'mean',
                                  Seeds = 0) {
  
  parameterLists <- expand.grid(variablesPerSplitList, minLeafPopulationList, maxNodesList)
  # generate the list of all the possible parameter set combinations
  fullParameterSpace <- as.list(parameterLists)
  # generate the classifier in featureCollection format
  classifierList <- classifierListsGenerator(parameterSets = fullParameterSpace,
                                             randomNumber = randomNumber,
                                             nTrees = nTrees,
                                             bagFraction = bagFraction,
                                             Seed = Seeds)
  
  kList <- seq(0, foldsValue - 1)
  kFoldAssignmentFC <- ee$FeatureCollection(lapply(kList, function(n) {
    ee$Feature(ee$Geometry$Point(c(0, 0)), list(Fold = n))
  }))
  
  classDf <- data.frame(Mean_MAE = numeric(),
                        StDev_MAE = numeric(),
                        ModelName = character(),
                        numberOfTrees = numeric(),
                        variablesPerSplit = numeric(),
                        minLeafPopulation = numeric(),
                        bagFraction = numeric(),
                        maxNodes = numeric())
  
  for (rf in classifierList) {
    accuracy_feature <- ee$Feature(computeCVAccuracy(rf,
                                                     propertyOfInterest,
                                                     modelType = 'REGRESSION',
                                                     kFoldAssignmentFC = kFoldAssignmentFC,
                                                     cvFoldString = cvFoldString,
                                                     classProperty = classProperty,
                                                     extractedVariableTable = inputTrainTable))
    # extract the parameter information
    parameterDict <- rf$getInfo()$properties$PerClassifier$classifier
    parameterDF <- as.data.frame(parameterDict)
    # extract the metrics information
    metricDict <- accuracy_feature$getInfo()$properties
    metricDF <- as.data.frame(metricDict)
    
    resultDF <- cbind(metricDF, parameterDF)
    classDf <- rbind(classDf, resultDF)
  }
  
  # sort the grid search result by ascending of Mean_MAE
  classDfSorted <- classDf[order(classDf$Mean_MAE), ]
  
  return(head(classDfSorted, 1))
}

######

# Number of repetitions for RFE execution (> 1 for repeated cross-validation)
rep_rfe <- 1

# Number of folds in RFE (for cross-validation)
fold_rfe <- 10

# Number os subsets of covariates in the RFE 
size_rfe <- c(seq(10, 49, 2), seq(50, 120, 5))
# size_rfe<- c(1:ncol(samplesT)-1)

# Used models (see available models in caret documentation)
models <- c("rf")

# Metric used for model optimization
metric_otm <- "MAE"

# Number of tries for hyperparameter tuning in the RFE
tn_length_rfe <- 4

vars = stack$bandNames()$getInfo()

# Gerar uma lista ee para salvar as sementes
seedList <- seq(7, 50, 1)

# Definir a lista de variáveis dependentes
vpdList <- list(name_var)

########

for (seed in seedList) {
  for (vpd in vpdList) {
    # Para ambos os escaladores, mean e max, os mesmos dados de treinamento estão sendo usados para modelagem
    inputVariableTable <- ee$FeatureCollection(paste0('projects/ee-ramonssouza93/assets/Data/Gridsampled_Train_Table_', name_var, '_seed_', seed))
    
    samplesT <- ee_as_sf(inputVariableTable)
    
    samplesT <- samplesT %>% st_drop_geometry()
    
    samplesT <- samplesT %>% select(all_of(vars), all_of(name_var))
    
    # Escrever no diretório local
    write.table(samplesT, paste0('./results/tables/', name_var,'_seed_' ,seed, '_samples_training.csv'), row.names = FALSE, sep = ";", dec = ',')
    
    X <- samplesT %>% select(vars)
    
    y <- samplesT %>% select(name_var)
    
    set.seed(seed)
    
    # Recursive Feature Elimination
    # Define the control using a random forest selection function
    rfe_ctrl <- rfeControl(functions = rfFuncs, # random forest
                           method = "repeatedcv", # repeated cv
                           repeats = rep_rfe, # number of repeats
                           number = fold_rfe) # number of folds
    
  
    
    formu <- as.formula(paste(name_var, "~ ."))
    
    
    set.seed(seed)
    
    rfe_fit <- rfe(form = formu,
                   data = samplesT,
                   sizes = size_rfe,
                   method = models,
                   metric = metric_otm,
                   # trControl = model_ctrl,
                   # tuneLength = tn_length_rfe,
                   rfeControl = rfe_ctrl,
                   maximize = ifelse(metric_otm %in% c("RMSE", "MAE"),
                                     FALSE, TRUE))
    
    
    # Escrever no diretório local
    write.table(rfe_fit$results, paste0('results/tables/RFE/', vpd, '_rfe_fit_results_Seed_', seed, '.csv'), row.names = FALSE, sep = ";", dec = ',')
    
    features_otm <- samplesT %>% select(rfe_fit$optVariables)
    
    # Escrever no diretório local
    write.table(features_otm, paste0('results/tables/propertyOfInterest/', vpd, '_propertyOfInterest_Seed_', seed, '.csv'), row.names = FALSE, sep = ";", dec = ',')
    
    propertyOfInterest <- colnames(features_otm)
    
    cat(propertyOfInterest, "\n")
    
    topModelParameter <- gridSearchEarthEngine(inputTrainTable = inputVariableTable,
                                               propertyOfInterest = propertyOfInterest,
                                               classProperty = vpd,
                                               randomNumber = 48L,
                                               nTrees = 250L,
                                               Seeds = seed)
    
    # Escrever a tabela dos melhores parâmetros no diretório local
    write.table(topModelParameter, paste0('results/GridSearchResults/', vpd, '_global_Modeling_Grid_subsample_Grid_Search_Seed_', seed, '.csv'), row.names = FALSE, sep = ";", dec = ',')
    
    # Mostrar o progresso da busca em grade pelo número da semente
    cat(paste0('Busca em grade para a semente: ', seed, ' está concluída!'), "\n")
  }
}

# # Gerar uma lista ee para salvar as sementes
seedList <- seq(1, 50, 1)
# 
# # Definir a lista de variáveis dependentes
# vpdList <- list('C_log')

for (seed in seedList) {
  for (vpd in vpdList) {
    # Carregar os dados dos pontos com os covariáveis
    trainTable <- ee$FeatureCollection(paste0('projects/ee-ramonssouza93/assets/Data/Gridsampled_Train_', name_var))
    
    # Ler a tabela de parâmetros
    parameterTable <- read.table(paste0('results/GridSearchResults/', vpd, '_global_Modeling_Grid_subsample_Grid_Search_best.csv'), sep = ';', dec = ',', header=TRUE)
    
    
    # Extrair os parâmetros
    variablesPerSplitVal <- as.integer(parameterTable$variablesPerSplit[1]) # mtry
    minLeafPopulationVal <- as.integer(parameterTable$minLeafPopulation[1]) # minrow
    maxNodesVal <- as.integer(parameterTable$maxNodes[1]) # max depth
    cat('seed', seed, variablesPerSplitVal, minLeafPopulationVal, maxNodesVal, '\n')
    
    # Definir o classificador random forest
    rfClassifier <- ee$Classifier$smileRandomForest(numberOfTrees = 250,
                                                    variablesPerSplit = variablesPerSplitVal, # mtry
                                                    minLeafPopulation = minLeafPopulationVal, # minrow
                                                    maxNodes = maxNodesVal, # max depth
                                                    bagFraction = 0.632,
                                                    seed = seed)$setOutputMode('REGRESSION')
    
    # Ler a propriedade de interesse
    
    propertyOfInterest <- colnames(read.table(paste0('results/tables/propertyOfInterest/', vpd, '_propertyOfInterest_Seed_', seed, '.csv'), sep = ';', dec = ',', header=TRUE))
    
    # Treinar o classificador
    trainedClassifier <- rfClassifier$train(features = trainTable,
                                            classProperty = vpd,
                                            inputProperties = propertyOfInterest)
    
    # Obter a importância das variáveis
    dict <- trainedClassifier$explain()
    variable_importance <- ee$Feature(NULL, ee$Dictionary(dict$get('importance')))

    # Exportar a importância das variáveis para CSV
    variable_importance <- ee_as_sf(variable_importance)

    # Escrever no diretório local
    # write.table(variable_importance, paste0('./results/tables/varImp/RF_', vpd, '_', seed, '.csv'), row.names = FALSE, sep = ";", dec = ',')

    # Executar a previsão para gerar o mapa
    ImgMap <- stack$select(propertyOfInterest)$classify(trainedClassifier)
    
    # samples_t <- ImgMap$sampleRegions(
    #   collection = training,
    #   scale = scale_proj,
    #   tileScale = 16,
    #   geometries = TRUE
    # )
    # 
    # computeRFresidual <- function(feature) {
    #   o <- ee$Number(feature$get(name_var))
    #   p <- ee$Number(feature$get('classification'))
    #   
    #   r <- o$subtract(p)
    #   return(feature$set(list(residual = r)))
    # }
    # 
    # samples_t <- samples_t$map(computeRFresidual)
    # 
    # print(samples_t$first()$getInfo())
    # 
    # Map$addLayer(samples_t, list(NULL), name_var)
    # 
    # # Combine os redutores de média e desvio padrão para eficiência.
    # combinedReducer <- ee$Reducer$mean()$combine(
    #   reducer2 = ee$Reducer$stdDev(),
    #   sharedInputs = TRUE
    # )
    # 
    # # Estime a média global e o desvio padrão (SD) dos pontos.
    # stats <- samples_t$reduceColumns(
    #   reducer = combinedReducer,
    #   selectors = list('residual')
    # )
    # print(stats$getInfo())
    # 
    # soilImg <- samples_t$inverseDistance(
    #   range = 50 * 10000,
    #   propertyName = 'residual',
    #   mean = stats$get('mean'),
    #   stdDev = stats$get('stdDev'),
    #   gamma = 0.5
    # )$reproject(crs = "EPSG:4326", scale = scale_proj)
    # 
    # print(soilImg$bandNames()$getInfo())
    # 
    # ImgMap <- ImgMap$add(soilImg)
    # 
    # # Map$centerObject(studyarea, 4)
    # # print(soilImg$bandNames()$getInfo())
    # Map$addLayer(ImgMap$clip(studyarea), list(
    #   bands = list("classification"),
    #   min = 0.6,
    #   max = 3,
    #   palette=c('#ffffbf','#1a9641','#d7191c')), name_var)

    
    samples_test <- ImgMap$sampleRegions(
      collection = testing,
      scale = scale_proj,
      tileScale = 16
    )
    
    
    # Calculate the Pearson correlation coefficient
    correlation = samples_test$reduceColumns(ee$Reducer$pearsonsCorrelation(), selectors = c(vpd, 'classification'))$get('correlation')
    
    
    samples_test <- samples_test$reduceColumns(
      reducer = ee$Reducer$toList(2),
      selectors = c(vpd, 'classification')
    )$get('list')
    
    
    obs1 = ee$Array(samples_test)$transpose()$cut(list(0, -1))$project(list(1))
    pred1 = ee$Array(samples_test)$transpose()$cut(list(1, -1))$project(list(1))
    rmse = obs1$subtract(pred1)$pow(2)$reduce('mean', list(0))$sqrt()$get(list(0))
    mae = obs1$subtract(pred1)$abs()$reduce('mean', list(0))$get(list(0))
    
    # Calcular R-squared (R²)
    meanObs = obs1$reduce('mean', list(0))$get(list(0))
    totalSS = obs1$subtract(meanObs)$pow(2)$reduce('sum', list(0))$get(list(0))
    residualSS = obs1$subtract(pred1)$pow(2)$reduce('sum', list(0))$get(list(0))
    
    r2 = ee$Number(1)$subtract(residualSS$divide(totalSS))
    
    # Calcular a média das observações
    meanObserved <- obs1$reduce('mean', list(0))$get(list(0))
    
    # Calcular a média das previsões
    meanPredicted <- pred1$reduce('mean', list(0))$get(list(0))
    
    # Calcular a covariância entre as observações e as previsões
    covObsPred <- obs1$subtract(meanObserved)$multiply(pred1$subtract(meanPredicted))$reduce('sum', list(0))$get(list(0))
    
    # Calcular a variância das observações
    varObserved <- obs1$subtract(meanObserved)$pow(2)$reduce('sum', list(0))$get(list(0))
    
    # Calcular a variância das previsões
    varPredicted <- pred1$subtract(meanPredicted)$pow(2)$reduce('sum', list(0))$get(list(0))
    
    # Calcular o Lin's Concordance Correlation Coefficient (LCCC)
    lccc <- ee$Number(2)$multiply(correlation)$multiply(varObserved$sqrt())$multiply(varPredicted$sqrt())$divide(varObserved$add(varPredicted)$add((meanObserved$subtract(meanPredicted))$pow(2)))

    Metrics_testing <- ee$FeatureCollection(ee$Feature(NULL,list(mae = mae, rmse = rmse, r2 = r2, lccc = lccc)))
    
    Metrics_testing_df <- ee_as_sf(Metrics_testing)
    print(Metrics_testing_df)
    # write.table(Metrics_testing_df, paste0('results/tables/Testing/', vpd, '_metrics_testing_Seed_', seed, '.csv'), row.names = FALSE, sep = ";", dec = ',')
    
    # # Configurar a exportação da imagem
    # predictionExport <- ee$batch$Export$image$toAsset(image = ImgMap$exp()$subtract(1)$toFloat(),
    #                                                   description = paste0(vpd, '_Map_To_Asset_', seed),
    #                                                   assetId = paste0('projects/data-soil-ch/assets/Imgs/EnsambledMaps/Predicted_', vpd, '_Seed_', seed),
    #                                                   scale = scale_proj,
    #                                                   region = studyarea,
    #                                                   crs = 'EPSG:4326',
    #                                                   #crsTransform = c(0.008333333333333333, 0, -180, 0, -0.008333333333333333, 90),
    #                                                   maxPixels = 1e13)
    # 
    # # Iniciar a tarefa de exportação
    # predictionExport$start()

  }
}



# # Reverte log1p para o valor normal
# revertLog1p <- function(log1p_value) {
#   original_value <- expm1(log1p_value)
#   return(original_value)
# }
# 
# revertLog1p(2.7561333)


for (vpd in vpdList) {
  # Carrega as imagens previstas pelos modelos de conjunto
  # Define uma imagem vazia
  firstImage <- ee$Image(paste0('projects/data-soil-ch/assets/Imgs/EnsambledMaps/Predicted_', vpd, '_Seed_1'))$rename('Model_1')$toFloat()
  # Carrega as outras imagens e adiciona-as como bandas à primeira imagem acima
  modelList <- seq(2, 50, 1)
  for (ml in modelList) {
    perModelImage <- ee$Image(paste0('projects/data-soil-ch/assets/Imgs/EnsambledMaps/Predicted_', vpd, '_Seed_', ml))$rename(paste0('Model_', ml))$toFloat()
    firstImage <- firstImage$addBands(perModelImage)
  }
  
  # modelList <- seq(25, 50, 1)
  # for (ml in modelList) {
  #   perModelImage <- ee$Image(paste0('projects/data-soil-p/assets/Imgs/EnsambledMaps/Predicted_', vpd, '_Seed_', ml))$rename(paste0('Model_', ml))$toFloat()
  #   firstImage <- firstImage$addBands(perModelImage)
  # }
  # 
  # img1 <- ee$Image('projects/images-c/assets/Imgs/EnsambledMaps/Predicted_P_log_Seed_21')
  # img2 <- ee$Image('projects/images-c/assets/Imgs/EnsambledMaps/Predicted_P_log_Seed_22')
  # img3 <- ee$Image('projects/images-c/assets/Imgs/EnsambledMaps/Predicted_P_log_Seed_23')
  # img4 <- ee$Image('projects/images-c/assets/Imgs/EnsambledMaps/Predicted_P_log_Seed_24')
  # 
  # firstImage <- firstImage$addBands(img1)$addBands(img2)$addBands(img3)$addBands(img4)
  # 
  cat('The band names are:', firstImage$bandNames()$getInfo(), '\n')
  
  
  # Calcula as imagens de média e variação
  meanImage <- firstImage$reduce(ee$Reducer$mean())
  variImage <- firstImage$reduce(ee$Reducer$stdDev())$divide(meanImage)
  # Obtem o percentil 95%
  percentileImage <- firstImage$reduce(ee$Reducer$percentile(c(5, 95), c('lower', 'upper')))
  
  # Adiciona essas duas imagens aos ativos do GEE
  meanExport <- ee$batch$Export$image$toAsset(image = meanImage$toFloat(),
                                              description =  paste0('Predicted_', vpd, '_Ensambled_Mean'),
                                              assetId = paste0('projects/soil-data-sdm/assets/Predicted_', vpd, '_Ensambled_Mean'),
                                              scale = scale_proj,
                                              region = studyarea,
                                              crs = 'EPSG:4326',
                                              # crsTransform = c(0.008333333333333333, 0, -180, 0, -0.008333333333333333, 90),
                                              maxPixels = 1e13)
  
  # Inicia a tarefa de exportação
  meanExport$start()

  variExport <- ee$batch$Export$image$toAsset(image = variImage$toFloat(),
                                              description = paste0('Predicted_', vpd, '_Ensambled_Variation_Coefficient'),
                                              assetId = paste0('projects/soil-data-sdm/assets/Predicted_', vpd, '_Ensambled_Variation_Coefficient'),
                                              scale = scale_proj,
                                              region = studyarea,
                                              crs = 'EPSG:4326',
                                              # crsTransform = c(0.008333333333333333, 0, -180, 0, -0.008333333333333333, 90),
                                              maxPixels = 1e13)
  
  # Inicia a tarefa de exportação
  variExport$start()

  percentileExport <- ee$batch$Export$image$toAsset(image = percentileImage$toFloat(),
                                                    description = paste0('Predicted_', vpd, '_Ensambled_Percentile'),
                                                    assetId = paste0('projects/soil-data-sdm/assets/Predicted_', vpd, '_Ensambled_Percentile'),
                                                    scale = scale_proj,
                                                    region = studyarea,
                                                    crs = 'EPSG:4326',
                                                    # crsTransform = c(0.008333333333333333, 0, -180, 0, -0.008333333333333333, 90),
                                                    maxPixels = 1e13)
  
  # Inicia a tarefa de exportação
  percentileExport$start()
}

# Imprime a informação de que a exportação está sendo executada no Google Earth Engine
cat('Export is running on Google Earth Engine!\nPlease check it on the Google Earth Engine UI.\n')
