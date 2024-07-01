library(reshape2) # Para a função melt
library(ggplot2)
library(scales)
library(dplyr)
library(forcats) # Para fct_reorder
library(tidyr)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(forcats)
library(GD)
library(colorspace)

setwd('F:\\dsm_caatinga\\results\\tables\\varImp')

# Defina suas variáveis
vars <- c('C', 'N', 'CN', 'P')

all_df <- data.frame()

names <- read.table("F:\\dsm_caatinga\\names_rasters.csv", sep = ';', dec = ',', header = T)

for (name_var in vars) {
  current_files <- list.files(pattern = paste0('RF_', name_var, '_.*\\.csv$'), full.names = TRUE)
  for (f in current_files) {
    print(f)
    df <- read.table(f, sep = ';', dec = ',', header = T)
    df <- df %>% select(-c('geometry'))
    colunas_maior_soma <- names(sort(colSums(df), decreasing = TRUE)[1:5])
    colunas_maior_soma <- colunas_maior_soma[!is.na(colunas_maior_soma)]
    df <- df %>% select(colunas_maior_soma)
    df['file'] = data.frame(f)
    df['Variable'] <- name_var
    all_df <- bind_rows(all_df, df)
  }
}

all_df <- all_df %>% select(-c('file'))


# Usando dplyr para selecionar a coluna Variable
# e melt para transformar as outras colunas em linhas
reshaped_df <- all_df %>%
  select(Variable, everything()) %>%
  melt(id.vars = "Variable")

# Primeiro, certifique-se de que as colunas 'variable' em ambos os dataframes são do mesmo tipo
reshaped_df$variable <- as.character(reshaped_df$variable)
names$Variables <- as.character(names$Variables)

# Agora, você pode usar a função 'match' para encontrar os índices correspondentes
# e atualizar 'variable' em 'df1' com 'new_names' de 'df2'
reshaped_df$variable <- names$New_Variables[match(reshaped_df$variable, names$Variables)]

reshaped_df <- reshaped_df[complete.cases(reshaped_df$value), ]

names_u <- unique(reshaped_df$variable)

# Criando um vetor de índices com base nos nomes
indices <- match(names$New_Variables, names_u)

indices <- indices[!is.na(indices)]

# Ordenando df1 com base nos índices
names_u <- names_u[indices]

# Criando uma sequência de IDs
ids <- paste0("X", seq_along(names_u))

# Reverte log1p para o valor normal
revertLog1p <- function(log1p_value) {
  original_value <- expm1(log1p_value)
  return(original_value)
}

# Definir as cores para cada variável
colors <- c("C" = "#e65c00", "N" = "#38A800", "CN" = "#2B74F9", "P" = "#C1286C")

# Gerar gráficos para cada variável
plots <- list()

for (var in vars) {

  # Ler a tabela de parâmetros
  data_soil <- read.table(paste0("F:/dsm_caatinga/results/tables/Training/",var,"_log_samples_training.csv"), sep = ';', dec = ',', header = TRUE)

  colunas <- c("B","G","R","NIR","SWIR1","SWIR2","NDVI","NDWI","GNDVI","CFLUX","EVI","MSAVI2","TCW")
  
  # Converter as colunas para numérico
  data_soil[colunas] <- lapply(data_soil[colunas], function(x) as.numeric(as.character(x)))
  
  colunas_selecionadas <- names(data_soil %>% select(names$Variables[88:129]))
  
  # Converter as colunas para numérico
  data_soil[colunas_selecionadas] <- lapply(data_soil[colunas_selecionadas], function(x) as.character(as.character(x)))
  
  df <- all_df[grepl(paste0("^",var, "$"), all_df$Variable), ]
  
  df <- df[, colSums(is.na(df)) != nrow(df)]
  
  df <- df %>% select(-c('Variable'))
  
  names_c <- colnames(df)
  
  data <- data_soil %>% select(c(all_of(names_c), paste0(var, "_log")))
  
  # set optional parameters of optimal discretization
  discmethod <- c("equal","natural","quantile","geometric","sd")
  
  discitv <- c(3:8)
  
  # continuous_variable <- colnames(data)[-c(70:111,119)]
  classes <- sapply(data, class)
  classes <- classes[!names(classes) %in% paste0(var, "_log")]
  continuous_variable <- names(classes[classes %in% c("numeric")])
  print(continuous_variable)
  
  cgdm <- GD::gdm(as.formula(paste(paste0(var, "_log"), "~ .")),
                  continuous_variable = continuous_variable,
                  data = data,
                  discmethod = discmethod, discitv = discitv)
  
  # Primeiro, certifique-se de que as colunas 'variable' em ambos os dataframes são do mesmo tipo
  cgdm$Factor.detector$Factor$variable <- as.character(cgdm$Factor.detector$Factor$variable )
  
  # Agora, você pode usar a função 'match' para encontrar os índices correspondentes
  # e atualizar 'variable' em 'df1' com 'new_names' de 'df2'
  cgdm$Factor.detector$Factor$variable <- names$New_Variables[match(cgdm$Factor.detector$Factor$variable, names$Variables)]
  
  cgdm$Factor.detector$Factor$variable <- ids[match(cgdm$Factor.detector$Factor$variable, names_u)]
  fd <- cgdm$Factor.detector$Factor
  
  # fd <- cgdm$Factor.detector$Factor %>%
  #   filter(sig <=  0.05)
  
  # Adiciona uma nova coluna 'label' com asteriscos para sig <= 0.05
  fd$label <- ifelse(fd$sig <= 0.05, "*", "")
  
  # # Utiliza a nova coluna 'label' para os rótulos no ggplot
  # p <- ggplot(fd, aes(x = qv, y = reorder(variable, qv))) +
  #   geom_bar(stat = "identity", fill = colors[var]) +
  #   labs(x = "QV", y = "") +
  #   theme_minimal() +
  #   scale_fill_manual(values = colors) +
  #   scale_x_continuous(labels = label_number()) +
  #   theme(axis.title.x = element_text(size = 18, face = "bold"),
  #         axis.title.y = element_text(size = 18, face = "bold"),
  #         axis.text.x = element_text(size = 18, face = "bold"),
  #         axis.text.y = element_text(size = 18, face = "bold"),
  #         plot.margin = unit(c(1, 3, 1, 1), "lines"))
  
  p <- ggplot(fd, aes(x = qv, y = reorder(variable, qv))) +
    geom_bar(stat = "identity", fill = colors[var]) +
    geom_text(aes(label = label), vjust = 0.7, hjust = -0.1, color = "black", size = 10) + # Adiciona esta linha
    labs(x = "QV", y = "") +
    theme_minimal() +
    scale_fill_manual(values = colors) +
    scale_x_continuous(labels = label_number()) +
    theme(axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"),
          plot.margin = unit(c(1, 3, 1, 1), "lines"))
  
  plots[[var]] <- p
  
  inter_dect <- cgdm$Interaction.detector$Interaction %>% select(c("var1", "var2","qv12", "interaction"))
  
  melted_cormat <- melt(inter_dect, id.vars = c("var1", "var2","qv12", "interaction"), na.rm = TRUE)
  
  
  
  # Write the data frame to a file
  write.table(melted_cormat %>% arrange(desc(qv12)), file=paste0("F:\\dsm_caatinga\\results\\tables\\geodetector\\qv12_", var, ".csv"), 
              na = '',  row.names = FALSE, sep = ";" , dec = ",")
  
  # Primeiro, certifique-se de que as colunas 'variable' em ambos os dataframes são do mesmo tipo
  melted_cormat$var1 <- as.character(melted_cormat$var1)
  melted_cormat$var2 <- as.character(melted_cormat$var2)
  
  melted_cormat$var1 <- names$New_Variables[match(melted_cormat$var1, names$Variables)]
  melted_cormat$var2 <- names$New_Variables[match(melted_cormat$var2, names$Variables)]
  
  melted_cormat$var1 <- ids[match(melted_cormat$var1, names_u)]
  melted_cormat$var2 <- ids[match(melted_cormat$var2, names_u)]
  
  melted_cormat$var1
  
  var1 = unique(melted_cormat$var1)
  var2 = unique(melted_cormat$var2)
  
  mat = matrix(0,nrow = length(var1), ncol = length(var2))
  rownames(mat) = var1
  colnames(mat) = var2
  
  
  for (i in 1:nrow(melted_cormat))
  {
    mat[melted_cormat$var1[i], melted_cormat$var2[i]] = melted_cormat$qv12[i]
  }
  
  mat2 = matrix(0,nrow = length(var1), ncol = length(var2))
  rownames(mat2) = var1
  colnames(mat2) = var2
  
  for(i in 1:nrow(melted_cormat))
  {
    mat2[melted_cormat$var1[i], melted_cormat$var2[i]] = melted_cormat$interaction[i]
  }
  
  mat[mat <= 0] <- NA
  
  mat2[mat2 == "0"] <- ""
  
  # Invertendo a ordem das linhas
  mat <- mat[nrow(mat):1, ]
  
  # Invertendo a ordem das colunas
  mat <- mat[, ncol(mat):1]
  
  # Invertendo a ordem das linhas
  mat2 <- mat2[nrow(mat2):1, ]
  
  # Invertendo a ordem das colunas
  mat2 <- mat2[, ncol(mat2):1]
  
  w <- as.integer(length(colnames(mat))*0.8)
  h <- as.integer(length(colnames(mat))*0.6)
  
  # cell_fun = function(j, i, x, y, w, h, fill) {
  #   if (!is.na(mat2[i, j]) && mat2[i, j] == "Enhance, nonlinear") {
  #     grid.text(paste0(sprintf("%.2f", mat[i, j]), "\n***"), x, y, gp = gpar(fontsize = 14, fontface = "bold"))
  #   } else if (!is.na(mat2[i, j]) && mat2[i, j] == "Enhance, bi-") {
  #     grid.text(paste0(sprintf("%.2f", mat[i, j]), "\n**"), x, y, gp = gpar(fontsize = 14, fontface = "bold"))
  #   } else if (!is.na(mat2[i, j]) && mat2[i, j] == "Weaken, uni-") {
  #     grid.text(paste0(sprintf("%.2f", mat[i, j]), "\n*"), x, y, gp = gpar(fontsize = 14, fontface = "bold"))
  #   }
  # }
  
  # cell_fun = function(j, i, x, y, w, h, fill) {
  #   if (!is.na(mat2[i, j]) && mat2[i, j] == "Enhance, nonlinear") {
  #     grid.text(paste0(sprintf("%.2f", mat[i, j]), "\n*"), x, y, gp = gpar(fontsize = as.integer(length(colnames(mat))*0.8), fontface = "bold"))
  #   } else if (!is.na(mat2[i, j]) && mat2[i, j] == "Enhance, bi-") {
  #     grid.text(paste0(sprintf("%.2f", mat[i, j]), "\n**"), x, y, gp = gpar(fontsize = as.integer(length(colnames(mat))*0.8), fontface = "bold"))
  #   }
  # }
  
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(mat2[i, j]) && mat2[i, j] == "Enhance, nonlinear") {
      grid.text(paste0(sprintf("%.2f", mat[i, j]), "\n*"), x, y, gp = gpar(fontsize = 12, fontface = "bold"))
    } else if (!is.na(mat2[i, j]) && mat2[i, j] == "Enhance, bi-") {
      grid.text(paste0(sprintf("%.2f", mat[i, j]), "\n**"), x, y, gp = gpar(fontsize = 12, fontface = "bold"))
    }
  }
  
  # Ajuste das margens
  padding <- unit(c(2, 2, 2, 2), "cm") # top, right, bottom, left
  
  # hp <- ComplexHeatmap::Heatmap(mat, cluster_rows = F, cluster_columns = F,
  #                               column_dend_side = "bottom",
  #                               col = hcl.colors(5, palette = "Earth"),
  #                               heatmap_legend_param = list(
  #                                 title = "Interaction Value",
  #                                 legend_height = unit(as.integer(length(colnames(mat))*0.4), "cm"),
  #                                 legend_width = unit(as.integer(length(colnames(mat))*0.2), "cm"),
  #                                 title_position = "leftcenter-rot",
  #                                 legend_gp = gpar(col = "black", fontsize = as.integer(length(colnames(mat))*0.6), fontface = "bold"), # Ajuste o valor de fontsize conforme necessário para a legenda
  #                                 labels_gp = gpar(col = "black", fontsize = as.integer(length(colnames(mat))*0.6), fontface = "bold"),
  #                                 title_gp = gpar(col = "black", fontsize = as.integer(length(colnames(mat))*0.6), fontface = "bold")
  #                               ),
  #                               cell_fun = cell_fun, na_col = "white",
  #                               row_names_side = "left",
  #                               # Ajustando o tamanho da fonte das etiquetas das linhas e colunas
  #                               row_names_gp = gpar(fontsize = as.integer(length(colnames(mat))*0.9), fontface = "bold"), # Aumentado para 20
  #                               column_names_gp = gpar(fontsize = as.integer(length(colnames(mat))*0.9), fontface = "bold") # Aumentado para 20
  # )
  
  hp <- ComplexHeatmap::Heatmap(mat, cluster_rows = F, cluster_columns = F,
                                column_dend_side = "bottom",
                                col = hcl.colors(5, palette = "Earth"),
                                heatmap_legend_param = list(
                                  title = "Interaction Value",
                                  legend_height = unit(6, "cm"),
                                  legend_width = unit(3, "cm"),
                                  title_position = "leftcenter-rot",
                                  legend_gp = gpar(col = "black", fontsize = 10, fontface = "bold"), # Ajuste o valor de fontsize conforme necessário para a legenda
                                  labels_gp = gpar(col = "black", fontsize = 10, fontface = "bold"),
                                  title_gp = gpar(col = "black", fontsize = 10, fontface = "bold")
                                ),
                                cell_fun = cell_fun, na_col = "white",
                                row_names_side = "left",
                                # Ajustando o tamanho da fonte das etiquetas das linhas e colunas
                                row_names_gp = gpar(fontsize = 14, fontface = "bold"), # Aumentado para 20
                                column_names_gp = gpar(fontsize = 14, fontface = "bold") # Aumentado para 20
  )
  # lgd_list = list(
  #   Legend( labels = c("Enhance, nonlinear", "Enhance, bi-", "Weaken, uni-"), title = "",  labels_gp = gpar(col = "black", fontsize = 20, fontface = "bold"),
  #           graphics = list(
  #             function(x, y, w, h) grid.text("***", x = x, y = y,
  #                                            gp = gpar(fill = "black", fontsize = 20, fontface = "bold")),
  #             function(x, y, w, h) grid.text("**", x = x, y = y,
  #                                            gp = gpar(fill = "black", fontsize = 20, fontface = "bold")),
  #             function(x, y, w, h) grid.text("*", x = x, y = y,
  #                                            gp = gpar(fill = "black", fontsize = 20, fontface = "bold")))
  #   ))
  
  # lgd_list = list(
  #   Legend( labels = c("Enhance, nonlinear", "Enhance, bi-"), title = "",  labels_gp = gpar(col = "black", fontsize = as.integer(length(colnames(mat))*0.7), fontface = "bold"),
  #           graphics = list(
  #             function(x, y, w, h) grid.text("*", x = x, y = y,
  #                                            gp = gpar(fill = "black", fontsize = as.integer(length(colnames(mat))*0.7), fontface = "bold")),
  #             function(x, y, w, h) grid.text("**", x = x, y = y,
  #                                            gp = gpar(fill = "black", fontsize = as.integer(length(colnames(mat))*0.7), fontface = "bold")))
  #   ))
  
  lgd_list = list(
    Legend( labels = c("Enhance, nonlinear", "Enhance, bi-"), title = "",  labels_gp = gpar(col = "black", fontsize = 14, fontface = "bold"),
            graphics = list(
              function(x, y, w, h) grid.text("*", x = x, y = y,
                                             gp = gpar(fill = "black", fontsize = 14, fontface = "bold")),
              function(x, y, w, h) grid.text("**", x = x, y = y,
                                             gp = gpar(fill = "black", fontsize = 14, fontface = "bold")))
    ))
  # tiff(paste("H:\\Meu Drive\\Projetos\\dsm_caatinga\\results\\", var,"_Interaction2.tiff",sep=""), units="in", width= w, height= h, res=600, compression = "lzw")
  tiff(paste("F:\\dsm_caatinga\\results\\", var,"_Interaction2.tiff",sep=""), units="in", width= 13, height= 9, res=600, compression = "lzw")
  
  # Agora vamos criar o boxplot
  draw(hp, annotation_legend_list = lgd_list, ht_gap = unit(1, "cm"), merge_legend = TRUE, padding = padding)
  dev.off()
}

figure <- ggarrange(plotlist = plots,
                    labels = c("C", "N", "C/N", "P"),
                    ncol = 2, nrow = 2,
                    font.label = list(size = 26, face = "bold"),
                    label.x = 0.5, label.y = 1.02,
                    common.legend = TRUE, legend = "bottom")


tiff(paste("F:\\dsm_caatinga\\results\\Factor_detector.tiff",sep=""), units="in", width=15, height=20, res=600, compression = "lzw")
figure + theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5, "mm"))
dev.off()
