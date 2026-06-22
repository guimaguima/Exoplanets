# install.packages(c("tidyverse", "corrplot", "car"))

library(tidyverse)
library(corrplot)
library(car)        
library(skimr)
library(broom)
library(patchwork)

sink('eda.txt')

kepler_data <- read_csv("data/cumulative.csv")

kepler_clean <- kepler_data %>%
  select(
    koi_disposition, 
    koi_period, koi_impact, koi_duration, koi_depth, 
    koi_prad, koi_teq, koi_insol, 
    koi_steff, koi_slogg, koi_srad
  ) %>%
  drop_na() %>%
  filter(koi_disposition %in% c("CONFIRMED", "FALSE POSITIVE")) %>%
  mutate(
    koi_disposition = as.factor(koi_disposition),
    koi_disposition = fct_relevel(koi_disposition, "CONFIRMED")
  )

print(skim(kepler_clean))

cat("\n--- Features Transformadas ---\n")

kepler_scaled <- kepler_clean %>%
  mutate(
    log_koi_period   = log10(koi_period + 1),
    log_koi_prad     = log10(koi_prad + 1),
    log_koi_insol    = log10(koi_insol + 1),
    log_koi_depth    = log10(koi_depth + 1),
    log_koi_srad     = log10(koi_srad + 1),
    log_koi_teq      = log10(koi_teq + 1),
    log_koi_impact   = log10(koi_impact + 1),
    log_koi_duration = log10(koi_duration + 1)
  ) %>%
  mutate(across(
    c(log_koi_period, log_koi_prad, log_koi_insol, log_koi_depth, 
      log_koi_srad, log_koi_teq, log_koi_impact, log_koi_duration, 
      koi_steff, koi_slogg), 
    ~ ( . - median(., na.rm = TRUE) ) / mad(., na.rm = TRUE), 
    
    .names = "{.col}_z"
  )) %>%
  select(koi_disposition, ends_with("_z"))


print(skim(kepler_scaled))

# Análise Multivariada de Colinearidade (VIF)
cat("\n--- Colinearidade (VIF) ---\n")

glm_dummy <- glm(koi_disposition ~ ., data = kepler_scaled, family = binomial())
vif_values <- car::vif(glm_dummy)
print(vif_values)


cat("\n--- Remover koi_teq ---\n")
kepler_scaled_small <- kepler_scaled %>%
  select(
    koi_disposition,
    log_koi_period_z, 
    log_koi_prad_z,
    log_koi_insol_z,
    log_koi_depth_z,
    log_koi_srad_z,
    # log_koi_teq_z,
    log_koi_impact_z, 
    log_koi_duration_z,
    koi_steff_z, 
    koi_slogg_z
    
    # Remover koi_teq "temperatura de equilíbrio do planeta" (653 VIF)
    # É proporcional e deve ser calculado 
    # usando koi_steff "temperatura da estrela", koi_srad "raio da estrela" e talvez outras variáveis
  )

glm_dummy <- glm(koi_disposition ~ ., data = kepler_scaled_small, family = binomial())
vif_values <- car::vif(glm_dummy)
print(vif_values)

cat("\n--- Remover koi_insol ---\n")
kepler_scaled_small <- kepler_scaled %>%
  select(
    koi_disposition,
    log_koi_period_z, 
    log_koi_prad_z,
    # log_koi_insol_z,
    log_koi_depth_z,
    log_koi_srad_z,
    # log_koi_teq_z,
    log_koi_impact_z, 
    log_koi_duration_z,
    koi_steff_z, 
    koi_slogg_z
    
    # Remover koi_insol "fluxo de insolação" (69 VIF)
    # É proporcional e deve ser calculado 
    # usando koi_steff "temperatura da estrela", koi_srad "raio da estrela" e talvez outras variáveis
  )

glm_dummy <- glm(koi_disposition ~ ., data = kepler_scaled_small, family = binomial())
vif_values <- car::vif(glm_dummy)
print(vif_values)

cat("\n--- Remover koi_depth ---\n")
kepler_scaled_small <- kepler_scaled %>%
  select(
    koi_disposition,
    log_koi_period_z, 
    log_koi_prad_z,
    # log_koi_insol_z,
    # log_koi_depth_z,
    log_koi_srad_z,
    # log_koi_teq_z,
    log_koi_impact_z, 
    log_koi_duration_z,
    koi_steff_z, 
    koi_slogg_z
    
    # Remover koi_depth "profundidade da curva de luz" (69 VIF)
    # É proporcional e deve ser calculado 
    # usando koi_prad "Raio do Planeta" e koi_srad "Raio da Estrela" e talvez outras variáveis
  )

glm_dummy <- glm(koi_disposition ~ ., data = kepler_scaled_small, family = binomial())
vif_values <- car::vif(glm_dummy)
print(vif_values)

cat("\n--- Remover koi_slogg ---\n")
kepler_scaled_small <- kepler_scaled %>%
  select(
    koi_disposition,
    log_koi_period_z, 
    log_koi_prad_z,
    # log_koi_insol_z,
    # log_koi_depth_z,
    log_koi_srad_z,
    # log_koi_teq_z,
    log_koi_impact_z, 
    log_koi_duration_z,
    koi_steff_z, 
    # koi_slogg_z
    
    # Remover koi_slogg "Gravidade da Estrela" (34 VIF)
    # É proporcional a koi_srad "Raio da Estrela" e koi_steff "Temperatura da Estrela"
  )

glm_dummy <- glm(koi_disposition ~ ., data = kepler_scaled_small, family = binomial())
vif_values <- car::vif(glm_dummy)
print(vif_values)

features <- c(
  "log_koi_period_z", 
  "log_koi_prad_z", 
  "log_koi_srad_z", 
  "log_koi_impact_z", 
  "log_koi_duration_z", 
  "koi_steff_z"
)

formula_todas_interacoes <- as.formula(
  paste("koi_disposition ~ (", paste(features, collapse = " + "), ")^2")
)

glm_triagem <- glm(
  formula_todas_interacoes, 
  data = kepler_scaled_small, 
  family = binomial(link = "logit")
)


interacoes_ranking <- tidy(glm_triagem) %>%
  filter(str_detect(term, ":")) %>%
  arrange(p.value) %>%
  mutate(
    is_significant = p.value < 0.05,
    effect_size = abs(estimate)
  )

cat("\n--- RANKING DAS INTERAÇÕES (Do maior poder explicativo para o menor) ---\n")
print(interacoes_ranking, n = 15)

pares_possiveis <- combn(features, 2, simplify = FALSE)

criar_grafico_interacao <- function(par_vars) {
  var1 <- par_vars[1]
  var2 <- par_vars[2]
  
  ggplot(kepler_scaled_small, aes(x = .data[[var1]], y = .data[[var2]], color = koi_disposition)) +
    geom_point(alpha = 0.15, size = 0.8) +
    geom_density_2d(alpha = 0.9, linewidth = 0.5) +
    scale_color_manual(values = c("CONFIRMED" = "#2c7bb6", "FALSE POSITIVE" = "#d7191c")) +
    labs(
      subtitle = paste(var1, "X", var2),
      x = var1,
      y = var2
    ) +
    theme_minimal() +
    theme(
      plot.subtitle = element_text(size = 9, face = "bold"),
      axis.title = element_text(size = 8),
      legend.position = "none" # Remove a legenda de todos para não poluir
    )
}

lista_graficos <- map(pares_possiveis, criar_grafico_interacao)

painel_amostra <- wrap_plots(lista_graficos[1:6], ncol = 3) + 
  plot_annotation(
    title = "Matriz de Interação (Amostra de 6 pares)",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(painel_amostra)

painel_amostra <- wrap_plots(lista_graficos[7:15], ncol = 3) + 
  plot_annotation(
    title = "2 Matriz de Interação (Amostra de 6 pares)",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(painel_amostra)



#pdf("Matriz_Interacoes_Kepler.pdf", width = 10, height = 8)
#print(wrap_plots(lista_graficos[1:4], ncol = 2))
#print(wrap_plots(lista_graficos[5:8], ncol = 2))
#print(wrap_plots(lista_graficos[9:12], ncol = 2))
#print(wrap_plots(lista_graficos[13:15], ncol = 2))
#dev.off()

sink()
