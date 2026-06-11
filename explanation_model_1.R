#install.packages("corrplot")
#install.packages("HDInterval")
#install.packages("tidybayes")
#install.packages("pROC")
#install.packages("caret")
library(tidyverse)
library(rstanarm)
library(bayesplot)
library(corrplot)
library(HDInterval)
library(tidybayes)
library(ggplot2)
library(pROC)
library(caret)

cat("\n============================\n")
cat("\n--- Modelo 1: Atributos Físicos Totais ---\n")
cat("\n============================\n")

color_scheme_set("mix-blue-red")

planetas_raw <- read_csv('bayes/cumulative.csv')

colunas_para_remover <- c(
  'rowid', 'kepid', 'kepoi_name', 'kepler_name',
  'koi_pdisposition', 'koi_score',
  'koi_fpflag_nt', 'koi_fpflag_ss', 'koi_fpflag_co', 'koi_fpflag_ec',
  'koi_tce_delivname'
)

planetas <- planetas_raw %>%
  select(-any_of(colunas_para_remover)) %>%
  filter(koi_disposition %in% c("CONFIRMED", "FALSE POSITIVE")) %>%
  drop_na(koi_steff, koi_slogg, koi_srad, koi_period, koi_depth, koi_prad, koi_model_snr) %>%
  mutate(
    koi_disposition = as.factor(koi_disposition)
  )



preditores <- c(
  'periodo_orbital_dias', 'parametro_impacto', 'duracao_transito_hrs', 
  'profundidade_transito_ppm', 'raio_planetario_terra', 'temp_equilibrio_k', 'insolacao', 
  'snr_modelo', 'num_planetas_sistema', 'ascensao_reta', 'declinacao', 'magnitude_kepler', 
  'temp_estrela_k', 'gravidade_estrela_log', 'raio_estrela_sol'
)

planetas <- planetas %>% 
  drop_na(any_of(preditores))


planetas <- planetas %>%
  rename(
    periodo_orbital_dias = koi_period,
    insolacao = koi_insol,
    duracao_transito_hrs = koi_duration,
    profundidade_transito_ppm = koi_depth,
    raio_planetario_terra = koi_prad,
    parametro_impacto = koi_impact,
    temp_equilibrio_k = koi_teq,
    snr_modelo = koi_model_snr,
    num_planetas_sistema = koi_tce_plnt_num,
    magnitude_kepler = koi_kepmag,
    temp_estrela_k = koi_steff,
    gravidade_estrela_log = koi_slogg,
    raio_estrela_sol = koi_srad,
    ascensao_reta = ra,
    declinacao = dec
  )%>% 
  drop_na(any_of(preditores))


cat("\n Mapa de Calor: Correlação dos Dados (Preditores) ")
matriz_cor_dados <- cor(planetas %>% select(all_of(preditores)), use = "complete.obs")

corrplot(matriz_cor_dados, 
         method = "color", 
         type = "full", 
         tl.col = "black", 
         tl.cex = 0.7,      
         addCoef.col = "black" 
)



f_modelo <- as.formula(paste("koi_disposition ~", paste(preditores, collapse = " + ")))

set.seed(42)
modelo_bayesiano <- stan_glm(
  formula = f_modelo,
  data = planetas,
  family = binomial(link = "logit"),
  prior = student_t(df = 3, location = 0, scale = 2.5, autoscale = TRUE),
  prior_intercept = student_t(df = 3, location = 0, scale = 5, autoscale = TRUE),
  chains = 4,         
  cores = 4,          
  iter = 4000,
  seed = 42,
  refresh = 500       
)


cat("\n--- Resumo da Posteriori ---\n")

mcmc_trace(
  modelo_bayesiano, 
  pars = preditores
) + ggplot2::ggtitle("Trace Plots de Preditores Chave")


N_por_cadeia <- 1000 
limite_acf <- 1.96 / sqrt(N_por_cadeia)

mcmc_acf(
  modelo_bayesiano, 
  pars = preditores
) + 
  ggplot2::geom_hline(
    yintercept = c(-limite_acf, limite_acf), 
    linetype = "dashed", 
    color = "#1c73b8", 
    alpha = 0.7,
    linewidth = 0.8
  ) +
  ggplot2::ggtitle("Função de Autocorrelação (ACF)", subtitle = "Envelope pontilhado azul indica o Intervalo de Confiança de 95% para ruído branco")
cat("\n--- Sumário e Intervalos ---\n")

cat("\n--- Intervalos de Maior Densidade a Posteriori (HPD - 95%) ---\n")

matriz_amostras <- as.matrix(modelo_bayesiano)[, preditores]

tabela_hpd <- hdi(matriz_amostras, credMass = 0.95)

modelo_bayesiano %>%
  gather_draws(!!!syms(preditores)) %>%
  ggplot(aes(y = .variable, x = .value)) +
  stat_halfeye(.width = 0.95, point_interval = mode_hdi, 
               fill = "#1c73b8", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    title = "Posterior: Intervalos HPD (95%)",
    subtitle = "Linha vermelha no zero indica efeito nulo. Pontos representam a Moda (MAP).",
    x = "Estimativa do Coeficiente",
    y = "Preditores"
  )

modelo_bayesiano %>%
  gather_draws(!!!syms(preditores)) %>%
  ggplot(aes(y = .variable, x = .value)) +
  stat_halfeye(.width = 0.95, point_interval = mode_hdi, 
               fill = "#1c73b8", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  coord_cartesian(xlim = c(-0.12, 0.12)) +
  theme_minimal() +
  labs(
    title = "Posterior: Efeitos Menores (Zoom HPD 95%)",
    subtitle = "Variáveis com grandes efeitos foram cortadas visualmente para permitir a leitura",
    x = "Estimativa do Coeficiente (Zoom)",
    y = "Preditores"
  )


print(summary(modelo_bayesiano, digits = 3))

tabela_odds_ratio <- modelo_bayesiano %>%
  gather_draws(!!!syms(preditores)) %>%
  
  mode_hdi(.width = 0.95) %>% 
  
  mutate(
    Odds_Ratio = exp(.value),
    HPD_Inferior_OR = exp(.lower),
    HPD_Superior_OR = exp(.upper)
  ) %>%
  
  select(.variable, Odds_Ratio, HPD_Inferior_OR, HPD_Superior_OR) %>%
  rename(Preditores = .variable) %>%
  arrange(desc(Odds_Ratio)) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

print(as.data.frame(tabela_odds_ratio))


cat("\n====================================================================\n")
cat(" AVALIAÇÃO DE DESEMPENHO PREDITIVO (USANDO A MODA / MAP) \n")
cat("====================================================================\n")

matriz_probabilidades <- posterior_epred(modelo_bayesiano)

calcular_moda_continua <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

probs_preditas_moda <- apply(matriz_probabilidades, 2, calcular_moda_continua)

niveis_roc <- c("FALSE POSITIVE", "CONFIRMED")
roc_obj <- roc(
  response = planetas$koi_disposition, 
  predictor = probs_preditas_moda, 
  levels = niveis_roc
)

plot(roc_obj, main = "Curva ROC - Desempenho na Base de Treino (Moda)", 
     col = "#1c73b8", lwd = 3, print.auc = TRUE, print.auc.y = 0.4)

corte_otimo_info <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"), best.method = "youden")
corte_ideal <- corte_otimo_info$threshold[1]

cat(paste("\n-> Ponto de Corte (Threshold) Ótimo encontrado pela ROC:", round(corte_ideal, 4), "\n"))

classificacao_final <- ifelse(probs_preditas_moda >= corte_ideal, "CONFIRMED", "FALSE POSITIVE")
classificacao_final <- factor(classificacao_final, levels = c("CONFIRMED", "FALSE POSITIVE"))

referencia_real <- factor(planetas$koi_disposition, levels = c("CONFIRMED", "FALSE POSITIVE"))

matriz_confusao <- confusionMatrix(classificacao_final, referencia_real, positive = "CONFIRMED")

cat("\n--- Matriz de Confusão ---\n")
print(matriz_confusao$table)

cat("\n--- Métricas Principais ---\n")
cat(sprintf("Acurácia Geral:    %.2f%%\n", matriz_confusao$overall["Accuracy"] * 100))
cat(sprintf("Sensibilidade:     %.2f%% (Taxa de acerto nos Confirmados)\n", matriz_confusao$byClass["Sensitivity"] * 100))
cat(sprintf("Especificidade:    %.2f%% (Taxa de acerto nos Falsos Positivos)\n", matriz_confusao$byClass["Specificity"] * 100))
cat(sprintf("F1-Score:          %.4f\n", matriz_confusao$byClass["F1"]))

