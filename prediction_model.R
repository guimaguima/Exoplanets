library(tidyverse)
library(tidymodels)
library(rstanarm)
library(bayesplot)
library(corrplot)
library(HDInterval)
library(tidybayes)
library(ggplot2)
library(pROC)
library(caret)

color_scheme_set("mix-blue-red")

planetas_raw <- read_csv('data/cumulative.csv')

preditores_mod2 <- c(
  'periodo_orbital_dias', 'parametro_impacto', 'duracao_transito_hrs', 
  'profundidade_transito_ppm', 'raio_planetario_terra', 'temp_equilibrio_k', 
  'snr_modelo', 'num_planetas_sistema', 
  'temp_estrela_k', 'gravidade_estrela_log', 'raio_estrela_sol'
)

planetas <- planetas_raw %>%
  filter(koi_disposition %in% c("CONFIRMED", "FALSE POSITIVE")) %>%
  rename(
    periodo_orbital_dias = koi_period,
    duracao_transito_hrs = koi_duration,
    profundidade_transito_ppm = koi_depth,
    raio_planetario_terra = koi_prad,
    parametro_impacto = koi_impact,
    temp_equilibrio_k = koi_teq,
    snr_modelo = koi_model_snr,
    num_planetas_sistema = koi_tce_plnt_num,
    temp_estrela_k = koi_steff,
    gravidade_estrela_log = koi_slogg,
    raio_estrela_sol = koi_srad
  ) %>%
  drop_na(any_of(preditores_mod2)) %>%
  mutate(koi_disposition = factor(koi_disposition, levels = c("FALSE POSITIVE", "CONFIRMED")))

set.seed(42)
split <- initial_split(planetas, prop = 0.70, strata = koi_disposition)
dados_treino <- training(split)
dados_teste  <- testing(split)

f_modelo2 <- as.formula(paste("koi_disposition ~", paste(preditores_mod2, collapse = " + ")))

set.seed(42)
modelo_bayesiano_2 <- stan_glm(
  formula = f_modelo2,
  data = dados_treino,
  family = binomial(link = "logit"),
  prior = student_t(df = 3, location = 0, scale = 2.5, autoscale = TRUE),
  prior_intercept = student_t(df = 3, location = 0, scale = 5, autoscale = TRUE),
  chains = 4,         
  cores = 4,          
  iter = 4000,
  seed = 42,
  refresh = 500       
)

# Resumo posteriori
mcmc_trace(
  modelo_bayesiano_2, 
  pars = preditores_mod2
) + ggplot2::ggtitle("Trace Plots de Preditores Chave")

N_por_cadeia <- 1000 
limite_acf <- 1.96 / sqrt(N_por_cadeia)

mcmc_acf(
  modelo_bayesiano_2, 
  pars = preditores_mod2
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

matriz_amostras <- as.matrix(modelo_bayesiano_2)[, preditores_mod2]

tabela_hpd <- hdi(matriz_amostras, credMass = 0.95)

modelo_bayesiano_2 %>%
  gather_draws(!!!syms(preditores_mod2)) %>%
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

modelo_bayesiano_2 %>%
  gather_draws(!!!syms(preditores_mod2)) %>%
  ggplot(aes(y = .variable, x = .value)) +
  stat_halfeye(.width = 0.95, point_interval = mode_hdi, 
               fill = "#1c73b8", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  coord_cartesian(xlim = c(-0.04, 0.04)) +
  theme_minimal() +
  labs(
    title = "Posterior: Efeitos Menores (Zoom HPD 95%)",
    subtitle = "Variáveis com grandes efeitos foram cortadas visualmente para permitir a leitura",
    x = "Estimativa do Coeficiente (Zoom)",
    y = "Preditores"
  )

print(summary(modelo_bayesiano_2, digits = 3))

preditores_mod2_intercept <- c(preditores_mod2, "(Intercept)")
tabela_odds_ratio <- modelo_bayesiano_2 %>%
  gather_draws(!!!syms(preditores_mod2_intercept)) %>%
  
  mode_hdi(.width = 0.95) %>% 
  
  mutate(
    Odds_Ratio = exp(.value),
    HPD_Inferior_OR = exp(.lower),
    HPD_Superior_OR = exp(.upper)
  ) %>%
  
  select(.variable, Odds_Ratio, HPD_Inferior_OR, HPD_Superior_OR) %>%
  rename(Preditores = .variable) %>%
  arrange(desc(Odds_Ratio)) %>%
  mutate(
    across(
      where(is.numeric),
      ~ format(round(., 4), scientific = FALSE)
    )
  )

print(as.data.frame(tabela_odds_ratio))

# Teste
matriz_prob_treino <- posterior_epred(modelo_bayesiano_2)
matriz_prob_teste <- posterior_epred(modelo_bayesiano_2, newdata = dados_teste)

calcular_moda_continua <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

probs_treino_media <- apply(matriz_prob_treino, 2, mean)
probs_treino_mediana <- apply(matriz_prob_treino, 2, median)
probs_treino_moda <- apply(matriz_prob_treino, 2, calcular_moda_continua)

probs_teste_media <- apply(matriz_prob_teste, 2, mean)
probs_teste_mediana <- apply(matriz_prob_teste, 2, median)
probs_teste_moda <- apply(matriz_prob_teste, 2, calcular_moda_continua)

avaliar_estimador <- function(nome_estimador, prob_treino, prob_teste, y_treino, y_teste) {
  
  roc_treino <- roc(response = y_treino, predictor = prob_treino, levels = c("FALSE POSITIVE", "CONFIRMED"), quiet = TRUE)
  coords_roc <- coords(roc_treino, "best", ret = c("threshold"), best.method = "youden")
  corte_ideal <- coords_roc$threshold[1]
  
  roc_teste <- roc(response = y_teste, predictor = prob_teste, levels = c("FALSE POSITIVE", "CONFIRMED"), quiet = TRUE)
  auc_teste <- as.numeric(auc(roc_teste))
  
  classificacao <- ifelse(prob_teste >= corte_ideal, "CONFIRMED", "FALSE POSITIVE")
  classificacao <- factor(classificacao, levels = c("FALSE POSITIVE", "CONFIRMED"))
  y_teste_fator <- factor(y_teste, levels = c("FALSE POSITIVE", "CONFIRMED"))
  
  cm <- confusionMatrix(classificacao, y_teste_fator, positive = "CONFIRMED")
  
  metricas <- data.frame(
    Metodo = nome_estimador,
    Threshold = round(corte_ideal, 4),
    AUC_Teste = round(auc_teste, 4),
    Acuracia = round(cm$overall["Accuracy"], 4),
    Sensibilidade = round(cm$byClass["Sensitivity"], 4),
    Especificidade = round(cm$byClass["Specificity"], 4),
    F1_Score = round(cm$byClass["F1"], 4),
    row.names = NULL
  )
  
  list(metricas = metricas, matriz = cm$table, roc = roc_teste)
}

resultado_media <- avaliar_estimador("Média", probs_treino_media, probs_teste_media, dados_treino$koi_disposition, dados_teste$koi_disposition)
resultado_mediana <- avaliar_estimador("Mediana", probs_treino_mediana, probs_teste_mediana, dados_treino$koi_disposition, dados_teste$koi_disposition)
resultado_moda <- avaliar_estimador("Moda", probs_treino_moda, probs_teste_moda, dados_treino$koi_disposition, dados_teste$koi_disposition)

tabela_comparativa <- bind_rows(resultado_media$metricas, resultado_mediana$metricas, resultado_moda$metricas)
print(tabela_comparativa %>% as.data.frame())

cat("\n--- Matriz de Confusão: Média ---\n")
print(resultado_media$matriz)

cat("\n--- Matriz de Confusão: Mediana ---\n")
print(resultado_mediana$matriz)

cat("\n--- Matriz de Confusão: Moda ---\n")
print(resultado_moda$matriz)

plot(resultado_media$roc, col = "#1c73b8", lwd = 2, main = "Comparação ROC - Conjunto de Teste")
plot(resultado_mediana$roc, col = "#d95f02", lwd = 2, add = TRUE, lty = 2)
plot(resultado_moda$roc, col = "#1b9e77", lwd = 2, add = TRUE, lty = 3)

legend("bottomright", 
       legend = c(
         paste("Média (AUC:", round(auc(resultado_media$roc), 3), ")"),
         paste("Mediana (AUC:", round(auc(resultado_mediana$roc), 3), ")"),
         paste("Moda (AUC:", round(auc(resultado_moda$roc), 3), ")")
       ), 
       col = c("#1c73b8", "#d95f02", "#1b9e77"), 
       lwd = 2, lty = c(1, 2, 3))
