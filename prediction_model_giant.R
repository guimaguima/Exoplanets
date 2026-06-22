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
library(matrixStats)

color_scheme_set("mix-blue-red")

planetas_raw <- read_csv('data/cumulative.csv')

preditores_mod <- c(
  'periodo_orbital_dias', 'parametro_impacto', 'duracao_transito_hrs', 
  'profundidade_transito_ppm', 'raio_planetario_terra', 'temp_equilibrio_k', 'insolacao', 
  'snr_modelo', 'num_planetas_sistema', 'ascensao_reta', 'declinacao', 'magnitude_kepler', 
  'temp_estrela_k', 'gravidade_estrela_log', 'raio_estrela_sol'
)

planetas <- planetas_raw %>%
  filter(koi_disposition %in% c("CONFIRMED", "FALSE POSITIVE")) %>%
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
  ) %>%
  drop_na(any_of(preditores_mod)) %>%
  mutate(koi_disposition = factor(koi_disposition, levels = c("FALSE POSITIVE", "CONFIRMED")))

set.seed(42)
split <- initial_split(planetas, prop = 0.70, strata = koi_disposition)
dados_treino <- training(split)
dados_teste  <- testing(split)

f_modelo <- as.formula(paste("koi_disposition ~", paste(preditores_mod, collapse = " + ")))

ajustar_modelo_priori <- function(priori_coefs, priori_interc) {
  stan_glm(
    formula = f_modelo,
    data = dados_treino,
    family = binomial(link = "logit"),
    prior = priori_coefs,
    prior_intercept = priori_interc,
    chains = 4,         
    cores = 4,          
    iter = 4000,
    seed = 42,
    refresh = 500
  )
}

mod_t_original <- ajustar_modelo_priori(
  priori_coefs = student_t(df = 3, location = 0, scale = 2.5, autoscale = TRUE),
  priori_interc = student_t(df = 3, location = 0, scale = 5, autoscale = TRUE)
)

mod_t_vaga <- ajustar_modelo_priori(
  priori_coefs = student_t(df = 3, location = 0, scale = 100, autoscale = TRUE),
  priori_interc = student_t(df = 3, location = 0, scale = 100, autoscale = TRUE)
)

mod_normal <- ajustar_modelo_priori(
  priori_coefs = normal(location = 0, scale = 2.5, autoscale = TRUE),
  priori_interc = normal(location = 0, scale = 5, autoscale = TRUE)
)

mod_laplace <- ajustar_modelo_priori(
  priori_coefs = laplace(location = 0, scale = 2.5, autoscale = TRUE),
  priori_interc = normal(location = 0, scale = 5, autoscale = TRUE)
)

mod_normal_rest <- ajustar_modelo_priori(
  priori_coefs = normal(location = 0, scale = 0.5, autoscale = TRUE),
  priori_interc = normal(location = 0, scale = 5, autoscale = TRUE)
)

mod_cauchy <- ajustar_modelo_priori(
  priori_coefs = cauchy(location = 0, scale = 2.5, autoscale = TRUE),
  priori_interc = cauchy(location = 0, scale = 5, autoscale = TRUE)
)

mod_cauchy_vaga <- ajustar_modelo_priori(
  priori_coefs = cauchy(location = 0, scale = 100, autoscale = TRUE),
  priori_interc = cauchy(location = 0, scale = 100, autoscale = TRUE)
)

# Análises posteriori

plotar_diagnosticos_mcmc <- function(modelo, nome_modelo, preditores, N_por_cadeia = 1000) {
  cat(paste0("\n--- Gerando Gráficos de Diagnóstico MCMC: ", nome_modelo, " ---\n"))
  
  # Trace Plot
  p_trace <- mcmc_trace(modelo, pars = preditores) + 
    ggplot2::ggtitle(paste("Trace Plots -", nome_modelo))
  print(p_trace)
  
  # Autocorrelação (ACF)
  limite_acf <- 1.96 / sqrt(N_por_cadeia)
  p_acf <- mcmc_acf(modelo, pars = preditores) + 
    ggplot2::geom_hline(
      yintercept = c(-limite_acf, limite_acf), 
      linetype = "dashed", color = "#1c73b8", alpha = 0.7, linewidth = 0.8
    ) +
    ggplot2::ggtitle(paste("Função de Autocorrelação (ACF) -", nome_modelo), 
                     subtitle = "Envelope pontilhado azul indica o IC de 95% para ruído branco")
  print(p_acf)
}

calc_moda_continua <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

calcular_moda_discreta <- function(x) {
  tabela_freq <- table(x)
  valor_mais_comum <- names(tabela_freq)[which.max(tabela_freq)]
  return(as.numeric(valor_mais_comum))
}

resumir_posteriori <- function(modelo, nome_modelo, preditores) {
  cat(paste0("\n--- Sumário do Modelo: ", nome_modelo, " ---\n"))
  print(summary(modelo, digits = 3))
  
  cat(paste0("\n--- Intervalos HPD (95%): ", nome_modelo, " ---\n"))
  matriz_amostras <- as.matrix(modelo)[, preditores]
  print(hdi(matriz_amostras, credMass = 0.95))
  
  # Gráficos HPD
  base_plot <- modelo %>%
    gather_draws(!!!syms(preditores)) %>%
    ggplot(aes(y = .variable, x = .value)) +
    stat_halfeye(.width = 0.95, point_interval = mode_hdi, fill = "#1c73b8", alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    theme_minimal()
  
  p_hpd_completo <- base_plot + 
    labs(title = paste("Posterior: Intervalos HPD (95%) -", nome_modelo),
         subtitle = "Linha vermelha no zero indica efeito nulo. Pontos representam a Moda (MAP).",
         x = "Estimativa do Coeficiente", y = "Preditores")
  
  p_hpd_zoom <- base_plot + 
    coord_cartesian(xlim = c(-0.04, 0.04)) +
    labs(title = paste("Posterior: Efeitos Menores (Zoom HPD 95%) -", nome_modelo),
         subtitle = "Variáveis com grandes efeitos cortadas visualmente",
         x = "Estimativa do Coeficiente (Zoom)", y = "Preditores")
  
  print(p_hpd_completo)
  print(p_hpd_zoom)
  
  cat(paste0("\n--- Tabela de Odds Ratio: ", nome_modelo, " ---\n"))
  preditores_intercept <- c(preditores, "(Intercept)")
  tabela_odds_ratio <- modelo %>%
    gather_draws(!!!syms(preditores_intercept)) %>%
    mode_hdi(.width = 0.95) %>% 
    mutate(
      Odds_Ratio = exp(.value),
      HPD_Inferior_OR = exp(.lower),
      HPD_Superior_OR = exp(.upper)
    ) %>%
    select(.variable, Odds_Ratio, HPD_Inferior_OR, HPD_Superior_OR) %>%
    rename(Preditores = .variable) %>%
    arrange(desc(Odds_Ratio)) %>%
    mutate(across(where(is.numeric), ~ format(round(., 4), scientific = FALSE)))
  
  print(as.data.frame(tabela_odds_ratio))
}

avaliar_desempenho_classificacao <- function(modelo, nome_modelo, dados_treino, dados_teste) {
  
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
      Metodo = nome_estimador, Threshold = round(corte_ideal, 4), AUC_Teste = round(auc_teste, 4),
      Acuracia = round(cm$overall["Accuracy"], 4), Sensibilidade = round(cm$byClass["Sensitivity"], 4),
      Especificidade = round(cm$byClass["Specificity"], 4), F1_Score = round(cm$byClass["F1"], 4)
    )
    list(metricas = metricas, matriz = cm$table, roc = roc_teste)
  }
  
  matriz_prob_treino <- posterior_epred(modelo)
  matriz_prob_teste <- posterior_epred(modelo, newdata = dados_teste)
  
  resultado_media <- avaliar_estimador("Média", apply(matriz_prob_treino, 2, mean), apply(matriz_prob_teste, 2, mean), dados_treino$koi_disposition, dados_teste$koi_disposition)
  resultado_mediana <- avaliar_estimador("Mediana", apply(matriz_prob_treino, 2, median), apply(matriz_prob_teste, 2, median), dados_treino$koi_disposition, dados_teste$koi_disposition)
  resultado_moda <- avaliar_estimador("Moda", apply(matriz_prob_treino, 2, calcular_moda_continua), apply(matriz_prob_teste, 2, calcular_moda_continua), dados_treino$koi_disposition, dados_teste$koi_disposition)
  
  cat(paste0("\n--- Comparação de Desempenho no Teste: ", nome_modelo, " ---\n"))
  print(bind_rows(resultado_media$metricas, resultado_mediana$metricas, resultado_moda$metricas) %>% as.data.frame())
  
  plot(resultado_media$roc, col = "#1c73b8", lwd = 2, main = paste("ROC - Teste -", nome_modelo))
  plot(resultado_mediana$roc, col = "#d95f02", lwd = 2, add = TRUE, lty = 2)
  plot(resultado_moda$roc, col = "#1b9e77", lwd = 2, add = TRUE, lty = 3)
  legend("bottomright", 
         legend = c(paste("Média (AUC:", round(auc(resultado_media$roc), 3), ")"),
                    paste("Mediana (AUC:", round(auc(resultado_mediana$roc), 3), ")"),
                    paste("Moda (AUC:", round(auc(resultado_moda$roc), 3), ")")), 
         col = c("#1c73b8", "#d95f02", "#1b9e77"), lwd = 2, lty = c(1, 2, 3))
}

analisar_preditiva_zy <- function(modelo, nome_modelo, dados_teste) {
  cat(paste0("\n--- Z|Y preditiva: ", nome_modelo, " ---\n"))
  
  Z_preditiva <- posterior_predict(modelo, newdata = dados_teste)
  y_real <- as.numeric(dados_teste$koi_disposition) - 1
  
  color_scheme_set("mix-blue-red")
  p_ppc <- ppc_bars(y = y_real, yrep = Z_preditiva) +
    labs(title = paste("Distribuição Preditiva a Posteriori (Z | Y) -", nome_modelo),
         subtitle = "Barras escuras: Dados Reais | Barras claras: Incerteza Simulada (Z)",
         x = "Classe (0 = Falso Positivo, 1 = Confirmado)", y = "Contagem") +
    theme_minimal()
  print(p_ppc)
  
  extrair_estatisticas <- function(contagens, nome_classe) {
    cat(paste("\n---", nome_classe, "---\n"))
    cat("Ponto Preto (Mediana):", median(contagens), "\n")
    cat("Moda:", calcular_moda_discreta(contagens), "\n")
    cat("Média:", mean(contagens), "\n")
    intervalo <- quantile(contagens, probs = c(0.025, 0.975))
    cat("Limite Inferior (2.5%):", intervalo[1], "\n")
    cat("Limite Superior (97.5%):", intervalo[2], "\n")
  }
  
  extrair_estatisticas(rowSums(Z_preditiva == 1), "Classe 1 (CONFIRMED)")
  extrair_estatisticas(rowSums(Z_preditiva == 0), "Classe 0 (FALSE POSITIVE)")
}

avaliar_lpd_teste <- function(modelo, nome_modelo, dados_teste) {
  cat(paste0("\n--- Log Predictive Density (LPD) no Teste: ", nome_modelo, " ---\n"))
  
  # Gera a matriz S x N de log-verossimilhanças no teste
  matriz_log_lik <- log_lik(modelo, newdata = dados_teste)
  
  S <- nrow(matriz_log_lik) # Número de amostras MCMC
  
  # Calcula o log da média das probabilidades para cada observação
  # Usamos logSumExp para estabilidade numérica matemática
  lpd_por_obs <- apply(matriz_log_lik, 2, function(x) logSumExp(x) - log(S))
  
  # Soma total para o modelo
  lpd_total <- sum(lpd_por_obs)
  
  cat("LPD Total:", round(lpd_total, 4), "\n")
  cat("(Valores mais altos / menos negativos são melhores)\n")
  
  return(lpd_total)
}

plotar_diagnosticos_mcmc(mod_t_original, "Modelo t_3(0, 2.5)", preditores_mod)
resumir_posteriori(mod_t_original, "Modelo t_3(0, 2.5)", preditores_mod)
avaliar_desempenho_classificacao(mod_t_original, "Modelo t_3(0, 2.5)", dados_treino, dados_teste)
analisar_preditiva_zy(mod_t_original, "Modelo t_3(0, 2.5)", dados_teste)
avaliar_lpd_teste(mod_t_original, "Modelo t_3(0, 2.5)", dados_teste)

plotar_diagnosticos_mcmc(mod_t_vaga, "Modelo t_3(0, 100)", preditores_mod)
resumir_posteriori(mod_t_vaga, "Modelo t_3(0, 100)", preditores_mod)
avaliar_desempenho_classificacao(mod_t_vaga, "Modelo t_3(0, 100)", dados_treino, dados_teste)
analisar_preditiva_zy(mod_t_vaga, "Modelo t_3(0, 100)", dados_teste)
avaliar_lpd_teste(mod_t_vaga, "Modelo t_3(0, 100)", dados_teste)

plotar_diagnosticos_mcmc(mod_normal, "Modelo Normal", preditores_mod)
resumir_posteriori(mod_normal, "Modelo Normal", preditores_mod)
avaliar_desempenho_classificacao(mod_normal, "Modelo Normal", dados_treino, dados_teste)
analisar_preditiva_zy(mod_normal, "Modelo Normal", dados_teste)
avaliar_lpd_teste(mod_normal, "Modelo Normal", dados_teste)

plotar_diagnosticos_mcmc(mod_laplace, "Modelo Laplace", preditores_mod)
resumir_posteriori(mod_laplace, "Modelo Laplace", preditores_mod)
avaliar_desempenho_classificacao(mod_laplace, "Modelo Laplace", dados_treino, dados_teste)
analisar_preditiva_zy(mod_laplace, "Modelo Laplace", dados_teste)
avaliar_lpd_teste(mod_laplace, "Modelo Laplace", dados_teste)

plotar_diagnosticos_mcmc(mod_normal_rest, "Modelo Normal restritivo (0.5 VAR)", preditores_mod)
resumir_posteriori(mod_normal_rest, "Modelo Normal restritivo (0.5 VAR)", preditores_mod)
avaliar_desempenho_classificacao(mod_normal_rest, "Modelo Normal restritivo (0.5 VAR)", dados_treino, dados_teste)
analisar_preditiva_zy(mod_normal_rest, "Modelo Normal restritivo (0.5 VAR)", dados_teste)
avaliar_lpd_teste(mod_normal_rest, "Modelo Normal restritivo (0.5 VAR)", dados_teste)

plotar_diagnosticos_mcmc(mod_cauchy, "Modelo Cauchy", preditores_mod)
resumir_posteriori(mod_cauchy, "Modelo Cauchy", preditores_mod)
avaliar_desempenho_classificacao(mod_cauchy, "Modelo Cauchy", dados_treino, dados_teste)
analisar_preditiva_zy(mod_cauchy, "Modelo Cauchy", dados_teste)
avaliar_lpd_teste(mod_cauchy, "Modelo Cauchy", dados_teste)

plotar_diagnosticos_mcmc(mod_cauchy_vaga, "Modelo Cauchy Vago", preditores_mod)
resumir_posteriori(mod_cauchy_vaga, "Modelo Cauchy Vago", preditores_mod)
avaliar_desempenho_classificacao(mod_cauchy_vaga, "Modelo Cauchy Vago", dados_treino, dados_teste)
analisar_preditiva_zy(mod_cauchy_vaga, "Modelo Cauchy Vago", dados_teste)
avaliar_lpd_teste(mod_cauchy_vaga, "Modelo Cauchy Vago", dados_teste)