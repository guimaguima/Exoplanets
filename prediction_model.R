#remove.packages(c("rstan", "StanHeaders", "BH"))
#install.packages("rstanarm")
#install.packages("tidymodels")
# install.packages(c("pROC", "caret", "tidymodels"))
library(tidyverse)
library(tidymodels)
library(rstanarm)
library(bayesplot)
library(pROC)
library(caret)

color_scheme_set("mix-blue-red")

planetas_raw <- read_csv('bayes/cumulative.csv')

colunas_para_remover <- c(
  'rowid', 'kepid', 'kepoi_name', 'kepler_name',
  'koi_pdisposition', 'koi_score',
  'koi_fpflag_nt', 'koi_fpflag_ss', 'koi_fpflag_co', 'koi_fpflag_ec',
  'koi_tce_delivname'
)

preditores <- c(
  'koi_period', 'koi_time0bk', 'koi_impact', 'koi_duration', 
  'koi_depth', 'koi_prad', 'koi_teq', 'koi_insol', 
  'koi_model_snr', 'koi_tce_plnt_num', 'ra', 'dec', 'koi_kepmag', 
  'koi_steff', 'koi_slogg', 'koi_srad', 
  'koi_depth_snr_weighted'
)


planetas_processados <- planetas_raw %>%
  select(-any_of(colunas_para_remover)) %>%
  filter(koi_disposition %in% c("CONFIRMED", "FALSE POSITIVE", "CANDIDATE")) %>%
  drop_na(any_of(preditores)) %>%
  mutate(
    koi_disposition = factor(koi_disposition, levels = c("FALSE POSITIVE", "CONFIRMED", "CANDIDATE")),
    koi_depth_snr_weighted = koi_depth * koi_model_snr 
  ) %>%
  mutate(across(all_of(preditores), ~ as.numeric(scale(.))))


dados_candidatos <- planetas_processados %>% filter(koi_disposition == "CANDIDATE")

dados_modelagem <- planetas_processados %>% 
  filter(koi_disposition != "CANDIDATE") %>%
  mutate(koi_disposition = droplevels(koi_disposition))

set.seed(42)
split <- initial_split(dados_modelagem, prop = 0.70, strata = koi_disposition)
dados_treino <- training(split)
dados_teste  <- testing(split)

f_modelo <- as.formula(paste("koi_disposition ~", paste(preditores, collapse = " + ")))

set.seed(42)
modelo_preditivo <- stan_glm(
  formula = f_modelo,
  data = dados_treino,
  family = binomial(link = "logit"),
  prior = student_t(df = 3, location = 0, scale = 2.5, autoscale = TRUE),
  prior_intercept = student_t(df = 3, location = 0, scale = 5, autoscale = TRUE),
  chains = 4,         
  cores = 4,          
  iter = 2000,
  seed = 42,
  refresh = 500       
)


cat("\n--- Avaliando Modelo no Conjunto de Teste ---\n")

matriz_prob_teste <- posterior_epred(modelo_preditivo, newdata = dados_teste)
probs_teste_media <- apply(matriz_prob_teste, 2, mean)

roc_obj <- roc(response = dados_teste$koi_disposition, predictor = probs_teste_media, levels = c("FALSE POSITIVE", "CONFIRMED"))

plot(roc_obj, main = "Curva ROC - Conjunto de Teste", col = "#1c73b8", lwd = 2)
auc_val <- auc(roc_obj)
legend("bottomright", legend = paste("AUC =", round(auc_val, 3)), col = "#1c73b8", lwd = 2)

corte_otimo_info <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"), best.method = "youden")
corte_ideal <- corte_otimo_info$threshold[1]

cat(paste("\nCorte Ideal encontrado pela curva ROC:", round(corte_ideal, 4), "\n"))

predicoes_finais_teste <- ifelse(probs_teste_media >= corte_ideal, "CONFIRMED", "FALSE POSITIVE")
predicoes_finais_teste <- factor(predicoes_finais_teste, levels = c("FALSE POSITIVE", "CONFIRMED"))

matriz_confusao <- confusionMatrix(predicoes_finais_teste, dados_teste$koi_disposition, positive = "CONFIRMED")
print(matriz_confusao)


cat("\n--- Realizando predições nos exoplanetas candidatos ---\n")

matriz_prob_cand <- posterior_epred(modelo_preditivo, newdata = dados_candidatos)

resultados_candidatos <- dados_candidatos %>%
  mutate(
    prob_media = apply(matriz_prob_cand, 2, mean),
    hpd_inferior = apply(matriz_prob_cand, 2, quantile, probs = 0.025),
    hpd_superior = apply(matriz_prob_cand, 2, quantile, probs = 0.975),
    
    classificacao_bayesiana = case_when(
      hpd_inferior >= corte_ideal ~ "CONFIRMADO (Alta Certeza)",
      hpd_superior < corte_ideal ~ "FALSO POSITIVO (Alta Certeza)",
      TRUE ~ "INCERTO (Abstenção)"
    )
  ) %>%
  select(koi_disposition, prob_media, hpd_inferior, hpd_superior, classificacao_bayesiana, everything())

cat("\nResumo da Classificação dos Candidatos com Corte de", round(corte_ideal, 3), ":\n")
print(table(resultados_candidatos$classificacao_bayesiana))

resultados_plot <- resultados_candidatos %>% arrange(prob_media) %>% mutate(id = row_number())

ggplot(resultados_plot, aes(x = id, y = prob_media, color = classificacao_bayesiana)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_linerange(aes(ymin = hpd_inferior, ymax = hpd_superior), alpha = 0.2) +
  geom_hline(yintercept = corte_ideal, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    title = "Classificação Bayesiana de Exoplanetas Candidatos",
    subtitle = paste("Intervalo HPD 95% | Linha de Corte Otimizada (ROC):", round(corte_ideal, 3)),
    x = "Candidatos (Ordenados por Probabilidade)",
    y = "Probabilidade de ser Confirmado",
    color = "Classificação"
  ) +
  theme(legend.position = "bottom")