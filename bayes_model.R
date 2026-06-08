#remove.packages(c("rstan", "StanHeaders", "BH"))
#install.packages("rstanarm")
#install.packages("tidymodels")
library(tidyverse)
library(tidymodels)
library(rstanarm)
library(bayesplot)

planetas_raw <- read_csv('bayes/cumulative.csv')

colunas_para_remover <- c(
  'rowid', 'kepid', 'kepoi_name', 'kepler_name',
  'koi_pdisposition', 'koi_score',
  'koi_fpflag_nt', 'koi_fpflag_ss', 'koi_fpflag_co', 'koi_fpflag_ec',
  'koi_tce_delivname'
)

eps <- 1e-8

planetas <- planetas_raw %>%
  select(-any_of(colunas_para_remover)) %>%
  filter(koi_disposition %in% c("CONFIRMED", "FALSE POSITIVE")) %>%
  drop_na(koi_steff, koi_slogg, koi_srad, koi_period, koi_depth, koi_prad, koi_model_snr) %>%
  mutate(
    koi_disposition = as.factor(koi_disposition),
    koi_period_rel_err = koi_period_err1 / (abs(koi_period) + eps),
    koi_depth_rel_err = koi_depth_err1 / (abs(koi_depth) + eps),
    koi_prad_rel_err = koi_prad_err1 / (abs(koi_prad) + eps),
    koi_depth_snr_weighted = koi_depth * koi_model_snr
  )

dados_pca <- planetas %>% select(koi_steff, koi_slogg, koi_srad)
pca_estelar <- prcomp(dados_pca, center = TRUE, scale. = TRUE)
planetas$stellar_pca_index <- pca_estelar$x[, 1]

preditores <- c('koi_period', 'koi_period_err1', 'koi_time0bk', 'koi_time0bk_err1',
       'koi_impact', 'koi_impact_err1', 'koi_impact_err2', 'koi_duration',
       'koi_duration_err1', 'koi_depth', 'koi_depth_err1', 'koi_prad',
       'koi_prad_err1', 'koi_prad_err2', 'koi_teq', 'koi_insol',
       'koi_insol_err1', 'koi_insol_err2', 'koi_model_snr', 'koi_tce_plnt_num',
       'koi_steff_err1', 'koi_steff_err2', 'koi_slogg_err1', 'koi_slogg_err2',
       'koi_srad_err1', 'koi_srad_err2', 'ra', 'dec', 'koi_kepmag',
       'koi_period_rel_err', 'koi_depth_rel_err', 'koi_prad_rel_err',
       'koi_depth_snr_weighted', 'stellar_pca_index')

planetas <- planetas %>% drop_na(any_of(preditores))

set.seed(42)
split <- initial_split(planetas, prop = 0.70, strata = koi_disposition)
planetas_treino <- training(split)
planetas_teste  <- testing(split)

planetas_treino <- planetas_treino %>%
  mutate(across(all_of(preditores), ~ as.numeric(scale(.))))

planetas_teste <- planetas_teste %>%
  mutate(across(all_of(preditores), ~ as.numeric(scale(.))))

f_modelo <- as.formula(paste("koi_disposition ~", paste(preditores, collapse = " + ")))

modelo_bayesiano <- stan_glm(
  formula = f_modelo, 
  data = planetas_treino, 
  family = binomial(link = "logit"), 
  chains = 4,                        
  iter = 2000,                       
  seed = 42,
  refresh = 500                      
)

print(summary(modelo_bayesiano, digits = 3))

plot(modelo_bayesiano, pars = "beta", ci_level = 0.95)

plot(modelo_bayesiano, plotfun = "trace", pars = "beta")

plot(modelo_bayesiano, plotfun = "acf", pars = "beta")

previsoes_prob <- posterior_epred(modelo_bayesiano, newdata = planetas_teste)

cat("\nProbabilidade média de o primeiro planeta do conjunto de teste ser CONFIRMED:", 
    mean(previsoes_prob[, 1]), "\n")

