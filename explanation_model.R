install.packages("corrplot")
library(tidyverse)
library(rstanarm)
library(bayesplot)
library(corrplot)

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
    koi_disposition = as.factor(koi_disposition),
    koi_depth_snr_weighted = koi_depth * koi_model_snr 
  )

preditores <- c(
  'koi_period', 'koi_time0bk', 'koi_impact', 'koi_duration', 
  'koi_depth', 'koi_prad', 'koi_teq', 'koi_insol', 
  'koi_model_snr', 'koi_tce_plnt_num', 'ra', 'dec', 'koi_kepmag', 
  'koi_steff', 'koi_slogg', 'koi_srad', 
  'koi_depth_snr_weighted'
)

planetas <- planetas %>% 
  drop_na(any_of(preditores)) %>%
  mutate(across(all_of(preditores), ~ as.numeric(scale(.))))

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
  iter = 2000,
  seed = 42,
  refresh = 500       
)



cat("\n--- Resumo da Posteriori ---\n")
print(summary(modelo_bayesiano, digits = 3))

mcmc_intervals(
  as.matrix(modelo_bayesiano), 
  prob = 0.5, prob_outer = 0.95, 
  point_est = "median"
) + ggplot2::ggtitle("Posterior: Intervalos de Confiança (95%) e Efeitos Principais")


mcmc_trace(
  as.matrix(modelo_bayesiano), 
  pars = c("koi_period", "koi_steff", "koi_srad", "koi_prad")
) + ggplot2::ggtitle("Trace Plots de Preditores Chave")


N_por_cadeia<- 1000 
limite_acf <- 1.96 / sqrt(N_por_cadeia)

mcmc_acf(
  as.matrix(modelo_bayesiano), 
  pars = c("koi_period", "koi_depth", "koi_teq")
) + 

  ggplot2::geom_hline(
    yintercept = c(-limite_acf, limite_acf), 
    linetype = "dashed", 
    color = "#1c73b8", 
    alpha = 0.7,
    linewidth = 0.8
  ) +
  ggplot2::ggtitle("Função de Autocorrelação (ACF)", subtitle = "Envelope pontilhado azul indica o Intervalo de Confiança de 95% para ruído branco")

cat("\n Mapa de Calor: Correlação dos Dados (Preditores) ")
matriz_cor_dados <- cor(planetas %>% select(all_of(preditores)), use = "complete.obs")

corrplot(matriz_cor_dados, 
         method = "color", 
         type = "upper", 
         tl.col = "black", 
         tl.cex = 0.7,      
         addCoef.col = "black" 
)

pp_check(modelo_bayesiano, plotfun = "bars") + 
  ggplot2::ggtitle("Posterior Predictive Check (O modelo simula dados parecidos com a base?)")
