library(tidyverse)
library(data.table)
library(stringr)
library(stringi)
library(sf)
library(mapview)

library(MASS)
library(AER)

library(stargazer)

sf_use_s2(F)

"%ni%" <- Negate("%in%")

# wd to the source root of the section
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path))) 

grid_cctv_prot_long_w_lag <- readRDS('data/cctv_protests_data_long.RDS')

grid_cctv_prot_long_w_lag <- grid_cctv_prot_long_w_lag %>% 
  filter(year %in% 2018:2021) %>% 
  filter(ADMIN_L5 == 'Центральный административный округ') %>% 
  filter(!is.na(reg_crime_per_1k_pop_lag1))


grid_cctv_prot_long_w_lag <- grid_cctv_prot_long_w_lag %>% 
  mutate(
    protests_before_2y_only = protests_before_2y - protests_before_1y,
    protests_after_2y_only = protests_after_2y - protests_after_1y
  )

# Regression function
make_reg_tab <- function(dat, formula, nb = T, error_cluserization = T) {
  
  if (nb) {
    mod <- glm.nb(formula, 
                  data = dat,
                  control=glm.control(maxit=150))
  } else {
    mod <- glm(formula, 
               data = dat,
               family = binomial(link = cloglog))
  }
  
  result_log_l <- round(logLik(mod))
  result_AIC <- round(AIC(mod))
  result_N <- length(mod$fitted.values)
  
  if (error_cluserization) {
    result_coef <- exp(coef(mod))
    
    result_vcov <- vcovCL(mod, 
                          cluster = grid_cctv_prot_long_w_lag$district, 
                          type = "HC0")
    result_mod <- coeftest(mod, vcov = result_vcov)
    result_se <- sqrt(diag(result_vcov))
    
  } else {
    result_mod <- mod
    
    beta <- coef(mod)
    se <- sqrt(diag(vcov(mod)))
    
    result_coef <- exp(beta)
    result_se   <- exp(beta) * se 
    
  }
  
  results <- list(mod,
                  result_mod,
                  result_coef,
                  result_se,
                  result_log_l, 
                  result_AIC, 
                  result_N)
  
  return(results)
}

# ----------------------


# ----------------------
# Compute the regression
# Spatial lag hex reg

# Negative binomial. Whole city

# 0. Protest (dummy) before 1 year of cctv (count) installation in a hex. Naive
# m0 <- make_reg_tab(grid_cctv_prot_long_w_lag, 
#                    'cctv ~ protests_before_1y_dummy')

m1 <- make_reg_tab(dat = grid_cctv_prot_long_w_lag, 
                   formula = 'cctv ~ protests_before_2y_only + as.factor(date_w_parts) +  cctv_ver2017_dummy + reg_crime_per_1k_pop_lag1 + lag + population + ln_rent_price + district')


# 2. Protest (dummy) before cctv (count) installation in a hex. With controls 
m2 <- make_reg_tab(grid_cctv_prot_long_w_lag, 
                   'cctv ~ protests_before_1y + as.factor(date_w_parts) +  cctv_ver2017_dummy + reg_crime_per_1k_pop_lag1 + lag + population + ln_rent_price + district')


# 3. protests_after_1y_dummy. With controls 
m3 <- make_reg_tab(grid_cctv_prot_long_w_lag, 
                   'cctv ~ protests_after_1y + as.factor(date_w_parts) +  cctv_ver2017_dummy + reg_crime_per_1k_pop_lag1 + lag + population + ln_rent_price + district')

# 4. protests_after_2y_dummy. With controls 
m4 <- make_reg_tab(grid_cctv_prot_long_w_lag, 
                   'cctv ~ protests_after_2y_only + as.factor(date_w_parts) +  cctv_ver2017_dummy + reg_crime_per_1k_pop_lag1 + lag + population + ln_rent_price + district')

m_pack <- list(m1,m2,m3,m4)


stargazer(
  lapply(m_pack, function(x){x[[2]]}),
  coef = lapply(m_pack, function(x){x[[3]]}),
  se = lapply(m_pack, function(x){x[[4]]}),
  p = lapply(m_pack, function(x){x[[2]][,4]}),
  type = 'html',
  out = 'output/reg_results/m2_NegBin_wCrime_CAO_FIN.doc',
  title = 'Negative Binomial Regression. CAO', 
  dep.var.caption = 'Number of New Public Place CCTV',
  digits = 2,
  covariate.labels = c("Protests before (2 year, dummy)", "Protests before (1 year, dummy)", "Protests after (1 year, dummy)", "Protests after (2 year, dummy)", 'CCTV established by 2017 (dummy)', 'Crimes before 1 year in a district'),
  omit.stat=c("all"),
  omit = c('date_w_parts','district', 'population', 'rent_price'),
  omit.labels = c('Half-year', 'district', 'District population', 'District rent price'),
  omit.yes.no = c('Yes','No'),
  add.lines=
    list(
      c("Observations", sapply(m_pack, function(x){x[[6]]})),
      c("AIC", sapply(m_pack, function(x){x[[5]]})),
      c("Log Likelihood", round(sapply(m_pack, function(x){x[[2]]}),2)))
)



stargazer(
  lapply(m_pack, function(x){x[[2]]}),
  coef = lapply(m_pack, function(x){x[[3]]}),
  se = lapply(m_pack, function(x){x[[4]]}),
  p = lapply(m_pack, function(x){x[[2]][,4]}),
  type = 'latex',
  out = 'output/reg_results/m2_NegBin_wCrime_CAO_FIN.tex',
  title = 'Negative Binomial Regression. CAO', 
  dep.var.caption = 'Number of New Public Place CCTV',
  digits = 2,
  covariate.labels = c("Protests before (2 year, dummy)", "Protests before (1 year, dummy)", "Protests after (1 year, dummy)", "Protests after (2 year, dummy)", 'CCTV established by 2017 (dummy)', 'Crimes before 1 year in a district'),
  omit.stat=c("all"),
  omit = c('date_w_parts','district', 'population', 'rent_price'),
  omit.labels = c('Half-year', 'district', 'District population', 'District rent price'),
  omit.yes.no = c('Yes','No'),
  add.lines=
    list(
      c("Observations", sapply(m_pack, function(x){x[[6]]})),
      c("AIC", sapply(m_pack, function(x){x[[5]]})),
      c("Log Likelihood", round(sapply(m_pack, function(x){x[[2]]}),2)))
)

