library(data.table)
library(tidyverse)
library(mice)
library(lfe)
library(sandwich)
library(lmtest)
library(stargazer)


"%ni%" <- Negate("%in%")
se <- function(x) { sd(x)/sqrt(length(x)) }

# wd to the source root of the section
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path))) 

# Load the data for imputation
clearance <- readRDS('data/cctv_crimes_data_balanced.RDS')

clearance_full <- readRDS('data/cctv_crimes_data_unbalanced.RDS')
setDT(clearance_full)
clearance_full[,population:=population/1000]


# Make reg fun
make_reg_tab_w_imp <- function(dat, imp_formulas, reg_formula, m = 20) {
  
  # Imputation
  imp <- mice(
    data = dat,
    m = m, method = "pmm", # predictive mean matching
    maxit = 50,                
    seed = 42,
    formulas = imp_formulas
  )
  
  # Check
  # plot(imp)
  # stripplot(imp, pch = 20, cex = 1.2)
  # densityplot(imp)
  
  # Make reg
  mice_mod <- lapply(1:m, function(i) {
    lm(reg_formula,
       data = mice::complete(data = imp, i))
  })
  
  pool.mod <- mice::pool(mice_mod)
  smry <- summary(pool.mod)
  
  b <- smry$estimate
  names(b) <- smry$term
  se <- smry$std.error
  t_stat = smry$statistic
  p_val <- smry$p.value
  
  r2 <- round(pool.r.squared(pool.mod),2)
  r2 <- paste0(r2[1],' (',round(r2[3]-r2[1],2),')')

  aic <-  paste0(
    round(mean(sapply(mice_mod, AIC))),
    ' (', round(se(sapply(mice_mod, AIC)),2),')'
  )
  
  N <- nrow(dat)
  
  return(
    list(
      mice_mod[[1]],
      b,
      se,
      t_stat,
      p_val,
      r2,
      N,
      aic,
      imp
    )
  )
}

make_reg_tab <- function(dat, formula, clust_err = T) {
  
  result <- lm(formula, 
               data = dat)
  
  res_summ <- summary(result)
  
  N <- length(result$fitted.values)
  r2 <- round(res_summ[["r.squared"]],2)
  r2_adj <- round(res_summ[["adj.r.squared"]],2)
  resid_std_er <- paste0(
    round(res_summ[["sigma"]],2),' (df = ',
    res_summ[["fstatistic"]][["dendf"]],')'
  )
  f_stat <- paste0(
    round(res_summ[["fstatistic"]][["value"]],2),'*** (df = ',
    res_summ[["fstatistic"]][["numdf"]],'; ',
    res_summ[["fstatistic"]][["dendf"]],')'
  )
  
  if (clust_err) {
    result_orig <- result
    result <- coeftest(result, 
                       vcov = vcovCL, 
                       type = 'HC0', 
                       cluster = dat[['district']])
  } 
  
  result_df <- tidy(result)

  b <- result_df$estimate
  names(b) <- result_df$term
  se <- result_df$std.error
  t_stat <- result_df$statistic
  p_val <- result_df$p.value
  
  r2 <- round(summary(result_orig)$r.squared,2)
  aic <- round(AIC(result_orig))
  N <- nrow(dat)
  
  imp <- NA
  
  return(
    list(
      result,
      b,
      se,
      t_stat,
      p_val,
      r2,
      N,
      aic,
      imp
    )
  )
}

# ----------------------------------------------------------------------

# Clearance in general 
imputation_controls_full <- c('district','AO','year2020_plus', "population","ln_rent_price","year",'murder_rate_2013','reg_crime_per_1k_pop',
                         'clearance_percent', 
                         'cctv_public_per_1k_pop','cctv_courtyard_per_1k_pop','cctv_ent_per_1k_pop',
                         'cctv_public',"cctv_courtyard","cctv_entrance")

df <- clearance[,
                c(imputation_controls_full), with = F]
df <- df[,
         district := as.factor(district)]
df <- df[,
         AO := as.factor(AO)]
df <- df[,
         year := as.factor(year)]

# ---------------------------------------------------------
# Balanced imputed panel

m1.1 <- make_reg_tab_w_imp(
  dat = df, 
  imp_formulas = list(
    as.formula('clearance_percent ~ year + cctv_public + cctv_courtyard + cctv_entrance + population + ln_rent_price + district'),
    as.formula('reg_crime_per_1k_pop ~ year + cctv_public + cctv_courtyard + cctv_entrance + population + ln_rent_price + district')
  ),
  reg_formula = 'clearance_percent ~ (year2020_plus*cctv_public_per_1k_pop) + population + ln_rent_price + reg_crime_per_1k_pop + murder_rate_2013 + district'
)

m1.2 <- make_reg_tab_w_imp(
  dat = df, 
  imp_formulas = list(
    as.formula('clearance_percent ~ year + cctv_public + cctv_courtyard + cctv_entrance + population + ln_rent_price + district'),
    as.formula('reg_crime_per_1k_pop ~ year + cctv_public + cctv_courtyard + cctv_entrance + population + ln_rent_price + district')
  ), 
  reg_formula = 'clearance_percent ~ (year2020_plus*cctv_courtyard_per_1k_pop) + population + ln_rent_price + reg_crime_per_1k_pop + murder_rate_2013 + district'
)

m1.3 <- make_reg_tab_w_imp(
  dat = df, 
  imp_formulas = list(
    as.formula('clearance_percent ~ year + cctv_public + cctv_courtyard + cctv_entrance + population + ln_rent_price + district'),
    as.formula('reg_crime_per_1k_pop ~ year + cctv_public + cctv_courtyard + cctv_entrance + population + ln_rent_price + district')
  ), 
  reg_formula = 'clearance_percent ~ (year2020_plus*cctv_ent_per_1k_pop) + population + ln_rent_price + reg_crime_per_1k_pop + murder_rate_2013 + district'
)

# ----------------------------------------------------------------------
# Save imputed data for clearance_percent and reg_crime_per_1k_pop to the main df
df_filled <- df

df_filled[is.na(reg_crime_per_1k_pop),
          reg_crime_per_1k_pop := rowMeans(m1.1[[9]][["imp"]][["reg_crime_per_1k_pop"]])]

df_filled[,clearance_percent := clearance_percent]
df_filled[is.na(clearance_percent),
          clearance_percent := rowMeans(m1.1[[9]][["imp"]][["clearance_percent"]])]

saveRDS(df_filled,
        'data/produced/cctv_crimes_data_imputed.RDS')
# ----------------------------------------------------------------------

df_filled[,.N,by=district]
df_filled[,.(mean=mean(clearance_percent),
             se = se(clearance_percent)),
          by=year2020_plus]
df_filled[,.(mean=mean(reg_crime_per_1k_pop),
             se = se(reg_crime_per_1k_pop)),
          by=year2020_plus]



# ----------------------------------------------------------------------

# Thefts (unbalanced)

# Compute thefts reg
m2.1 <- make_reg_tab(formula = 'clearance_thefts_percent ~ (year2020_plus*cctv_public_per_1k_pop) + population + ln_rent_price + reg_crime_per_1k_pop + murder_rate_2013 + district',
                       dat = clearance_full)


m2.2 <- make_reg_tab(formula = 'clearance_thefts_percent ~ (year2020_plus*cctv_courtyard_per_1k_pop) + population +  ln_rent_price + reg_crime_per_1k_pop + murder_rate_2013 + district',
                       dat = clearance_full)

m2.3 <- make_reg_tab(formula = 'clearance_thefts_percent ~ (year2020_plus*cctv_ent_per_1k_pop) + population + ln_rent_price + reg_crime_per_1k_pop + murder_rate_2013 + district',
                       dat = clearance_full)
  

# ----------------------------------------------------------------------
# Pub. spaces


m3.1 <- make_reg_tab(formula = 'clearance_pub_sp_percent ~ (year2020_plus*cctv_public_per_1k_pop) + population + ln_rent_price + reg_crime_per_1k_pop + murder_rate_2013 + district',
                       dat = clearance_full)

m3.2 <- make_reg_tab(formula = 'clearance_pub_sp_percent ~ (year2020_plus*cctv_courtyard_per_1k_pop) + population +  ln_rent_price + reg_crime_per_1k_pop + murder_rate_2013 + district',
                       dat = clearance_full)

m3.3 <- make_reg_tab(formula = 'clearance_pub_sp_percent ~ (year2020_plus*cctv_ent_per_1k_pop) + population + ln_rent_price + reg_crime_per_1k_pop + murder_rate_2013 + district',
                       dat = clearance_full)


# ----------------------------------------------------------------------

m_pack <- list(m1.1,m1.2,m1.3,
               m2.1,m2.2,m2.3,
               m3.1,m3.2,m3.3)


stargazer(
  lapply(m_pack, function(x){x[[1]]}),
  coef = lapply(m_pack, function(x){x[[2]]}),
  se = lapply(m_pack, function(x){x[[3]]}),
  t = lapply(m_pack, function(x){x[[4]]}),
  p = lapply(m_pack, function(x){x[[5]]}),
  type = 'text',
  out = 'output/m1_clearance_rate_imp.doc',
  title = 'Clearance Rate (balanced)', 
  dep.var.caption = 'Clearance rate',
  digits = 2,
  omit.stat=c("all"),
  omit = c('district'),
  omit.labels = c('District'),
  omit.yes.no = c('Yes','No'),
  add.lines=
    list(
      c("Observations", sapply(m_pack, function(x){x[[7]]})),
      c("R2", sapply(m_pack, function(x){x[[6]]})),
      c("AIC", sapply(m_pack, function(x){x[[8]]}))
      )
)

stargazer(
  lapply(m_pack, function(x){x[[1]]}),
  coef = lapply(m_pack, function(x){x[[2]]}),
  se = lapply(m_pack, function(x){x[[3]]}),
  t = lapply(m_pack, function(x){x[[4]]}),
  p = lapply(m_pack, function(x){x[[5]]}),
  type = 'latex',
  out = 'output/m1_clearance_rate_imp.tex',
  title = 'Clearance Rate (balanced)', 
  dep.var.caption = 'Clearance rate',
  digits = 2,
  omit.stat=c("all"),
  omit = c('district'),
  omit.labels = c('District'),
  omit.yes.no = c('Yes','No'),
  add.lines=
    list(
      c("Observations", sapply(m_pack, function(x){x[[7]]})),
      c("R2", sapply(m_pack, function(x){x[[6]]})),
      c("AIC", sapply(m_pack, function(x){x[[8]]}))
    ))


# ----------------------------------------------------------------------

