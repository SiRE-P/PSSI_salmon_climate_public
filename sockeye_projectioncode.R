###Projection code for multi-stage salmon model###
###Works with cmdstanr::CmdStanMCMC or rstan::stanfit###

library(dplyr)
library(ggplot2)
library(tidyr)
library(tidybayes)
library(stringr)
library(boot)

# library(ggdist)
# library(stringr)
# library(purrr)
# library(posterior)

# --- NOTES for running projection model --- #

# Use Median spawner abundance from 1991-2020 period for projections
# standardized hindcast data is for the relative to period across the data that is available (~1981-2024?)
# projection data is standardized to the original hindcast data
# Use "..._complete" from model output data, includes raw and imputed data
# Use parameter coefficients from model output 
# run each year of 1991-2020 for base period to compare to
# Then use Jan's standardized projection data to run future scenario period 

# Use coefficients from model output each year of base scenario (1991-2020) and projection period (2041-2070)
# keep Osoyoos and wenatchee discharge data set to no effect (zero)

# uncertainty:
#  parameter - interations from model (e.g. b_ coefficients for env effects)
#  annual variability
#  scenario - SSP4.5 v SSP 8.5


##### FUNCTIONS #####

# --- extract stan model fit --- #

# extractor (draws-first arrays)
extract_draws <- function(fit, variable_list = "all") {
  out <- list()
  
  if (inherits(fit, "stanfit")) {
    # rstan path
    if(length(variable_list) == 1 && "all" %in% variable_list){
      variable_list <- names(rstan::extract(fm, permuted = TRUE))
    }
    avail <- names(rstan::extract(fit, permuted = TRUE))
    have  <- intersect(variable_list, avail)
    miss  <- setdiff(variable_list, avail)
    if (length(have)) out <- rstan::extract(fit, pars = have)
  }
  
  if (inherits(fit, "CmdStanMCMC")) {
    # cmdstanr path
    if(length(variable_list) == 1 && "all" %in% variable_list){
      variable_list <- fit$metadata()$stan_variables
    }
    da <- fit$draws(variables = variable_list, format = "draws_array")
    rvs <- posterior::as_draws_rvars(da)
    have <- intersect(variable_list, names(rvs))
    for (v in have) {out[[v]] <- posterior::draws_of(rvs[[v]])} 
    miss <- setdiff(variable_list, names(out))
    
  }
  
  if (length(miss)) message("Not saved in fit (skipping): ", paste(miss, collapse = ", "))
  return(out)
}



# --- output projection function --- #

## will need to change indexing if we're doing closed loop and fishing included

proj1 <- function(N_pop, years, N_fw_ages, N_mar_ages, N_iter, spawner_abd, post, ann.dev = NULL){
  
  
  ##### NEW MODEL HERE
  
  # --- model parameters --- #
  
  # dimension names
  iter <- 1:N_iter
  year <- years
  pop <- cuid_dat$pop_name_short
  age_fw <- c("2", "3")
  age <- c("1.1", "1.2", "1.3", "2.1", "2.2", "2.3")
  
  arr1 <- array(0,dim = c(iter = N_iter, year = length(years), pop = N_pop), dimnames = list(iter=iter,year=year,pop=pop))
  arrfw <- array(0,dim = c(iter = N_iter, year = length(years), age_fw = N_fw_ages, pop = N_pop), dimnames = list(iter=iter,year=year,age_fw=age_fw,pop=pop))
  arrage <- array(0,dim = c(iter = N_iter, year = length(years), age = N_fw_ages * N_mar_ages, pop = N_pop), dimnames = list(iter=iter,year=year,age=age,pop=pop))
  
  
  # BH model parameters - naming dimensions
  a_fry <- arr1
  b_fry <- arr1
  
  holdover <- arr1
  age2_survival <- arr1
  a_returns <- arrage
  #b_returns <- arrage
  
  # population parameters
  spawners <- arr1
  fry <- arr1
  smolts_age <- arrfw
  holdover_fry <- arr1
  presmolts <- arr1
  smolts <- arr1
  smolts_by <- arr1
  returns_age <- arrage
  returns <- arr1
  returns_by <- arr1
  
  # process deviations - added to mu values
  if (!is.null(ann.dev)){
    tdev.years <- ann.dev - 1981 + 1
    
    ta_fry_z_mean <- apply(post$a_fry_log_year_z[,,tdev.years], c(1,2), mean, na.rm = TRUE)
    ta_fry_z_sd <- apply(post$a_fry_log_year_z[,,tdev.years], c(1,2), sd, na.rm = TRUE)
    
    tholdover_z_mean <- apply(post$holdover_logit_z[,tdev.years,], c(1,2), mean, na.rm = TRUE)
    tholdover_z_sd <- apply(post$holdover_logit_z[,tdev.years,], c(1,2), sd, na.rm = TRUE)
    
    tage2sur_z_mean <- apply(post$age2_survival_logit_z[,tdev.years,], c(1,3), mean, na.rm = TRUE)
    tage2sur_fry_z_sd <- apply(post$age2_survival_logit_z[,tdev.years,], c(1,3), sd, na.rm = TRUE)
    
    ta_ret_z_mean <- apply(post$a_returns_logit_year_z[,,,,tdev.years], c(1,2,3,4), mean, na.rm = TRUE)
    ta_ret_z_sd <- apply(post$a_returns_logit_year_z[,,,,tdev.years], c(1,2,3,4), sd, na.rm = TRUE)
  }
  
  # model
  for(p in 1:N_pop){
    print(paste0("pop = ", p))
    
    reg <- cuid_datx$major_watershed_id[which(cuid_datx$cu_id == p)]
    
    # --- Env data and BH parameters --- #
    
    ## smolts
    for(y in 1:length(years)){
      tyear <- years[y]
      
      if((y - 2) >= 1){
        # smolt_age2 cov
        tFTwre1 <- FTwre$value[which(FTwre$year == (tyear - 2) & FTwre$p == p)] # incubation winter
        tFTsre1 <- FTsre$value[which(FTsre$year == (tyear - 1) & FTsre$p == p)] # first summer
        tFTwre2 <- FTwre$value[which(FTwre$year == (tyear - 1) & FTwre$p == p)] # second winter
        
        # smolt_age3 cov
        tFTsre2 <- FTsre$value[which(FTsre$year == tyear & FTsre$p == p)] # second summer
        tFTwre3 <- FTwre$value[which(FTwre$year == tyear & FTwre$p == p)] # third winter
        
      } else {
        tFTwre1 <- NA # 0 or NAs instead? Or maybe the lowest year of values to get returns for the first years
        tFTsre1 <- NA
        tFTwre2 <- NA
        
        tFTsre2 <- NA
        tFTwre3 <- NA
      }
      
      # a_smolts setup
      log_fw <- post$a_fry_log_year_mu[,p] + post$a_fry_log_year_sigma[,p] * 
        (if(!is.null(ann.dev)) rnorm(N_iter, mean = ta_fry_z_mean, sd = ta_fry_z_sd) else 0)
      # log_fw <- (post$a_fry_log_year_mu[,p])
      log_fw <- log_fw + post$b_FTwre[,p,1] * tFTwre1 +
        post$b_FTsre[,p,1] * tFTsre1 +
        post$b_FTwre[,p,2] * tFTwre2
      a_fry[,y,p] <- exp(log_fw)
      
      b_fry[,y,p] <- post$b_fry[,p]
      
      # fry to smolts estimation
      # holdover <- inv.logit(post$holdover_logit_mu[,p] + post$holdover_logit_sigma[,p])
      holdover[,y,p] <- inv.logit(post$holdover_logit_mu[,p] + post$holdover_logit_sigma[,p] * 
                                    (if(!is.null(ann.dev)) rnorm(N_iter, mean = tholdover_z_mean, sd = tholdover_z_sd) else 0))
      
      # age2_survival <- inv.logit(post$age2_survival_logit_mu[,p] + post$age2_survival_logit_sigma[,p] +
      age2_survival[,y,p] <- inv.logit(post$age2_survival_logit_mu[,p] + post$age2_survival_logit_sigma[,p] * 
                                         (if(!is.null(ann.dev)) rnorm(N_iter, mean = tage2sur_z_mean, sd = tage2sur_fry_z_sd) else 0) +
                                       post$b_FTsre[,p,2] * tFTsre2 + 
                                       post$b_FTwre[,p,3] * tFTwre3)
    }
    
    ## returns
    for(y in 1:length(years)){
      tyear <- years[y]
      
      for(fw_age in 1:N_fw_ages){
        for(mar_age in 1:N_mar_ages){
          
          idx <- (fw_age - 1) * N_mar_ages + mar_age
          
          if(y - mar_age >= 1){
            tTEMPcs <- TEMPcs$value[which(TEMPcs$year == (tyear - mar_age) & TEMPcs$r == reg)] # ocean entry
            tMLDcs <- MLDcs$value[which(MLDcs$year == (tyear - mar_age) & MLDcs$r == reg)] # ocean entry
          } else {
            tTEMPcs <- NA
            tMLDcs <- NA
          }
          
          if(y - (mar_age - 1) >= 1) {
            tTEMPwoo1 <- TEMPwoo$value[which(TEMPwoo$year == (tyear - (mar_age - 1)))] # winter1 
            tTEMPsoo1 <- TEMPsoo$value[which(TEMPsoo$year == (tyear - (mar_age - 1)))] # summer1
          } else {
            tTEMPwoo1 <- NA
            tTEMPsoo1 <- NA
          }
          
          if(y - (mar_age - 2) >= 1) {
            tTEMPwoo2 <- TEMPwoo$value[which(TEMPwoo$year == (tyear - (mar_age - 2)))] # winter2
            tTEMPsoo2 <- TEMPsoo$value[which(TEMPsoo$year == (tyear - (mar_age - 2)))] # summer2
          } else {
            tTEMPwoo2 <- NA
            tTEMPsoo2 <- NA
          }
          
          if(y - (mar_age - 3) >= 1) {
            tTEMPwoo3 <- TEMPwoo$value[which(TEMPwoo$year == (tyear - (mar_age - 3)))] # winter3
          } else {
            tTEMPwoo3 <- NA
          }
          
          tFTum <- FTum$value[which(FTum$year == (tyear) & FTum$p == p)] # upstream migration
          tFDum <- FDum$value[which(FDum$year == (tyear) & FDum$p == p)]
          
          # a_returns setup
          log_mar <- post$a_returns_logit[,idx,p] + post$a_returns_logit_year_sigma[,p] * 
            (if(!is.null(ann.dev)) rnorm(N_iter, mean = ta_ret_z_mean, sd = ta_ret_z_sd) else 0) 
          log_mar <- log_mar +
            post$b_TEMPcs[,p] * tTEMPcs + # needs to be the stnd data from Jan
            post$b_MLDcs[,p] * tMLDcs + # needs to be the stnd data from Jan
            post$b_TEMPwoo[,p,1] * tTEMPwoo1 +
            (if(mar_age >=2) post$b_TEMPsoo[,p,1] * tTEMPsoo1 else 0) +
            (if(mar_age >=2) post$b_TEMPwoo[,p,2] * tTEMPwoo2 else 0) +
            (if(mar_age >=3) post$b_TEMPsoo[,p,2] * tTEMPsoo2 else 0) +
            (if(mar_age >=3) post$b_TEMPwoo[,p,3] * tTEMPwoo3 else 0) +
            post$b_FTum[,p] * tFTum + 
            post$b_FDum[,p] * tFDum
          
          a_returns[,y,idx,p] <- inv.logit(log_mar)
          #b_returns[,y,idx,p] <- post$b_returns[,p]
        }
      }
    }
    
    # --- Population BH model --- #
    for(y in 1:length(years)){
      tyear <- years[y]
      
      # brood year (for indexing smolts)
      by <- y - 2
      byear <- tyear - 2
      
      # spawners - will need to change if we want a closed loop; spawners will need to be grabbed from returns estimation
      # Sji <- spawner_abd[by,p] # brood year spawners if we're doing closed loop
      if(by >= 1){
        Sji <- spawner_abd[1,p]
        spawners[,y,p] <- Sji
      } else {
        Sji <- spawner_abd[1,p]
        spawners[,y,p] <- Sji
      }
      
      ### spawner to smolts
      fry[,y,p] <- (a_fry[,y,p] * Sji) / (1 + post$b_fry[,p] * (Sji / 1e6))
      
      smolts_age[,y,1,p] <- fry[,y,p] * (1 - holdover[,y,p])
      holdover_fry[,y,p] <- fry[,y,p] * holdover[,y,p]
      
      if (y > 1){
        smolts_age[,y,2,p] <- holdover_fry[,y-1,p] * age2_survival[,y-1,p]
      } else {
        smolts_age[,y,2,p] <- holdover_fry[,1,p] * age2_survival[,1,p]
      }
      
      # presmolts for a given year = out migrating smolt_age1[by = y-2] + smolt_age2[by = y-3] + holdover_fry[by = y-2]
      presmolts[,y,p] <- smolts_age[,y,1,p] + holdover_fry[,y,p] + smolts_age[,y,2,p]
      
      # total smolts by smolt out-migration year
      smolts[,y,p] <- smolts_age[,y,1,p] + smolts_age[,y,2,p]
      
      ###  smolt to returns
      for(fw_age in 1:N_fw_ages){
        for(mar_age in 1:N_mar_ages){
          
          sy <- y - mar_age # define smolt migration year
          idx <- (fw_age - 1) * N_mar_ages + mar_age # define index for age class
          
          if(sy >= 1){
            tsmo <- smolts_age[,sy,fw_age,p]
          } else {
            tsmo <- smolts_age[,1,fw_age,p]
          }
          
          #returns_age[,y,idx,p] <- (a_returns[,y,idx,p] * tsmo) / (1 + post$b_returns[,p] * (tsmo / 1e6))
          returns_age[,y,idx,p] <- (a_returns[,y,idx,p] * tsmo)
          
          # calculate total returns by return year
          returns[,y,p] <- returns[,y,p] + returns_age[,y,idx,p];
          
          
        }
      }
    }
    
    for(y in 1:length(years)){
      # total smolt production by brood year
      smolts_by[,y,p] <- smolts_by[,y,p] +
        (if(y+2 <= length(years)) smolts_age[,y + 2,1,p] else smolts_age[,length(years), 1,p]) + 
        (if(y+3 <= length(years)) smolts_age[,y + 3,2,p] else smolts_age[,length(years), 2,p])
      
      for(fw_age in 1:N_fw_ages){
        for(mar_age in 1:N_mar_ages){
          
          idx <- (fw_age - 1) * N_mar_ages + mar_age # define index for age class
          
          # calculate returns by brood year
          if((y + (fw_age + 1 + mar_age)) <= length(years)){
            returns_by[,y,p] <- returns_by[,y,p] + returns_age[,y + (fw_age + 1 + mar_age), idx, p]
          } else {
            returns_by[,y,p] <- returns_by[,y,p] + returns_age[,length(years), idx, p]
          }
        }
      }
    }
    
  }
  
  
  # --- Store outputs --- #
  # combine indices
  out_smo <- as.data.frame.table(a_fry, responseName = "a_fry") |>
    full_join(as.data.frame.table(b_fry, responseName = "b_fry")) |>
    full_join(as.data.frame.table(age2_survival, responseName = "age2_survival")) |>
    full_join(as.data.frame.table(smolts_age, responseName = "smo_abd")) |>
    left_join(cuid_dat |> select(pop_name_short, domain), by = join_by("pop" == "pop_name_short")) |>
    mutate(year = as.numeric(as.character(year)),
           age_fw = as.numeric(as.character(age_fw)),
           brood_year = year - age_fw) |> # add brood year
    left_join(as.data.frame.table(spawners, responseName = "spa_abd") |> mutate(year = as.numeric(as.character(year))), 
              by = join_by("iter"=="iter", "pop"=="pop", "brood_year" == "year"))
    
  out_ret <- as.data.frame.table(a_returns, responseName = "a_returns") |>
    #full_join(as.data.frame.table(b_returns, responseName = "b_returns")) |>
    full_join(as.data.frame.table(returns_age, responseName = "ret_abd")) |>
    left_join(cuid_dat |> select(pop_name_short, domain), by = join_by("pop" == "pop_name_short")) |>
    mutate(year = as.numeric(as.character(year)),
           age_fw = as.numeric(substr(as.character(age),1,1)) + 1,
           age_mar = as.numeric(substr(as.character(age),3,3)),
           brood_year = year - (age_fw + age_mar)) |> # add brood year
    left_join(as.data.frame.table(spawners, responseName = "spa_abd") |> mutate(year = as.numeric(as.character(year))), 
              by = join_by("iter"=="iter", "pop"=="pop", "brood_year" == "year"))
  
  out_tot <- as.data.frame.table(spawners, responseName = "spa_abd") |>
    full_join(as.data.frame.table(smolts, responseName = "smo_abd")) |>
    full_join(as.data.frame.table(returns, responseName = "ret_abd")) |>
    left_join(cuid_dat |> select(pop_name_short, domain), by = join_by("pop" == "pop_name_short"))
  
  out_tot_by <- as.data.frame.table(spawners, responseName = "spa_abd") |>
    full_join(as.data.frame.table(smolts_by, responseName = "smo_abd")) |>
    full_join(as.data.frame.table(returns_by, responseName = "ret_abd")) |>
    left_join(cuid_dat |> select(pop_name_short, domain), by = join_by("pop" == "pop_name_short"))
  
  
  out_list <- list(out_smo = out_smo, out_ret = out_ret, out_tot = out_tot, out_tot_by = out_tot_by)
  return(out_list)
  
  ##### END NEW MODEL
}



##### READ DATA ##### 

year_start <- 1981

# historical data
env_dat <- read.csv("./data/env_data.csv") #PCIC observation driven data

## Current Historical time window baseline ~1981-2024
# Should use base case of 1991-2020 - present

# Jans code to process historical hindcast data

cuname_select <- c("Chilko Lake", "Great Central Lake","Sproat Lake","Tahltan Lake","Tatsamenie Lake","Quesnel Lake","Osoyoos Lake","Babine Lake (spawning channel)","Fraser Lake","Chilliwack Lake","Shuswap Lake","Francois Lake","Wenatchee Lake")

cu_reg <- data.frame(cu_name=c("Fraser Lake","Francois Lake","Chilko Lake","Chilliwack Lake","Cultus Lake","Shuswap Lake","Quesnel Lake","Osoyoos Lake","Wenatchee Lake","Great Central Lake","Sproat Lake","Babine Lake","Babine Lake (spawning channel)","Tahltan Lake","Tahltan Lake (enhanced)","Tatsamenie Lake","Tatsamenie Lake (enhanced)","Chutine Lake"), major_watershed=c("Fraser","Fraser","Fraser","Fraser","Fraser","Fraser","Fraser","Columbia","Columbia","Somass","Somass","Skeena","Skeena","Stikine","Stikine","Taku","Taku","Taku"))

#list CU id, set  to 3 or 4 domains 
cuid_dat <- cu_reg |> 
  filter(cu_name %in% cuname_select) |> 
  mutate(domain = case_when(major_watershed %in% c("Columbia", "Somass") ~ "South",
                            major_watershed == "Fraser" ~ "Fraser",
                            major_watershed  %in% c("Stikine", "Taku", "Skeena") ~ "North")) |>
  mutate(domain = factor(domain, levels = c("South", "Fraser", "North"))) |>
  arrange(domain, major_watershed, cu_name) |>
  mutate(cu_id = 1:n()) |> 
  mutate(cu_name = factor(cu_name, levels = cu_name)) |> 
  mutate(major_watershed = factor(major_watershed, levels = unique(major_watershed))) |> 
  mutate(domain_id = as.numeric(domain)) |> 
  mutate(major_watershed_id = as.numeric(major_watershed)) |> 
  mutate(pop_name_short = str_trim(str_remove(cu_name, "\\b[Ll]ake\\b"))) |> 
  mutate(pop_name_short = ifelse(pop_name_short == "Babine  (spawning channel)", "Babine", pop_name_short))

# Renaming columns for plotting and colouring
#in levelx reduce Babine Lake (spawning channel) to just Babine Lake
levelx=c("Osoyoos Lake","Wenatchee Lake","Great Central Lake","Sproat Lake","Chilko Lake","Chilliwack Lake","Francois Lake","Fraser Lake","Shuswap Lake","Quesnel Lake","Babine Lake","Tahltan Lake","Tatsamenie Lake")
levelxs=c("Osoyoos","Wenatchee","Great Central","Sproat","Chilko","Chilliwack","Francois","Fraser","Shuswap","Quesnel","Babine","Tahltan","Tatsamenie")
mcol = c("Osoyoos Lake"="hotpink","Wenatchee Lake"="pink2","Great Central Lake"="springgreen4","Sproat Lake"="springgreen2","Chilko Lake"="skyblue4","Chilliwack Lake"="skyblue1","Francois Lake"="blue1","Fraser Lake"="blue4","Shuswap Lake"="cyan4","Quesnel Lake"="cyan2","Babine Lake (spawning channel)"="gold1","Tahltan Lake"="purple4","Tatsamenie Lake"="orange3")
mcolx = c("Osoyoos Lake"="hotpink","Wenatchee Lake"="pink2","Great Central Lake"="springgreen4","Sproat Lake"="springgreen2","Chilko Lake"="skyblue4","Chilliwack Lake"="skyblue1","Francois Lake"="blue1","Fraser Lake"="blue4","Shuswap Lake"="cyan4","Quesnel Lake"="cyan2","Babine Lake"="gold1","Tahltan Lake"="purple4","Tatsamenie Lake"="orange3")
cus <- c("Osoyoos Lake","Wenatchee Lake","Great Central Lake","Sproat Lake","Chilko Lake","Chilliwack Lake","Francois Lake","Fraser Lake","Shuswap Lake","Quesnel Lake","Babine Lake","Tahltan Lake","Tatsamenie Lake")
cu_reg <- data.frame(cu_name=c("Fraser Lake","Francois Lake","Chilko Lake","Chilliwack Lake","Cultus Lake","Shuswap Lake","Quesnel Lake","Osoyoos Lake","Wenatchee Lake","Great Central Lake","Sproat Lake","Babine Lake","Babine Lake (spawning channel)","Tahltan Lake","Tahltan Lake (enhanced)","Tatsamenie Lake","Tatsamenie Lake (enhanced)","Chutine Lake"), major_watershed=c("Fraser","Fraser","Fraser","Fraser","Fraser","Fraser","Fraser","Columbia","Columbia","Somass","Somass","Skeena","Skeena","Stikine","Stikine","Taku","Taku","Taku"))
cuname_select <- c("Chilko Lake", "Great Central Lake","Sproat Lake","Tahltan Lake","Tatsamenie Lake","Quesnel Lake","Osoyoos Lake","Babine Lake (spawning channel)","Fraser Lake","Chilliwack Lake","Shuswap Lake","Francois Lake","Wenatchee Lake")
domainx<-c("South","Fraser","North")
names(domainx)<-c(1,2,3)

cuid_datx <- cuid_dat
cuid_datx$cu_name <- as.character(cuid_datx$cu_name)
cuid_datx$cu_name[cuid_datx$cu_name=="Babine Lake (spawning channel)"] = "Babine Lake"
cuid_datx$cu_name <- factor(cuid_datx$cu_name, levels=levelx)

major_watershed<-unique(cuid_datx$major_watershed)


#relabel "Babine Lake" to "... (spawning channel)" to match smo_dat naming
env_dat$cu[env_dat$cu=="Babine Lake"] = "Babine Lake (spawning channel)"

env_select <- env_dat |>
  filter(cu %in% cuname_select) |> 
  select(year, cu_name = cu, 
         FTsre = FTsreMean, 
         FTwre = FTwreMean,
         TEMPsoo = TEMPsooMean,
         TEMPwoo = TEMPwooMean,
         SSTsoo = SSTsooMean,
         SSTwoo = SSTwooMean,
         FTum = FTumMean, 
         FDum = FDumMean) |> 
  mutate(cu_name = factor(cu_name, levels = cuid_dat$cu_name)) |> 
  left_join(cuid_dat |> select(cu_name, major_watershed))

#Replace Wenatchee lake rearing data with Osoyoos data
for(r in 1:nrow(env_select)){if(env_select[r,2]=="Wenatchee Lake"){env_select[r,3]=env_select[r+1,3];env_select[r,4]=env_select[r+1,4]}}

env_long <- env_select |> 
  gather(key = variable, value = value, -cu_name, -year, -major_watershed) |> 
  group_by(variable) |> 
  mutate(value_stnd = (value - mean(value, na.rm = TRUE))/sd(value, na.rm = TRUE)) |> #Global stdn
  group_by(variable, cu_name, major_watershed) |> 
  mutate(value_stnd_cu = (value - mean(value, na.rm = TRUE))/sd(value, na.rm = TRUE)) |> #CU stdn
  left_join(cuid_dat |> select(cu_name, cu_id)) |> 
  mutate(year_id = year - year_start + 1) |> 
  ungroup()


coastal_dat <- read.csv("./data/env_data_region.csv")
coastal_dat$region <- factor(coastal_dat$region, levels=c("Columbia","Fraser","Somass","QCS","Hecate","Skeena","Stikine","Taku","Kodiak","OpenOcean"))

coastal_long <- coastal_dat |>
  filter(!region=="OpenOcean")|>
  select(-X)|>
  gather(key = variable, value = value, -region, -year, -area) |> 
  #group_by(variable, region) |> #regional stdn
  group_by(variable) |> #global stnd
  mutate(value_stnd = (value - mean(value, na.rm = TRUE))/sd(value, na.rm = TRUE)) |> 
  mutate(year_id = year - year_start + 1) |> 
  ungroup()|>
  mutate(region_id=as.numeric(region))|>
  select(-value) |> spread(key = variable, value = value_stnd)


# --- projection environmental data --- #

proj_dir <- "./"
files2 <- list.files(proj_dir, pattern = "projected245")

# filename <- files2[i]
# var_name <- strsplit(filename, "_")[[1]][1]
# tdat <- read.csv(paste0(proj_dir, filename))
# head(tdat)


# --- Imputed hindcast environmental data --- #

# read model output

###Updated model and model fit without b_returns
fit <- readRDS(file="./data/sockeye_climate_model_full.rds")

post <- extract_draws(fit)

# Jan's code to process imputed data for model output from Stan code 
# TEMPcs and MLDcs you would want to use the TEMPstnd and MLDstnd which have be normalized a second time in stan. The TEMPcs_complete and MLDcs_complete still are imputed and then integrated over the migration route, but are not used in the model.

FTwre <- spread_draws(fit, FTwre_complete[y,p]) %>% 
  group_by(y,p) %>% 
  summarize(value = mean(FTwre_complete, na.rm=TRUE)) %>% 
  mutate(year = y + (year_start-1), 
         pop = factor(cuid_datx$pop_name_short[p], levels = levelxs)) %>% 
  ungroup() #%>%
  # select(pop, year, value)

FTsre <- spread_draws(fit, FTsre_complete[y,p]) %>% 
  group_by(y,p) %>% 
  summarize(value = mean(FTsre_complete, na.rm=TRUE)) %>% 
  mutate(year = y + (year_start-1), 
         pop = factor(cuid_datx$pop_name_short[p], levels = levelxs)) %>% 
  ungroup() #%>%
  # select(pop, year, value)

FTum <- spread_draws(fit, FTum_complete[y,p]) %>% 
  group_by(y,p) %>% 
  summarize(value = mean(FTum_complete, na.rm=TRUE)) %>% 
  mutate(year = y + (year_start-1), 
         pop = factor(cuid_datx$pop_name_short[p], levels=levelxs)) %>% 
  ungroup() #%>%
  # select(pop, year, value)

FDum <- spread_draws(fit, FDum_complete[y,p]) %>% 
  group_by(y,p) %>% 
  summarize(value = mean(FDum_complete, na.rm=TRUE)) %>% 
  mutate(year = y + (year_start-1),
         pop = factor(cuid_datx$pop_name_short[p], levels=levelxs)) %>% 
  ungroup() #%>% 
  # select(pop, year, value)

# TEMPcs <- spread_draws(fit, TEMPcs_complete[y,r]) %>% 
#   group_by(y,r) %>% 
#   summarize(value = mean(TEMPcs_complete, na.rm=TRUE)) %>% 
#   mutate(year = y + (year_start-1), 
#          watershed = factor(major_watershed[r], levels = c("Columbia","Somass","Fraser","Skeena","Stikine","Taku"))) %>% 
#   ungroup() %>% 
#   select(watershed, year, value)

# MLDcs <- spread_draws(fit, MLDcs_complete[y,r]) %>% 
#   group_by(y,r) %>% 
#   summarize(value = mean(MLDcs_complete, na.rm=TRUE)) %>% 
#   mutate(year = y + (year_start-1),
#          watershed = factor(major_watershed[r], levels=c("Columbia","Somass","Fraser","Skeena","Stikine","Taku"))) %>% 
#   ungroup() %>% 
#   select(watershed, year, value)

TEMPcs <- spread_draws(fit, TEMPstnd[y,r]) %>% 
  group_by(y,r) %>% 
  summarize(value = mean(TEMPstnd,na.rm=TRUE)) %>% 
  mutate(year = y + (year_start-1),
         watershed = factor(major_watershed[r], levels=c("Columbia","Somass","Fraser","Skeena","Stikine","Taku"))) %>% 
  ungroup() #%>% 
  # select(watershed, year, value)

MLDcs <- spread_draws(fit, MLDstnd[y,r]) %>% 
  group_by(y,r) %>% 
  summarize(value = mean(MLDstnd, na.rm=TRUE)) %>% 
  mutate(year = y + (year_start-1), 
         watershed = factor(major_watershed[r], levels=c("Columbia","Somass","Fraser","Skeena","Stikine","Taku"))) %>% 
  ungroup() #%>% 
  # select(watershed, year, value)

TEMPwoo <- spread_draws(fit, TEMPwoo_complete[y]) %>% 
  group_by(y) %>% 
  summarize(value = mean(TEMPwoo_complete, na.rm=TRUE)) %>% 
  mutate(year = y + (year_start-1), 
         openocean = "Open Ocean") %>% 
  ungroup() #%>% 
  # select(openocean, year, value)

TEMPsoo <- spread_draws(fit, TEMPsoo_complete[y]) %>% 
  group_by(y) %>% 
  summarize(value = mean(TEMPsoo_complete, na.rm=TRUE)) %>% 
  mutate(year = y + (year_start-1), 
         openocean = "Open Ocean") %>% 
  ungroup() #%>% 
  # select(openocean, year, value)


# --- Salmon data --- #
se_adult <- read.csv("./data/adult_table_public.csv", header = TRUE)
se_juv <- read.csv("./data/juvenile_table_public.csv", header = TRUE)

#escapement #check: brood_year = adult_return_year ???
esc_dat <- se_adult %>%
  select(cu_name, adult_return_year, total_escape) %>%
  rename(brood_year = adult_return_year) %>%
  # filter(!cu_name %in% cuname_enhanced) %>%
  filter(brood_year >= year_start & cu_name%in%cuname_select)

##combine data
model_dat <- esc_dat %>%
  filter(brood_year >= year_start) %>% 
  left_join(cuid_dat, by = join_by(cu_name==cu_name))


#cu domains
north <- c("Skeena","Stikine","Taku")
south <- c("Somass","Columbia")
# cuid_dat$domain <- ifelse(cuid_dat$region %in% north, "North",ifelse(cuid_dat$region%in%south,"South","Fraser"))
cuid_dat$domainx <- as.numeric(as.factor(cuid_dat$domain))
doname<-unique(cuid_dat$domain)
names(doname)<-unique(cuid_dat$domainx)

# escapement data
e <- model_dat %>%
  select(cu_id, brood_year, total_escape) %>%
  distinct() %>%
  na.omit()

#means of log(escapes) for priors 1991-2020 period
e_med <- e %>% 
  filter(brood_year %in% c(1991:2020)) %>%
  group_by(cu_id) %>%
  summarise(e_median = median(total_escape, na.rm = TRUE)) %>%
  pivot_wider(names_from = "cu_id", values_from = "e_median") %>%
  ungroup() %>%
  as.matrix()

##### --- Settings --- #####

# settings

#post <- extract_draws(fit)
# S      <- 800       # n_draws
S <- nrow(post$smo_init)
N_pop  <- length(unique(model_dat$cu_id)) # n populaitons

# escapement data
spawner_abd <- e_med  # [N_years, N_pop]

# ann deviation reference
ann.dev.years <- NULL
ann.dev.years <- 1991:2020

##### --- Running model --- #####

# --- model run 1 - base period --- #

# years
years <- 1991:2020

# MAKE SURE TO RUN HINDCAST IMPUTED ENVIRONMENTAL DATA SECTION!

out_base <- proj1(N_pop = N_pop, years = years, N_fw_ages = 2, N_mar_ages = 3, N_iter = S, spawner_abd = spawner_abd, post = post, ann.dev = ann.dev.years)


# --- model run 2 - projection SSP4.5 --- #

# years
years <- 2041:2070

# set up projection data

# proj_dir <- "./"
# list.files(proj_dir, pattern = "projected245")

FTwre <- read.csv(paste0(proj_dir,"FTwre_projected245_2041_2070.csv")) %>% select(-X) %>%
  pivot_longer(cols = c(-year), names_to = "pop", values_to = "value") %>%
  mutate(pop = gsub(" Lake", "", gsub("\\.", " ", pop))) %>%
  left_join(cuid_dat %>% select(cu_id, pop_name_short), by = join_by(pop == pop_name_short)) %>%
  rename(p = cu_id)

FTsre <- read.csv(paste0(proj_dir,"FTsre_projected245_2041_2070.csv")) %>% select(-X) %>%
  pivot_longer(cols = c(-year), names_to = "pop", values_to = "value") %>%
  mutate(pop = gsub(" Lake", "", gsub("\\.", " ", pop))) %>%
  left_join(cuid_dat %>% select(cu_id, pop_name_short), by = join_by(pop == pop_name_short)) %>%
  rename(p = cu_id)

TEMPcs <- read.csv(paste0(proj_dir,"TEMPcs_projected245_2041_2070.csv")) %>% select(-X) %>%
  pivot_longer(cols = c(-year), names_to = "reg", values_to = "value") %>%
  left_join(cuid_dat %>% select(major_watershed, major_watershed_id) %>% distinct(), by = join_by(reg == major_watershed)) %>%
  rename(r = major_watershed_id)

MLDcs <- read.csv(paste0(proj_dir,"MLDcs_projected245_2041_2070.csv")) %>% select(-X) %>%
  pivot_longer(cols = c(-year), names_to = "reg", values_to = "value") %>%
  left_join(cuid_dat %>% select(major_watershed, major_watershed_id) %>% distinct(), by = join_by(reg == major_watershed)) %>%
  rename(r = major_watershed_id)

TEMPwoo <- read.csv(paste0(proj_dir,"TEMP_OpenOcean_projected245_2041_2070.csv")) %>% select(year, TEMPwoo) %>%
  rename(value = TEMPwoo)

TEMPsoo <- read.csv(paste0(proj_dir,"TEMP_OpenOcean_projected245_2041_2070.csv")) %>% select(year, TEMPsoo) %>%
  rename(value = TEMPsoo)

FDum <- read.csv(paste0(proj_dir,"FDum_projected245_2041_2070.csv")) %>% select(-X) %>%
  pivot_longer(cols = c(-year), names_to = "pop", values_to = "value") %>%
  mutate(pop = gsub(" Lake", "", gsub("\\.", " ", pop)),
         value = replace_na(value, 0)) %>%
  left_join(cuid_dat %>% select(cu_id, pop_name_short), by = join_by(pop == pop_name_short)) %>%
  rename(p = cu_id)

FTum <- read.csv(paste0(proj_dir,"FTum_projected245_2041_2070.csv")) %>% select(-X) %>%
  pivot_longer(cols = c(-year), names_to = "pop", values_to = "value") %>%
  mutate(pop = gsub(" Lake", "", gsub("\\.", " ", pop)),
         value = replace_na(value, 0)) %>%
  left_join(cuid_dat %>% select(cu_id, pop_name_short), by = join_by(pop == pop_name_short)) %>%
  rename(p = cu_id)

# run model
out45 <- proj1(N_pop = N_pop, years = years, N_fw_ages = 2, N_mar_ages = 3, N_iter = S, spawner_abd = spawner_abd, post = post, ann.dev = ann.dev.years)


##### OUTPUT PLOTTING #####

# --- Data wrangling --- #

# Find means of values across iterations

tdatb <- out_base$out_tot |> mutate(scen = "base") |> 
  mutate(brood_year = as.numeric(as.character(year))) |>
  left_join(cuid_dat |> select(cu_name, pop_name_short), by = join_by(pop == pop_name_short)) |>
  group_by(pop, domain, scen, iter) |>
  summarize(spa_abd = mean(spa_abd, na.rm=TRUE),
            smo_abd = mean(smo_abd, na.rm=TRUE),
            ret_abd = mean(ret_abd, na.rm=TRUE)) |>
  ungroup() 
tdat45 <- out45$out_tot |> mutate(scen = "ssp45") |>   
  mutate(brood_year = as.numeric(as.character(year))) |>
  left_join(cuid_dat |> select(cu_name, pop_name_short), by = join_by(pop == pop_name_short)) |>
  group_by(pop, domain, scen, iter) |>
  summarize(spa_abd = mean(spa_abd, na.rm=TRUE),
            smo_abd = mean(smo_abd, na.rm=TRUE),
            ret_abd = mean(ret_abd, na.rm=TRUE)) |>
  ungroup()
# tdat85 <- out85$out_tot |> mutate(scen = "ssp85") |> 
#   mutate(brood_year = as.numeric(as.character(year))) |>
#   left_join(cuid_dat |> select(cu_name, pop_name_short), by = join_by(pop == pop_name_short)) |>
#   group_by(pop, domain, scen, iter) |>
#   summarize(spa_abd = mean(spa_abd, na.rm=TRUE),
#             smo_abd = mean(smo_abd, na.rm=TRUE),
#             ret_abd = mean(ret_abd, na.rm=TRUE)) |>
#   ungroup()

head(tdatb)
head(tdat45)


#####
# ssp45 v base  

tdat1 <- tdatb |>
  mutate(ret1 = ret_abd)
tdat2 <- tdat45 |>
  mutate(ret2 = ret_abd)

head(tdat1)
head(tdat2)

plot.dat45 <- tdat1 |>
  bind_cols(ret2 = tdat2$ret2) |>
  mutate(ret_diff = ret2 - ret1,
         ret_pch = ret_diff/ret1*100,
         rs1 = ret1/spa_abd,
         rs2 = ret2/spa_abd,
         rs_per = rs2 / rs1 * 100) |>
  left_join(cuid_dat, by = join_by("domain" == "domain", 
                                   "pop" == "pop_name_short")) 

head(plot.dat45)

#Projections-Figure####
plot.dat45$popx <- factor(plot.dat45$pop,levels=levelxs)
axis_range <- c(0.1,3000)
minorb <- c(seq(1,10,1),seq(20,90,10),seq(100,900,100),seq(1000,3000,1000))
majorb <- c(5,10,50,100,500,1000)

ggplot(plot.dat45, aes(x = popx, domain_id, y = rs_per,fill = domain)) +
  scale_fill_manual(values=c('darkorchid3','goldenrod3','seagreen3'))+
  stat_halfeye(color = "grey",size=2,normalize='xy')+
  #theme_bw(base_size=12)+
  theme_bw(base_size=14)+
  labs(title=NULL,fill=NULL,y="Productivity %",x=NULL)+
  theme(legend.position='bottom')+
  #scale_fill_discrete(guide = "none")+
  geom_hline(yintercept=100, linetype= 2)+
  coord_flip()+
  scale_y_log10(limits=axis_range, minor_breaks=minorb, breaks=majorb,expand=c(0,0))

