###This script prepares the siumlated data to run the adapted stan model "sockeye_climate_model_sim.stan" to assess model performance and parameter recovery###

library(tidyverse)
library(janitor)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(bayesplot)
library(sockeyeSim)

#functions####
env_to_matrix <- function(df, name, spread_name = "cu_name"){
  data_out <- df |> 
    select(all_of(c(spread_name, "year", name)))|> 
    arrange(year, spread_name) |> 
    spread(key = spread_name, value = name) |> 
    ungroup() |> 
    select(-year) |> 
    data.matrix()
  data_out[is.na(data_out)] <- -9999
  return(data_out)
}

missID_matrix <- function(df, name, spread_name = "cu_name") {
  data_out <- env_to_matrix(df = df, name = name, spread_name = spread_name)
  miss_id <- apply(data_out, 2, function(x) as.integer(x == -9999))
  return(miss_id)
}

pop_means <- function(df, name, spread_name = "cu_name"){
  data_out <- df |> 
    select(all_of(c(spread_name, "year", name)))|> 
    arrange(year, spread_name) |> 
    spread(key = spread_name, value = name) |> 
    ungroup() |> 
    select(-year) |> 
    data.matrix()|>
    colMeans(,na.rm=TRUE)
  return(data_out)
}

pop_sd <- function(df, name, spread_name = "cu_name"){
  data_out <- df |> 
    select(all_of(c(spread_name, "year", name)))|> 
    arrange(year, spread_name) |> 
    spread(key = spread_name, value = name) |> 
    ungroup() |> 
    select(-year) |> 
    summarise(across(everything(),sd, na.rm = TRUE))|>
    unlist()
  return(data_out)
}


cuname_select <- c("Chilko Lake", "Great Central Lake","Sproat Lake","Tahltan Lake","Tatsamenie Lake","Quesnel Lake","Osoyoos Lake","Babine Lake (spawning channel)","Fraser Lake","Chilliwack Lake","Shuswap Lake","Francois Lake","Wenatchee Lake")

cu_reg <- data.frame(cu_name=c("Fraser Lake","Francois Lake","Chilko Lake","Chilliwack Lake","Cultus Lake","Shuswap Lake","Quesnel Lake","Osoyoos Lake","Wenatchee Lake","Great Central Lake","Sproat Lake","Babine Lake","Babine Lake (spawning channel)","Tahltan Lake","Tahltan Lake (enhanced)","Tatsamenie Lake","Tatsamenie Lake (enhanced)","Chutine Lake"), major_watershed=c("Fraser","Fraser","Fraser","Fraser","Fraser","Fraser","Fraser","Columbia","Columbia","Somass","Somass","Skeena","Skeena","Stikine","Stikine","Taku","Taku","Taku"))

#list CU id, set  to 3 or 4 domains 
cuid_dat <- cu_reg |> 
  filter(cu_name %in% cuname_select) |> 
  # mutate(domain = case_when(major_watershed %in% c("Columbia") ~ "CRB",
  #                           major_watershed %in% c("Somass") ~ "WVI",
  #                           major_watershed == "Fraser" ~ "Fraser",
  #                           major_watershed  %in% c("Stikine", "Taku", "Skeena") ~ "North")) |>
  # mutate(domain = factor(domain, levels = c("CRB","WVI", "Fraser", "North"))) |>
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

year_start <- 1981

#Simulation data####
#sim_results <- sockeye_sim()
sim_results <- readRDS("./output/sim_results_20260630.rds")

sim.df <- sim_results$sim_data.df
sim.age.returns <- sim_results$age_prop_returns
sim.age.smolts <- sim_results$age_prop_smolts
#sim.parameters <- sim_results$env_slopes
#sim.mu <- sim_results$env_slopes_mu
sim.env <- sim_results$env.df
sim.coastal <- sim_results$coastal_env.df

#Salmon data####
all_totals <- sim.df|>
  mutate(cu_id = pop_id, cu_name = cuid_dat$cu_name[cu_id])|>
  mutate(year_id = year - year_start +1)|>
  filter(cu_name %in% cuname_select) |> 
  mutate(cu_name = factor(cu_name, levels = cuid_dat$cu_name))
  #mutate(presmolts = 0)

smolt_age_prop <- sim.age.smolts|>
  select(pop_id,year,fw,age_prop_obs,age_n_eff)|>
  spread(key = fw, value = age_prop_obs, fill = 0)|>
  mutate(cu_id = pop_id, year_id = year - year_start +1)|>
  select(cu_id,year_id,age_1_x = `1`, age_2_x = `2`,age_n_eff)
  
return_age_prop <- sim.age.returns|>
  mutate(age = paste("age", fw, mar, sep = "_")) |> 
  select(pop_id, year, age, age_prop_obs,age_n_eff) |> 
  spread(key = age, value = age_prop_obs, fill = 0) |> 
  mutate(cu_id = pop_id, year_id = year - year_start +1, n_eff = 1)|>
  select(cu_id,year_id,age_1_1:age_2_3,age_n_eff)

juv_summary <- all_totals |> left_join(smolt_age_prop) |> 
  #mutate(juveniles = ifelse(!is.na(smolts), smolts, presmolts)) |> 
  mutate(juveniles = smolts_obs) |>
  select(cu_name, year, juveniles, age_1_x, age_2_x) |> 
  gather(key = age, value = prop, age_1_x:age_2_x) |> 
  mutate(total = juveniles * prop) |> 
  mutate(total = total + 1) |> 
  group_by(cu_name, age) |> 
  summarise(mean_ln = mean(log(total), na.rm = TRUE), sd_ln = sd(log(total), na.rm = TRUE)) |> 
  ungroup()

mean_juv_ln <- juv_summary |> 
  select(-sd_ln) |> 
  spread(key = age, value = mean_ln)

sd_juv_ln <- juv_summary |> 
  select(-mean_ln) |> 
  mutate(sd_ln = ifelse(sd_ln ==0, 0.2, sd_ln)) |> 
  spread(key = age, value = sd_ln)

ret_summary <- all_totals |> left_join(return_age_prop) |> 
  select(cu_name, year, returns_obs, age_1_1:age_2_3) |> 
  gather(key = age, value = prop, age_1_1:age_2_3) |> 
  mutate(total = returns_obs * prop) |> 
  mutate(total = total + 1) |> 
  group_by(cu_name, age) |> 
  summarise(mean_ln = mean(log(total), na.rm = TRUE), sd_ln = sd(log(total), na.rm = TRUE)) |> 
  ungroup()

mean_ret_ln <- ret_summary |> 
  select(-sd_ln) |> 
  spread(key = age, value = mean_ln)


mean_ret_ln_array <- array(NA, dim = c(2, 3, nrow(cuid_dat)))
for(i in 1:nrow(cuid_dat)){
  mean_ret_ln_array[,,i] <- matrix(as.numeric(mean_ret_ln[i, -1]), nrow = 2, ncol = 3, byrow = TRUE)
}


sd_ret_ln <- ret_summary |> 
  select(-mean_ln) |> 
  mutate(sd_ln = ifelse(sd_ln ==0, 0.2, sd_ln)) |> 
  spread(key = age, value = sd_ln)

sd_ret_ln_array <- array(NA, dim = c(2, 3, nrow(cuid_dat)))
for(i in 1:nrow(cuid_dat)){
  sd_ret_ln_array[,,i] <- matrix(as.numeric(sd_ret_ln[i, -1]), nrow = 2, ncol = 3, byrow = TRUE)
}

#Env data####
env_select <- sim.env |>
  mutate(cu_name = cuid_dat$cu_name[pop_id])|>
  filter(cu_name %in% cuname_select)|> 
  select(year, cu_name, FTsre,FTwre,TEMPsoo,TEMPwoo,FTum,FDum)|> 
  mutate(cu_name = factor(cu_name, levels = cuid_dat$cu_name)) |> 
  #mutate(GCMTsre=FTsre,GCMTwre=FTwre,SSTsoo=TEMPsoo,SSTwoo=TEMPwoo,GCMTum=FTum,GCMDum=FDum)|>
  left_join(cuid_dat |> select(cu_name, major_watershed))|>
  mutate(year_id = year - year_start + 1)

#random NAs
# env_select2 <- env_select|>select(year:GCMDum,year_id)
# for(i in 1:100){
#   c=base::sample(3:8,1,replace=TRUE)
#   r=base::sample(1:nrow(env_select2),1,replace=TRUE)
#   env_select2[r,c] = NA}

rivers = unique(cuid_dat$major_watershed)
coastal_select <- sim.coastal |>
  mutate(river = rivers[river_id])|>
  select(year, river, TEMP=TEMPstnd,MLD=MLDstnd)|> 
  mutate(SST=TEMP,MLDbccm=MLD)|>
  mutate(year_id = year - year_start + 1)

# coastal_select2 <- coastal_select|>select(year:MLDbccm)
# for(i in 1:100){
#   c=base::sample(3:4,1,replace=TRUE)
#   r=base::sample(1:nrow(coastal_select2),1,replace=TRUE)
#   coastal_select2[r,c] = NA}

#stan data list####
model_data <- list(
  N_years = max(env_select$year_id),
  N_pop = length(cuid_dat$cu_name),
  N_rivers = max(cuid_dat$major_watershed_id),
  N_domains = max(cuid_dat$domain_id),
  N_spawner_obs = sum(!is.na(all_totals$spawners_obs)), 
  #N_presmolt_obs = sum(!is.na(all_totals$presmolts)),
  N_smolt_obs = sum(!is.na(all_totals$smolts_obs)),
  N_smolt_age_obs = nrow(smolt_age_prop),
  N_return_age_obs = nrow(return_age_prop),
  
  N_fw_ages = 2,
  N_mar_ages = 3,
  
  #pop IDs
  river_id = cuid_dat$major_watershed_id, 
  domain_id = cuid_dat$domain_id,
  
  #abundance observations
  #presmolt_obs = all_totals$presmolts[!is.na(all_totals$presmolts)],
  #presmolt_year = all_totals$year_id[!is.na(all_totals$presmolts)], 
  #presmolt_pop = all_totals$cu_id[!is.na(all_totals$presmolts)],

  smolt_obs = all_totals$smolts_obs[!is.na(all_totals$smolts_obs)],
  smolt_year = all_totals$year_id[!is.na(all_totals$smolts_obs)], 
  smolt_pop = all_totals$cu_id[!is.na(all_totals$smolts_obs)],
  
  spawner_obs = all_totals$spawners_obs[!is.na(all_totals$spawners_obs)],
  spawner_year = all_totals$year_id[!is.na(all_totals$spawners_obs)], 
  spawner_pop = all_totals$cu_id[!is.na(all_totals$spawners_obs)],
  
  #age class proportion observations
  smolt_age_prop = smolt_age_prop |> ungroup() |> select(age_1_x, age_2_x),
  smolt_age_n_eff = smolt_age_prop$age_n_eff,
  smolt_age_year = smolt_age_prop$year_id,
  smolt_age_pop = smolt_age_prop$cu_id,
  
  return_age_prop = return_age_prop |> ungroup() |> select(age_1_1:age_2_3),
  return_age_n_eff = return_age_prop$age_n_eff,
  #return_age_n_eff = rep(100, nrow(return_age_prop)),
  return_age_year = return_age_prop$year_id,
  return_age_pop = return_age_prop$cu_id,
  
  #initialization matrices
  s_priors_ln = mean_juv_ln |> select(-cu_name) |> data.matrix() |> t(), 
  s_priors_sd = sd_juv_ln |> select(-cu_name) |> data.matrix() |> t(),
  
  r_priors_ln = mean_ret_ln_array, 
  r_priors_sd = sd_ret_ln_array, 
  
  #Env data
  FTsre = env_to_matrix(env_select, "FTsre"),
  FTsre_missid = missID_matrix(env_select, "FTsre"),  
  N_FTsre_miss = sum(is.na(env_select$FTsre)),
  
  # GCMTsre = env_to_matrix(env_select, "GCMTsre"),
  # GCMTsre_missid = missID_matrix(env_select, "GCMTsre"),  
  # #N_GCMTsre_miss = sum(is.na(env_select$GCMTsre)),
  # GCMTsre_pop_mean = pop_means(env_select|>filter(year%in%2013:2024), "GCMTsre"),
  # GCMTsre_pop_sd = pop_sd(env_select|>filter(year%in%2013:2024), "GCMTsre"),
  
  FTwre = env_to_matrix(env_select, "FTwre"),
  FTwre_missid = missID_matrix(env_select, "FTwre"),  
  N_FTwre_miss = sum(is.na(env_select$FTwre)),
  
  # GCMTwre = env_to_matrix(env_select, "GCMTwre"),
  # GCMTwre_missid = missID_matrix(env_select, "GCMTwre"),  
  # #N_GCMTwre_miss = sum(is.na(env_select$GCMTwre)),
  # GCMTwre_pop_mean = pop_means(env_select|>filter(year%in%2013:2024), "GCMTwre"),
  # GCMTwre_pop_sd = pop_sd(env_select|>filter(year%in%2013:2024), "GCMTwre"),
  
  TEMP = env_to_matrix(coastal_select, "TEMP", spread_name = "river"),
  TEMP_missid = missID_matrix(coastal_select, "TEMP", spread_name = "river"),
  N_TEMP_miss = sum(is.na(coastal_select$TEMP)),
  
  # SST = env_to_matrix(coastal_select, "SST", spread_name = "river"),
  # SST_missid = missID_matrix(coastal_select, "SST", spread_name = "river"),
  # N_SST_miss = sum(is.na(coastal_select$SST)),
  
  MLD = env_to_matrix(coastal_select, "MLD", spread_name = "river"),
  MLD_missid = missID_matrix(coastal_select, "MLD", spread_name = "river"),
  N_MLD_miss = sum(is.na(coastal_select$MLD)),
  
  # MLDbccm = env_to_matrix(coastal_select, "MLDbccm", spread_name = "river"),
  # MLDbccm_missid = missID_matrix(coastal_select, "MLDbccm", spread_name = "river"),
  # N_MLDbccm_miss = sum(is.na(coastal_select$MLDbccm)),

  TEMPsoo = env_to_matrix(env_select, "TEMPsoo", spread_name = "cu_name")[,1],
  TEMPsoo_missid = missID_matrix(env_select, "TEMPsoo", spread_name = "cu_name")[,1],  
  N_TEMPsoo_miss = sum(is.na(env_select$TEMPsoo)),
  
  TEMPwoo = env_to_matrix(env_select, "TEMPwoo", spread_name = "cu_name")[,1],
  TEMPwoo_missid = missID_matrix(env_select, "TEMPwoo", spread_name = "cu_name")[,1],  
  N_TEMPwoo_miss = sum(is.na(env_select$TEMPwoo)),
  
  # SSTsoo = env_to_matrix(env_select, "SSTsoo", spread_name = "cu_name")[,1],
  # #SSTsoo_missid = missID_matrix(env_select, "SSTsoo", spread_name="region")[,1],  
  # #N_SSTsoo_miss = sum(is.na(env_select$SSTsoo)),
  
  # SSTwoo = env_to_matrix(env_select, "SSTwoo", spread_name = "cu_name")[,1],
  # #SSTwoo_missid = missID_matrix(env_select, "SSTwoo", spread_name="region")[,1],  
  # #N_SSTwoo_miss = sum(is.na(env_select$SSTwoo)),
  
  FTum = env_to_matrix(env_select, "FTum"),
  FTum_missid = missID_matrix(env_select, "FTum"),  
  N_FTum_miss = sum(is.na(env_select$FTum)),
  
  # GCMTum = env_to_matrix(env_select, "GCMTum"),
  # GCMTum_missid = missID_matrix(env_select, "GCMTum"),  
  # #N_GCMTum_miss = sum(is.na(env_select$GCMTum)),
  # GCMTum_pop_mean = pop_means(env_select|>filter(year%in%2013:2024), "GCMTum"),
  # GCMTum_pop_sd = pop_sd(env_select|>filter(year%in%2013:2024), "GCMTum"),
  
  FDum = env_to_matrix(env_select, "FDum"),
  FDum_missid = missID_matrix(env_select, "FDum"),  
  N_FDum_miss = sum(is.na(env_select$FDum)),
  
  # GCMDum = env_to_matrix(env_select, "GCMDum"),
  # GCMDum_missid = missID_matrix(env_select, "GCMDum"),  
  # #N_GCMDum_miss = sum(is.na(env_select$GCMDum)),
  # GCMDum_pop_mean = pop_means(env_select|>filter(year%in%2013:2024), "GCMDum"),
  # GCMDum_pop_sd = pop_sd(env_select|>filter(year%in%2013:2024), "GCMDum"),
  
  smolt_alpha1_prior_mean = 50,
  #smolt_alpha2_prior_mean = 1,
  #smolt_alpha1_year_sd_prior = 0.5,
  #smolt_alpha2_year_sd_prior = 0.5,
  
  return_alpha_ref_prior_mean = 0.015,
  
  smolt_beta_prior_mean = 7, 
  smolt_beta_prior_sd = 1,
  return_beta_prior_mean = 0.0001, #0.001 or 0.0001?
  return_beta_prior_sd = 1
)

#initialization function####
init_fun <- function() list(
  a_returns_logit_mu = rep(qlogis(model_data$return_alpha_ref_prior_mean),6), 
  a_fry_mu_global_mu = 4,
  a_fry_mu_global_sigma = 0.5,
  sigma_returns = c(4,1,2.5,4,2,3.5),
  b_fry_mu = 2,
  b_fry_sigma = 1,
  a_returns_logit_year_sigma = rep(0.5,model_data$N_pop),
  #b_returns_mu = log(model_data$return_beta_prior_mean),
  #b_returns_sigma = 0.5,
  log_spawner_kappa_mu = log(0.1), #0.05, 0.03
  log_spawner_kappa_sigma = 0.5,
  log_smolt_kappa_mu   = log(0.05),
  log_smolt_kappa_sigma = 0.5,
  b_FTwre_hyper_mu = 0,
  b_FTwre_hyper_sigma = 0.25,
  b_FTwre_sigma = matrix(0.5,nrow=model_data$N_domains,ncol=3),
  #b_FTwre_sigma = rep(0.5,3),
  b_FTsre_hyper_mu = 0,
  b_FTsre_hyper_sigma = 0.25,
  b_FTsre_sigma = matrix(0.5,nrow=model_data$N_domains,ncol=2),
  #b_FTsre_sigma = rep(0.5,2),
  b_TEMPcs_hyper_mu = 0,
  b_TEMPcs_hyper_sigma = 0.5,
  b_TEMPcs_sigma = rep(0.5,model_data$N_domains),
  b_MLDcs_hyper_mu = 0,
  b_MLDcs_hyper_sigma = 0.5,
  b_MLDcs_sigma = rep(0.5,model_data$N_domains),
  b_TEMPwoo_hyper_mu = 0,
  b_TEMPwoo_hyper_sigma = 0.5,
  b_TEMPwoo_sigma = matrix(0.5,nrow=model_data$N_domains,ncol=3),
  #b_TEMPwoo_sigma = rep(0.5,3),
  b_TEMPsoo_hyper_mu = 0,
  b_TEMPsoo_hyper_sigma = 0.5,
  b_TEMPsoo_sigma = matrix(0.5,nrow=model_data$N_domains,ncol=2),
  #b_TEMPsoo_sigma = rep(0.5,2),
  b_FTum_hyper_mu = 0,
  b_FTum_hyper_sigma = 0.5,
  b_FTum_sigma = rep(0.5,model_data$N_domains),
  b_FDum_hyper_mu = 0,
  b_FDum_hyper_sigma =0.5,
  b_FDum_sigma = rep(0.5,model_data$N_domains),
  #TEMP_SST_b = 1,
  #MLD_MLDbccm_b = rep(1,model_data$N_rivers),
  #MLD_TEMP_b = rep(-1,model_data$N_rivers),
  holdover_logit_mu = rep(qlogis(0.1),model_data$N_pop), #pop specific?
  holdover_logit_sigma = rep(0.5,model_data$N_pop),
  age2_survival_logit_mu = rep(qlogis(0.8),model_data$N_pop),
  age2_survival_logit_sigma = rep(0.5,model_data$N_pop)
)


#stan model####
model <- cmdstan_model("./sockeye_climate_model_sim.stan")

system.time(fit <- model$sample(data = model_data,  
                                init = function() init_fun(), 
                                chains = 4, 
                                parallel_chains = 4, 
                                iter_warmup = 1000,
                                save_warmup = TRUE,
                                iter_sampling = 2000,
                                thin = 2,
                                refresh = 50, 
                                max_treedepth = 11,
                                adapt_delta = 0.99,
                                seed = 9))

#assess convergence####
worst_Rhat <-fit$summary() %>% 
  as.data.frame() %>% 
  mutate(rhat  = round(rhat , 3)) %>% 
  arrange(desc(rhat)) %>% 
  filter(!is.na(rhat))
mcmc_trace(fit$draws(), pars = worst_Rhat$variable[1:20], np = nuts_params(fit))

#saveRDS(fit, "./output/sockeye_climate_model_MVN_env_sim3_20260630.rds")
