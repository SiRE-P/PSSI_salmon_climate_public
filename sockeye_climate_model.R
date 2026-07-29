library(tidyverse)
library(janitor)
library(cmdstanr)
library(tidybayes)
library(posterior)
library(bayesplot)

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

#load data####
se_adult <- read.csv("./data/adult_table.csv", header = TRUE)
se_juv <- read.csv("./data/juvenile_table.csv", header = TRUE)

cuname_enhanced <- c("Tahltan Lake (enhanced)", "Tatsamenie Lake (enhanced)")

all_cunames <- c("Chilko Lake", "Great Central Lake","Sproat Lake","Tahltan Lake","Tatsamenie Lake","Quesnel Lake","Osoyoos Lake","Babine Lake (spawning channel)","Francois Lake","Fraser Lake","Chilliwack Lake","Shuswap Lake","Wenatchee Lake")

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

#juvenile data####
juv_data <- se_juv |> 
  filter(cu_name %in% cuname_select) |> 
  mutate(cu_name = factor(cu_name, levels = cuid_dat$cu_name)) |> 
  select(cu_name, smolt_migration_year, smolt_abundance_age0, smolt_abundance_age1, smolt_abundance_age2, smolt_abundance_age3, smolt_abundance_total, presmolt_abundance) |> 
  filter(smolt_migration_year >= year_start) |> 
  group_by(cu_name, smolt_migration_year) |> 
  mutate(smolt_abundance_total = ifelse(is.na(smolt_abundance_total), sum(smolt_abundance_age0, smolt_abundance_age1, smolt_abundance_age2, smolt_abundance_age3, na.rm = TRUE), smolt_abundance_total)) |> 
  mutate(smolt_abundance_total = ifelse(smolt_abundance_total == 0, NA, smolt_abundance_total)) |> 
  mutate(smolt_abundance_total = ifelse(!is.na(smolt_abundance_total) &
                                          !is.na(presmolt_abundance) &
                                          smolt_abundance_total == presmolt_abundance,
                                        NA, smolt_abundance_total)) |>  
  filter(!is.na(smolt_abundance_total) | !is.na(presmolt_abundance)) |> 
  arrange(cu_name, smolt_migration_year) |> 
  rename(smolts = smolt_abundance_total, presmolts = presmolt_abundance) |> 
  mutate(smolts = smolts - rowSums(across(c(smolt_abundance_age0, smolt_abundance_age3)), na.rm = TRUE)) |> #remove age0s and 3s from totals
  select(-smolt_abundance_age0, -smolt_abundance_age3) |> 
  mutate(prop_age2 = smolt_abundance_age2 / (smolt_abundance_age1 + smolt_abundance_age2)) |> 
  mutate(smolt_abundance_age1 = ifelse(cu_name == "Great Central Lake" & prop_age2 == 0.1, NA, smolt_abundance_age1), #removing deterministic juv age classification
         smolt_abundance_age2 = ifelse(cu_name == "Great Central Lake" & prop_age2 == 0.1, NA, smolt_abundance_age2)) |> 
  mutate(smolt_abundance_age1 = ifelse(cu_name == "Sproat Lake" & (abs(prop_age2 - 0.03) < 1e-6 | abs(prop_age2 - 0.028) < 1e-6 | abs(prop_age2 - 0.03061224) < 1e-6), NA, smolt_abundance_age1),
         smolt_abundance_age2 = ifelse(cu_name == "Sproat Lake" & (abs(prop_age2 - 0.03) < 1e-6 | abs(prop_age2 - 0.028) < 1e-6 | abs(prop_age2 - 0.03061224) < 1e-6), NA, smolt_abundance_age2)) |> 
  select(-prop_age2)

juv_select_table <- data.frame(cu_name = all_cunames, type = c("smolts", "presmolts", "presmolts", "smolts", "smolts", rep("presmolts", 8)))

juv_total_abund <- juv_data |> 
  left_join(juv_select_table) |> 
  mutate(cu_name = factor(cu_name, levels = cuid_dat$cu_name)) |> 
  mutate(juv_abund = ifelse(type == "smolts", smolts, presmolts)) |> 
  select(cu_name, smolt_migration_year, juv_abund, type) |> 
  filter(!is.na(juv_abund))

juv_total_abund <- juv_total_abund |> 
  spread(key = type, value = juv_abund, fill = NA)

#adults####
adult_data <- se_adult%>%filter(!en_route_mort=="included") |> 
  filter(cu_name %in% cuname_select) |> 
  mutate(cu_name = factor(cu_name, levels = cuid_dat$cu_name)) |> 
  select(cu_name, adult_return_year, starts_with("total_returns"), -total_returns_total, -total_returns_other, -total_returns_comments) %>%
  pivot_longer(cols = c(-cu_name, -adult_return_year), names_to = "adult_age", values_to = "adult_abd") %>%
  mutate(smolt_age = as.numeric(substr(adult_age,15,15)),
         adult_age = as.numeric(substr(adult_age,17,17))) %>%
  filter(!cu_name %in% cuname_enhanced) %>%
  mutate(brood_year = adult_return_year-smolt_age-adult_age) %>%
  filter(smolt_age %in% c(1,2)) %>%
  filter(adult_age %in% c(1,2,3)) |>
  group_by(cu_name, adult_return_year) |> 
  arrange(cu_name, adult_return_year, smolt_age, adult_age) |> 
  filter(adult_return_year >= year_start) |> 
  group_by(cu_name, adult_return_year) |> 
  mutate(total_return = sum(adult_abd, na.rm = TRUE), n_ages = n()) |> 
  filter(total_return>0, n_ages == 6) |> 
  mutate(adult_abd = ifelse(is.na(adult_abd), 0, adult_abd)) |> 
  mutate(smolt_migration_year = adult_return_year - adult_age)

return_total <- adult_data |> 
  group_by(cu_name, adult_return_year) |> 
  summarise(adult_returns = sum(adult_abd, na.rm = TRUE)) |> 
  rename(year = adult_return_year)


#hatchery inputs and fishing####
fishing_hatchery <- se_adult |> 
  ungroup() |> 
  filter(!en_route_mort=="included") |> 
  mutate(hatchery_escape = ifelse(is.na(hatchery_escape), 0, hatchery_escape)) |> 
  mutate(hatchery_escape = ifelse(total_returns_total + hatchery_escape - total_escape < 0, -(total_returns_total + hatchery_escape - total_escape), hatchery_escape)) |> #check this one with Cam - Tahltan lake has negative fishing in some years, I assume it is because of hatchery, which is not separated in this pop
  mutate(fishing = total_returns_total + hatchery_escape - total_escape) |> 
  mutate(fishing_prop = fishing/total_returns_total, na.rm = TRUE) |> 
  #mutate(fishing_prop = ifelse(fishing_prop <0, 0, fishing_prop)) |> 
  select(adult_return_year, cu_name, hatchery_escape, total_escape, fishing, fishing_prop)

#compare smolt age between smolt and adult data####
juv_compare <- juv_data |> 
  select(cu_name, smolt_migration_year, smolt_abundance_age1, smolt_abundance_age2) |>
  #filter(cu_name %in% c("Chilko Lake", "Sproat Lake", "Great Central Lake", "Tahltan Lake", "Tatsamenie Lake")) |> 
  filter(cu_name %in% cuname_select) |> 
  mutate(prop_age2 = smolt_abundance_age2 / (smolt_abundance_age1 + smolt_abundance_age2)) |> 
  select(cu_name, smolt_migration_year, prop2_smoltdata = prop_age2) |> 
  left_join(adult_data |> 
              group_by(cu_name, smolt_migration_year) |> 
              mutate(smolt_year_total = sum(adult_abd), n_ages_smolt_year = n()) |> 
              filter(n_ages_smolt_year == 6) |> 
              group_by(cu_name, smolt_migration_year, smolt_age, smolt_year_total) |> 
              summarise(smolt_age_abd = sum(adult_abd)) |> 
              mutate(age_prop = smolt_age_abd/smolt_year_total) |> 
              filter(smolt_age == 2) |> 
              ungroup() |> 
              select(cu_name, smolt_migration_year, prop2_adultdata = age_prop))

#spawner totals####
spawner_total <- se_adult%>%filter(!en_route_mort=="included") |> 
  filter(cu_name %in% cuname_select) |> 
  mutate(cu_name = factor(cu_name, levels = cuid_dat$cu_name)) |> 
  select(cu_name, adult_return_year, total_escape) |> 
  rename(year = adult_return_year) |> 
  arrange(cu_name, year)

#all totals####
all_totals <- left_join(return_total, spawner_total) |> 
  left_join(juv_total_abund |> select(cu_name, year = smolt_migration_year, presmolts, smolts)) |> 
  left_join(fishing_hatchery |> rename(year = adult_return_year)) |> 
  select(cu_name, year, spawners = total_escape, presmolts, smolts, returns = adult_returns, fishing, hatchery_spawners = hatchery_escape, fishing_prop) |> 
  ungroup() |> 
  mutate(cu_name = factor(cu_name, levels = cuid_dat$cu_name))


all_totals <- all_totals |> 
  mutate(year_id = year - year_start + 1) |> 
  left_join(cuid_dat |> select(cu_name, cu_id))

#environmental data####
env_dat <- read.csv("./data/env_data.csv") 

#relabel "Babine Lake" to "... (spawning channel)" to match smo_dat naming
env_dat$cu[env_dat$cu=="Babine Lake"] = "Babine Lake (spawning channel)"

env_select <- env_dat |>
  filter(cu %in% cuname_select) |> 
  select(year, cu_name = cu, 
         FTsre = FTsreMean,
         GCMTsre = gcm_FTsreMean,
         FTwre = FTwreMean,
         GCMTwre = gcm_FTwreMean,
         TEMPsoo = TEMPsooMean,
         TEMPwoo = TEMPwooMean,
         SSTsoo = SSTsooMean,
         SSTwoo = SSTwooMean,
         FTum = FTumMean,
         GCMTum = gcm_FTumMean,
         FDum = FDumMean,
         GCMDum = gcm_FDumMean) |> 
  mutate(cu_name = factor(cu_name, levels = cuid_dat$cu_name)) |> 
  left_join(cuid_dat |> select(cu_name, major_watershed))

#Set to NA 2024 winter rearing temps for all as they lack the jan-march 2025 data to integrate
env_select[which(env_select$year==2024),'FTwre'] = NA
env_select[which(env_select$year==2024),'GCMTwre'] = NA

#Set to NA summer rearing temps for Osoyoos and Wenatchee lake 2024 because they are based on unreliable observation data and seem unreasonably high
env_select[which(env_select$cu_name%in%c("Osoyoos Lake","Wenatchee Lake")&env_select$year==2024),'FTsre'] = NA
env_select[which(env_select$cu_name%in%c("Osoyoos Lake","Wenatchee Lake")&env_select$year==2024),'GCMTsre'] = NA

env_long <- env_select |> 
  gather(key = variable, value = value, -cu_name, -year, -major_watershed) |> 
  group_by(variable) |> 
  mutate(value_stnd = (value - mean(value, na.rm = TRUE))/sd(value, na.rm = TRUE)) |> #Global stdn
  group_by(variable, cu_name, major_watershed) |> 
  mutate(value_stnd_cu = (value - mean(value, na.rm = TRUE))/sd(value, na.rm = TRUE)) |> #CU stdn
  left_join(cuid_dat |> select(cu_name, cu_id)) |> 
  mutate(year_id = year - year_start + 1) |> 
  ungroup()

#CU stdn
env_stnd_cu <- env_long |> 
  select(-value, -value_stnd) |> 
  spread(key = variable, value = value_stnd_cu)

#Global stdn
env_stnd_global <- env_long |> 
  select(-value, -value_stnd_cu) |> 
  spread(key = variable, value = value_stnd)

#Open Ocean stdn
env_stnd_ocean_global <- env_long |> 
  filter(variable %in% c("TEMPwoo","TEMPsoo","SSTwoo","SSTsoo")) |> 
  select(year, value, variable) |> 
  distinct() |> 
  group_by(variable) |> 
  mutate(value_stnd = (value - mean(value, na.rm = TRUE))/sd(value, na.rm = TRUE)) |>  
  select(-value) |> 
  spread(key = variable, value = value_stnd) |>
  mutate(region="Ocean")

#coastal data####
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

#prepare data for stan####
return_age_prop <- adult_data |> 
  group_by(cu_name, adult_return_year) |> 
  mutate(age_prop = adult_abd/total_return) |> 
  mutate(age = paste("age", smolt_age, adult_age, sep = "_")) |> 
  select(cu_name, adult_return_year, age, age_prop) |> 
  spread(key = age, value = age_prop, fill = 0) |> 
  left_join(cuid_dat |> select(cu_name, cu_id)) |> 
  mutate(year_id = adult_return_year - year_start + 1) |> 
  filter(year_id > 0)

smolt_ageprop_from_adults <- adult_data |> 
  group_by(cu_name, smolt_migration_year) |> 
  mutate(smolt_year_total = sum(adult_abd), n_ages_smolt_year = n()) |> 
  filter(n_ages_smolt_year == 6) |> 
  group_by(cu_name, smolt_migration_year, smolt_age, smolt_year_total) |> 
  summarise(smolt_age_abd = sum(adult_abd)) |> 
  mutate(age_prop = smolt_age_abd/smolt_year_total) |> 
  mutate(smolt_age = paste0("age_", smolt_age, "_x")) |> 
  select(cu_name, smolt_migration_year, smolt_age, age_prop) |> 
  spread(key = smolt_age, value = age_prop, fill = 0) |> 
  mutate(source = "return_age")

smolt_ageprop_from_juv <- juv_data |> 
  select(cu_name, smolt_migration_year, smolt_abundance_age1, smolt_abundance_age2) |> 
  #filter(cu_name %in% c("Chilko Lake", "Sproat Lake", "Great Central Lake", "Tahltan Lake", "Tatsamenie Lake")) |>
  filter(cu_name %in% cuname_select) |> 
  group_by(cu_name, smolt_migration_year) |> 
  mutate(total = smolt_abundance_age1 + smolt_abundance_age2) |> 
  mutate(age_1_x = `smolt_abundance_age1` / total, age_2_x = `smolt_abundance_age2`/total) |> 
  select(cu_name, smolt_migration_year, age_1_x, age_2_x) |> 
  mutate(source = "juv_age") |> 
  filter(!is.na(age_2_x))

smolt_age_prop <- bind_rows(smolt_ageprop_from_juv, smolt_ageprop_from_adults) |> 
  arrange(cu_name, smolt_migration_year, desc(source == "juv_age")) %>%
  group_by(cu_name, smolt_migration_year) |> 
  distinct(cu_name, smolt_migration_year, .keep_all = TRUE) |> 
  left_join(cuid_dat |> select(cu_name, cu_id)) |> 
  mutate(year_id = smolt_migration_year - year_start + 1) |> 
  filter(year_id > 0) |> 
  mutate(n_eff = ifelse(source == "return_age", 10, 100))

fishing_prop_wide <- all_totals |> 
  select(year, cu_name, fishing_prop) |> 
  mutate(fishing_prop = ifelse(fishing_prop == 0, 1e-6, fishing_prop)) |> 
  spread(key = cu_name, value = fishing_prop, fill = -9999)

hatchery_spawners_wide <- all_totals |> 
  select(year, cu_name, hatchery_spawners) |> 
  spread(key = cu_name, value = hatchery_spawners, fill = 0)

#prior age class abundances####
juv_summary <- all_totals |> left_join(smolt_age_prop) |> 
  mutate(juveniles = ifelse(!is.na(smolts), smolts, presmolts)) |> 
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
  select(cu_name, year, returns, age_1_1:age_2_3) |> 
  gather(key = age, value = prop, age_1_1:age_2_3) |> 
  mutate(total = returns * prop) |> 
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


#stan data list####
model_data <- list(
  N_years = max(env_stnd_global$year_id),
  N_pop = length(cuid_dat$cu_name),
  N_rivers = max(cuid_dat$major_watershed_id),
  N_domains = max(cuid_dat$domain_id),
  N_regions = max(coastal_long$region_id),
  A_regions = unique(coastal_long$area),
  N_spawner_obs = sum(!is.na(all_totals$spawners)), 
  N_presmolt_obs = sum(!is.na(all_totals$presmolts)),
  N_smolt_obs = sum(!is.na(all_totals$smolts)),
  #N_return_obs = sum(!is.na(all_totals$returns)),
  
  N_smolt_age_obs = nrow(smolt_age_prop),
  N_return_age_obs = nrow(return_age_prop),
  
  N_fw_ages = 2,
  N_mar_ages = 3,
  
  #pop IDs
  river_id = cuid_dat$major_watershed_id, 
  domain_id = cuid_dat$domain_id,
  
  fishing_prop = fishing_prop_wide |> select(-year) |> data.matrix(),
  fishing_missid = matrix(as.numeric(fishing_prop_wide |> select(-year) |> data.matrix() == -9999), ncol = length(cuid_dat$cu_name)),
  
  hatchery_escapes = hatchery_spawners_wide |> select(-year) |> data.matrix(),
  
  #abundance observations
  spawner_obs = all_totals$spawners[!is.na(all_totals$spawners)],
  spawner_year = all_totals$year_id[!is.na(all_totals$spawners)], 
  spawner_pop = all_totals$cu_id[!is.na(all_totals$spawners)],
  
  presmolt_obs = all_totals$presmolts[!is.na(all_totals$presmolts)],
  presmolt_year = all_totals$year_id[!is.na(all_totals$presmolts)], 
  presmolt_pop = all_totals$cu_id[!is.na(all_totals$presmolts)],
  
  smolt_obs = all_totals$smolts[!is.na(all_totals$smolts)],
  smolt_year = all_totals$year_id[!is.na(all_totals$smolts)], 
  smolt_pop = all_totals$cu_id[!is.na(all_totals$smolts)],
  
  #return_obs = all_totals$returns[!is.na(all_totals$returns)],
  #return_year = all_totals$year_id[!is.na(all_totals$returns)], 
  #return_pop = all_totals$cu_id[!is.na(all_totals$returns)],
  
  #age class proportion observations
  smolt_age_prop = smolt_age_prop |> ungroup() |> select(age_1_x, age_2_x),
  smolt_age_n_eff = smolt_age_prop$n_eff,
  smolt_age_year = smolt_age_prop$year_id,
  smolt_age_pop = smolt_age_prop$cu_id,
  
  return_age_prop = return_age_prop |> ungroup() |> select(age_1_1:age_2_3),
  return_age_n_eff = rep(100, nrow(return_age_prop)),
  return_age_year = return_age_prop$year_id,
  return_age_pop = return_age_prop$cu_id,
  
  #initialization matrices
  s_priors_ln = mean_juv_ln |> select(-cu_name) |> data.matrix() |> t(), 
  s_priors_sd = sd_juv_ln |> select(-cu_name) |> data.matrix() |> t(),
  
  r_priors_ln = mean_ret_ln_array, 
  r_priors_sd = sd_ret_ln_array, 
  
  #Env data
  FTsre = env_to_matrix(env_stnd_global, "FTsre"),
  FTsre_missid = missID_matrix(env_stnd_global, "FTsre"),  
  N_FTsre_miss = sum(is.na(env_stnd_global$FTsre)),
  
  GCMTsre = env_to_matrix(env_stnd_global, "GCMTsre"),
  GCMTsre_missid = missID_matrix(env_stnd_global, "GCMTsre"),  
  #N_GCMTsre_miss = sum(is.na(env_stnd_global$GCMTsre)),
  GCMTsre_pop_mean = pop_means(env_stnd_global|>filter(year%in%2013:2024), "GCMTsre"),
  GCMTsre_pop_sd = pop_sd(env_stnd_global|>filter(year%in%2013:2024), "GCMTsre"),
  
  FTwre = env_to_matrix(env_stnd_global, "FTwre"),
  FTwre_missid = missID_matrix(env_stnd_global, "FTwre"),  
  N_FTwre_miss = sum(is.na(env_stnd_global$FTwre)),
  
  GCMTwre = env_to_matrix(env_stnd_global, "GCMTwre"),
  GCMTwre_missid = missID_matrix(env_stnd_global, "GCMTwre"),  
  #N_GCMTwre_miss = sum(is.na(env_stnd_global$GCMTwre)),
  GCMTwre_pop_mean = pop_means(env_stnd_global|>filter(year%in%2013:2024), "GCMTwre"),
  GCMTwre_pop_sd = pop_sd(env_stnd_global|>filter(year%in%2013:2024), "GCMTwre"),
  
  TEMP = env_to_matrix(coastal_long, "TEMP", spread_name = "region"),
  TEMP_missid = missID_matrix(coastal_long, "TEMP", spread_name = "region"),
  N_TEMP_miss = sum(is.na(coastal_long$TEMP)),
  
  SST = env_to_matrix(coastal_long, "SST", spread_name = "region"),
  SST_missid = missID_matrix(coastal_long, "SST", spread_name = "region"),
  N_SST_miss = sum(is.na(coastal_long$SST)),
  
  MLD = env_to_matrix(coastal_long, "MLD", spread_name = "region"),
  MLD_missid = missID_matrix(coastal_long, "MLD", spread_name = "region"),
  N_MLD_miss = sum(is.na(coastal_long$MLD)),
  
  MLDbccm = env_to_matrix(coastal_long, "MLDbccm", spread_name = "region"),
  MLDbccm_missid = missID_matrix(coastal_long, "MLDbccm", spread_name = "region"),
  N_MLDbccm_miss = sum(is.na(coastal_long$MLDbccm)),

  TEMPsoo = env_to_matrix(env_stnd_ocean_global, "TEMPsoo", spread_name = "region")[,1],
  TEMPsoo_missid = missID_matrix(env_stnd_ocean_global, "TEMPsoo", spread_name = "region")[,1],  
  N_TEMPsoo_miss = sum(is.na(env_stnd_ocean_global$TEMPsoo)),
  
  TEMPwoo = env_to_matrix(env_stnd_ocean_global, "TEMPwoo", spread_name = "region")[,1],
  TEMPwoo_missid = missID_matrix(env_stnd_ocean_global, "TEMPwoo", spread_name = "region")[,1],  
  N_TEMPwoo_miss = sum(is.na(env_stnd_ocean_global$TEMPwoo)),
  
  SSTsoo = env_to_matrix(env_stnd_ocean_global, "SSTsoo", spread_name = "region")[,1],
  #SSTsoo_missid = missID_matrix(env_stnd_ocean_global, "SSTsoo", spread_name="region")[,1],  
  #N_SSTsoo_miss = sum(is.na(env_stnd_ocean_global$SSTsoo)),
  
  SSTwoo = env_to_matrix(env_stnd_ocean_global, "SSTwoo", spread_name = "region")[,1],
  #SSTwoo_missid = missID_matrix(env_stnd_ocean_global, "SSTwoo", spread_name="region")[,1],  
  #N_SSTwoo_miss = sum(is.na(env_stnd_ocean_global$SSTwoo)),
  
  FTum = env_to_matrix(env_stnd_cu, "FTum"),
  FTum_missid = missID_matrix(env_stnd_cu, "FTum"),  
  N_FTum_miss = sum(is.na(env_stnd_cu$FTum)),
  
  GCMTum = env_to_matrix(env_stnd_cu, "GCMTum"),
  GCMTum_missid = missID_matrix(env_stnd_cu, "GCMTum"),  
  #N_GCMTum_miss = sum(is.na(env_stnd_cu$GCMTum)),
  GCMTum_pop_mean = pop_means(env_stnd_cu|>filter(year%in%2013:2024), "GCMTum"),
  GCMTum_pop_sd = pop_sd(env_stnd_cu|>filter(year%in%2013:2024), "GCMTum"),
  
  FDum = env_to_matrix(env_stnd_cu, "FDum"),
  FDum_missid = missID_matrix(env_stnd_cu, "FDum"),  
  N_FDum_miss = sum(is.na(env_stnd_cu$FDum)),
  
  GCMDum = env_to_matrix(env_stnd_cu, "GCMDum"),
  GCMDum_missid = missID_matrix(env_stnd_cu, "GCMDum"),  
  #N_GCMDum_miss = sum(is.na(env_stnd_cu$GCMDum)),
  GCMDum_pop_mean = pop_means(env_stnd_cu|>filter(year%in%2013:2024), "GCMDum"),
  GCMDum_pop_sd = pop_sd(env_stnd_cu|>filter(year%in%2013:2024), "GCMDum"),
  
  smolt_alpha1_prior_mean = 50,
  #smolt_alpha2_prior_mean = 1,
  #smolt_alpha1_year_sd_prior = 0.5,
  #smolt_alpha2_year_sd_prior = 0.5,
  
  return_alpha_ref_prior_mean = 0.015,
  
  smolt_beta_prior_mean = 7, 
  smolt_beta_prior_sd = 1
  #return_beta_prior_mean = 0.0001, #0.001 or 0.0001?
  #return_beta_prior_sd = 1
)

#initialization function####
init_fun <- function() list(
  a_returns_logit_mu = rep(qlogis(model_data$return_alpha_ref_prior_mean),6), 
  a_fry_mu_global_mu = 4,
  a_fry_mu_global_sigma = 0.5,
  sigma_returns = c(4,1,2.5,4,2,3.5),
  fishing_logit_mu = rep(0,model_data$N_pop),
  sigma_fishing = c(1.5,0.8,2,1.5,0.8,1.5,1.2,1.2,1.2,1,0.5,4,0.8),
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
  TEMP_SST_b = 1,
  MLD_MLDbccm_b = rep(1,model_data$N_regions),
  MLD_TEMP_b = rep(-1,model_data$N_regions),
  holdover_logit_mu = rep(qlogis(0.1),model_data$N_pop), #pop specific?
  holdover_logit_sigma = rep(0.5,model_data$N_pop),
  age2_survival_logit_mu = rep(qlogis(0.8),model_data$N_pop),
  age2_survival_logit_sigma = rep(0.5,model_data$N_pop)
)

#stan model####
model <- cmdstan_model("./sockeye_climate_model.stan")

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

#saveRDS(fit, "./output/sockeye_climate_model_full.rds")

#assess convergence####
worst_Rhat <-fit$summary() %>% 
  as.data.frame() %>% 
  mutate(rhat  = round(rhat , 3)) %>% 
  arrange(desc(rhat)) %>% 
  filter(!is.na(rhat))
mcmc_trace(fit$draws(), pars = worst_Rhat$variable[1:20], np = nuts_params(fit))

#Should predominantly be Rhat below 1.01, 
#and ess_bulk above 400 (1/10th of the sample size of 4 chains x 1000 iterations)
# worst_Rhat %>% 
#   ggplot(aes(x = ess_bulk , y = rhat))+
#   geom_point()+
#   geom_hline(yintercept = 1.01, lty = 2)+
#   geom_vline(xintercept = 400, lty = 2)+
#   labs(title=NULL,x="bulk-ESS",y="R-hat")+
#   theme_bw(base_size = 14)

fit$summary() %>% 
  as.data.frame() %>% 
  filter(!is.na(rhat))%>%
  ggplot(aes(x = ess_bulk , y = rhat))+
  geom_point()+
  geom_hline(yintercept = 1.01, lty = 2)+
  geom_vline(xintercept = 400, lty = 2)+
  labs(title=NULL,x="bulk-ESS",y="R-hat")+
  theme_bw(base_size = 14)

