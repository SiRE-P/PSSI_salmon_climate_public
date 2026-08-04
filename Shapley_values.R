###Shapely value calculation from the model output to compare contributions of environmental effects, random effects and abundance data###

library(ggplot2)
library(dplyr)
library(tidyr)
library(GGally)
library(ggpubr)
library(ggrepel)
library(scales)
library(cmdstanr)
library(posterior)

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


fit <- readRDS("./output/sockeye_climate_model_full.rds")

# dims: [draw, year, pop]
rv <- as_draws_rvars(fit$draws())
D <- dim(draws_of(rv$b_fry))[1]

data <- model_data
ny   <- data$N_years
np   <- data$N_pop
nfw  <- data$N_fw_ages
nmar <- data$N_mar_ages

spawners0 <- draws_of(rv$spawners)
smo_init <- draws_of(rv$smo_init)
ret_init <- draws_of(rv$ret_init)

#Freshwater parameters
a_fry_mu <- draws_of(rv$a_fry_log_year_mu)
a_fry_sd <- draws_of(rv$a_fry_log_year_sigma)
a_fry_z  <- draws_of(rv$a_fry_log_year_z)

hold_mu <- draws_of(rv$holdover_logit_mu)
hold_sd <- draws_of(rv$holdover_logit_sigma)
hold_z  <- draws_of(rv$holdover_logit_z)

age2_mu <- draws_of(rv$age2_survival_logit_mu)
age2_sd <- draws_of(rv$age2_survival_logit_sigma)
age2_z  <- draws_of(rv$age2_survival_logit_z)

b_FTwre <- draws_of(rv$b_FTwre)
b_FTsre <- draws_of(rv$b_FTsre)

b_fry     <- draws_of(rv$b_fry)
#b_returns <- draws_of(rv$b_returns)

#Marine parameters
a_ret_logit <- draws_of(rv$a_returns_logit)
a_ret_year_sigma <- draws_of(rv$a_returns_logit_year_sigma)
a_ret_year_z  <- draws_of(rv$a_returns_logit_year_z)

b_TEMPcs  <- draws_of(rv$b_TEMPcs)
b_MLDcs   <- draws_of(rv$b_MLDcs)
b_TEMPwoo <- draws_of(rv$b_TEMPwoo)
b_TEMPsoo <- draws_of(rv$b_TEMPsoo)
b_FTum    <- draws_of(rv$b_FTum)
b_FDum    <- draws_of(rv$b_FDum)

#env covariates
FTsre_complete <- draws_of(rv$FTsre_complete)
FTwre_complete <- draws_of(rv$FTwre_complete)
TEMPstnd <- draws_of(rv$TEMPstnd)
MLDstnd <- draws_of(rv$MLDstnd)
TEMPwoo_complete <- draws_of(rv$TEMPwoo_complete)
TEMPsoo_complete <- draws_of(rv$TEMPsoo_complete)
FTum_complete <- draws_of(rv$FTum_complete)
FDum_complete <- draws_of(rv$FDum_complete)

#counterfactual functions####

median_over_years <- function(x) {
  stopifnot(length(dim(x)) == 3)
  
  d <- dim(x)[1]
  y <- dim(x)[2]
  p <- dim(x)[3]
  
  # mean across years: result is [d, p]
  m <- apply(x, c(1, 3), median)
  
  # replicate back across years -> [d, y, p]
  array(m, dim = c(d, 1, p))[, rep(1, y), ]
}



median_over_years1 <- function(x) {
  stopifnot(length(dim(x)) == 2)
  
  m <- apply(x, 1, median)
  matrix(m, nrow = nrow(x), ncol = ncol(x))
}


counter_returns <- function(spawner_var = TRUE, env_var = TRUE, random = TRUE, years = 1981:2024, cu_name = cuid_datx$cu_name, draws = D) {
  if(spawner_var == TRUE){
    spawners = spawners0
  } else{
    spawners = median_over_years(spawners0)
  }
  
  if(random == TRUE){
    random_fry_z = a_fry_z
    random_mar_z = a_ret_year_z
  } else{
    random_fry_z = array(0, dim = dim(a_fry_z))
    random_mar_z = array(0, dim = dim(a_ret_year_z))
    hold_z = array(0, dim = dim(hold_z))
    age2_z = array(0, dim = dim(age2_z))
  }
  
  if(env_var == FALSE){
    FTwre_complete = median_over_years(FTwre_complete)
    FTsre_complete = median_over_years(FTwre_complete)
    FTum_complete = median_over_years(FTum_complete)
    FDum_complete = median_over_years(FDum_complete)
    TEMPstnd = median_over_years(TEMPstnd)
    MLDstnd = median_over_years(MLDstnd)
    TEMPwoo_complete = median_over_years1(TEMPwoo_complete)
    TEMPsoo_complete = median_over_years1(TEMPsoo_complete)
  }
  
  
  a_fry      <- array(0, c(draws, ny, np))
  fryBH      <- array(0, c(draws, ny, np))
  holdover   <- array(0, c(draws, ny, np))
  age2surv   <- array(0, c(draws, ny, np))
  smolts_age <- array(0, c(draws, ny, nfw, np))
  smolts     <- array(0, c(draws, ny, np))
  ret <- array(0, c(draws, ny, nfw, nmar, np))
  returns_all    <- array(0, c(draws, ny, np))
  mar_surv <- array(0, c(draws, ny, nfw, nmar, np))
  
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = draws, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  for (d in 1:draws) {
    for (p in 1:np) {
      for (y in 1:ny) {
        log_fw <- a_fry_mu[d, p] + a_fry_sd[d, p] * random_fry_z[d, p, y]
        if (y - 2 >= 1) {
          log_fw <- log_fw +
            b_FTwre[d, p, 1] * FTwre_complete[d, y - 2, p] +
            b_FTsre[d, p, 1] * FTsre_complete[d, y - 1, p] +
            b_FTwre[d, p, 2] * FTwre_complete[d, y - 1, p]
        }
        a_fry[d,y,p] <- exp(log_fw)
        holdover[d, y, p] <- plogis(hold_mu[d, p] + hold_sd[d, p] * hold_z[d, y, p])
        age2surv[d, y, p] <- plogis(age2_mu[d, p] + age2_sd[d, p] * age2_z[d, y, p] +
                                      b_FTsre[d, p, 2] * FTsre_complete[d, y, p] +
                                      b_FTwre[d, p, 3] * FTwre_complete[d, y, p]
        )
        # spawners to fry BH (using exogenous spawners)
        by <- y - 2
        if (by >= 1) {
          
          fryBH[d,y,p] <- (a_fry[d, y, p] * spawners[d, by, p]) /
            (1 + b_fry[d, p] * (spawners[d, by, p] / 1e6))
          
          smolts_age[d, y, 1, p] <- fryBH[d,y,p] * (1 - holdover[d, y, p])
          
          if (y > 1) {
            smolts_age[d, y, 2, p] <- (fryBH[d,y-1,p] * holdover[d, y-1, p]) * age2surv[d, y - 1, p]
          } else {
            smolts_age[d, y, 2, p] <- smo_init[d, 2, p]
          }
          
        } else {
          smolts_age[d, y, , p] <- smo_init[d,, p]
        }
        
        smolts[d, y, p] <- sum(smolts_age[d, y, , p])
        
        # --------------------------------------------------------------
        # Smolt -> returns
        # --------------------------------------------------------------
        for (fw in 1:nfw) {
          for (mar in 1:nmar) {
            
            sy <- y - mar
            if (sy >= 1) {
              
              idx <- (fw - 1) * nmar + mar
              
              log_mar <- a_ret_logit[d, idx, p] + a_ret_year_sigma[d, p] * random_mar_z[d, fw, mar, p, y] +
                b_FTum[d, p] * FTum_complete[d, y, p] +
                b_FDum[d, p] * FDum_complete[d, y, p]
              
              # Environmental effects applied only when allowed (as in Stan)
              if (y - mar >= 1) {
                
                log_mar <- log_mar +
                  b_TEMPcs[d, p] * TEMPstnd[d, y - mar, data$river_id[p]] +
                  b_MLDcs[d, p] * MLDstnd[d, y - mar, data$river_id[p]] +
                  b_TEMPwoo[d, p, 1] * TEMPwoo_complete[d, y - (mar - 1)]
                
                if (mar >= 2) {
                  for (m in 2:mar) {
                    log_mar <- log_mar + b_TEMPsoo[d, p, m - 1] *
                      TEMPsoo_complete[d, y - (mar - m + 1)] +
                      b_TEMPwoo[d, p, m] * TEMPwoo_complete[d, y - (mar - m)]
                  }
                }
              }
              
              mar_surv[d,y,fw,mar,p] <- plogis(log_mar)
              #smolts to returns BH
              #ret[d,y,fw,mar,p] <- (mar_surv[d,y,fw,mar,p] * smolts_age[d, sy, fw, p]) / (1 + b_returns[d, p] * (smolts_age[d, sy, fw, p] / 1e6))
              ret[d,y,fw,mar,p] <- (mar_surv[d,y,fw,mar,p] * smolts_age[d, sy, fw, p])
              
              # Stan-style bounding
              ret[d,y,fw,mar,p] <- min(1e9, max(1e-9, ret[d,y,fw,mar,p] ))
              
            } else{
              ret[d,y,fw,mar,p] <- ret_init[d,fw,mar,p]
            }
            
            returns_all[d, y, p] <- returns_all[d, y, p] + ret[d,y,fw,mar,p] 
            
          }
        }
      }
    }
    setTxtProgressBar(pb, d)
  }
  close(pb)
  
  returns.df <- data.frame(returns=as.vector(matrix(returns_all,nrow=draws)),draw=rep(1:draws,length(years)*length(cu_name)),year=rep(rep(years,each=draws),length(cu_name)), cu_name = rep(cu_name, each=draws*length(years)))
  returns.df <- left_join(returns.df, cuid_dat)
  
  return(returns.df)
}

#counterfactuals####
null.df <- counter_returns(spawner_var = FALSE, env_var = FALSE, random = FALSE, draws = D)
S.df <- counter_returns(spawner_var = TRUE, env_var = FALSE, random = FALSE, draws = D)
E.df <- counter_returns(spawner_var = FALSE, env_var = TRUE, random = FALSE, draws = D)
R.df <- counter_returns(spawner_var = FALSE, env_var = FALSE, random = TRUE, draws = D)
ER.df <- counter_returns(spawner_var = FALSE, env_var = TRUE, random = TRUE, draws = D)
SR.df <- counter_returns(spawner_var = TRUE, env_var = FALSE, random = TRUE, draws = D)
SE.df <- counter_returns(spawner_var = TRUE, env_var = TRUE, random = FALSE, draws = D)
SER.df <- counter_returns(spawner_var = TRUE, env_var = TRUE, random = TRUE, draws = D)

##Shapley####
#combine effects and corrected effects
shapley.df <- SER.df|>select(draw, year, domain, cu_name, fSER = returns) |>
  left_join(S.df|>select(fS = returns, draw, year, domain, cu_name)) |>
  left_join(E.df|>select(fE = returns, draw, year, domain, cu_name)) |>
  left_join(R.df|>select(fR = returns, draw, year, domain, cu_name)) |>
  left_join(SE.df|>select(fSE = returns, draw, year, domain, cu_name)) |>
  left_join(SR.df|>select(fSR = returns, draw, year, domain, cu_name)) |>
  left_join(ER.df|>select(fER = returns, draw, year, domain, cu_name)) |> 
  left_join(null.df|>select(f0 = returns, draw, year, domain, cu_name)) 

#for each effect correct for the effect of median(spawners), that isolates the effect of spawners, random or env across all possible combinations for one effect
#they are then weighted to the sum of one and to avoid bias of possible combinations for each effect
#normalizing the absolute effect by the total effect produces relative effects
shapley.df <- shapley.df|>mutate(
  phi_spawner =
    (1/3) * (fS   - f0)  +
    (1/6) * (fSE  - fE)  +
    (1/6) * (fSR  - fR)  +
    (1/3) * (fSER - fER),
  phi_env =
    (1/3) * (fE   - f0)  +
    (1/6) * (fSE  - fS)  +
    (1/6) * (fER  - fR)  +
    (1/3) * (fSER - fSR),
  phi_random =
    (1/3) * (fR   - f0)  +
    (1/6) * (fSR  - fS)  +
    (1/6) * (fER  - fE)  +
    (1/3) * (fSER - fSE),
  total = abs(phi_spawner) + abs(phi_random) + abs(phi_env),
  prop_spawner = abs(phi_spawner) / total,
  prop_env = abs(phi_env) / total,
  prop_random = abs(phi_random) / total
) |> 
  filter(total > 0) |> 
  filter(year > 1984)

ggplot()+stat_pointinterval(data=shapley.df, aes(x = year, y = phi_env),colour='dodgerblue')+
  geom_hline(yintercept = 0, lty = 2, color = "grey")+
  facet_wrap(~cu_name, scales = "free_y")+
  theme_bw(base_size=10)+
  scale_x_continuous(breaks=c(1985,1990,2000,2010,2020))+
  scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M"))+
  labs(y = "Environmental contribution to recruitment", x = "")

ggplot()+stat_pointinterval(data=shapley.df, aes(x = year, y = phi_random),colour='goldenrod2')+
  geom_hline(yintercept = 0, lty = 2, color = "grey")+
  facet_wrap(~cu_name, scales = "free_y")+
  theme_bw(base_size=10)+
  scale_x_continuous(breaks=c(1985,1990,2000,2010,2020))+
  scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M"))+
  labs(y = "Random effect contribution to recruitment", x = "")

ggplot()+stat_pointinterval(data=shapley.df, aes(x = year, y = phi_spawner),colour='salmon')+
  geom_hline(yintercept = 0, lty = 2, color = "grey")+
  facet_wrap(~cu_name, scales = "free_y")+
  theme_bw(base_size=10)+
  scale_x_continuous(breaks=c(1985,1990,2000,2010,2020))+
  scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M"))+
  labs(y = "Spawner abundance contribution to recruitment", x = "")

# shapley.df |> 
#   select(year, cu_name, draw, phi_env, phi_spawner, phi_random) |> 
#   gather(key = type, value = phi, -year, -cu_name, -draw) |> 
#   ggplot(aes(x = year, y = phi, color = type))+
#   stat_summary(fun = median, geom = "point")+  geom_hline(yintercept = 0, lty = 2, color = "grey")+
#   facet_wrap(~cu_name, scales = "free_y")+
#   theme_bw()

#yearly effects, merge with adult data to indicate years with and without actual observations
#observed returns
se_adult <- read.csv("./data/adult_table.csv", header = TRUE)
adult_data <- se_adult%>%filter(!en_route_mort=="included") |>
  filter(cu_name %in% cuname_select, adult_return_year>=1981) |>
  mutate(cu_name = factor(cu_name, levels = cuid_dat$cu_name)) |>
  select(cu_name, year=adult_return_year, spawners=total_escape, returns_1.1=total_returns_1.1, returns_1.2=total_returns_1.2, returns_1.3=total_returns_1.3, returns_2.1=total_returns_2.1, returns_2.2=total_returns_2.2, returns_2.3=total_returns_2.3)|>
  mutate(returns_total=rowSums(across(starts_with('returns_')),na.rm=TRUE))
adult_data$cu_name <- as.character(adult_data$cu_name)
adult_data$cu_name[adult_data$cu=="Babine Lake (spawning channel)"] = "Babine Lake"
adult_data$cu_name <- factor(adult_data$cu_name, levels=levelx)
adult_data$returns_total[which(adult_data$returns_total==0)] = NA
adult_data$domain <- "observed"
adult_data$pop <- adult_data$cu_name

shapley.year <- shapley.df|>select(pop=cu_name,year,phi_spawner,phi_env,phi_random)|>pivot_longer(-c(pop,year),values_to = "phi", names_to = "eff")|>group_by(pop,year,eff)|>summarize(mean_phi = mean(phi,na.rm=TRUE))
shapley.year_obs <- shapley.year|>left_join(adult_data|>select(pop,year,returns_total))|>mutate(obs=ifelse(is.na(returns_total),"no","yes"))|>mutate(eff=recode(eff,'phi_env'='Environmental','phi_random'='Random','phi_spawner'='Spawner'))
ggplot()+geom_point(data=shapley.year_obs,aes(x=year,y=mean_phi,colour=eff,shape=obs),size=3)+
  scale_colour_manual(values=c("Environmental" = "dodgerblue", "Random" = "goldenrod2","Spawner"="salmon"))+
  scale_shape_manual(values=c("no"=1,"yes"=19))+
  facet_wrap(~pop, scales = "free_y", ncol=3)+
  theme_bw(base_size=10)+
  scale_x_continuous(breaks=c(1985,1990,2000,2010,2020))+
  theme(legend.position='bottom')


#Proportion explained over the years
shapley.df_var <- shapley.df|>select(pop = cu_name,year,phi_spawner,phi_env,phi_random)|>group_by(pop)|>
  summarise(
    var_spawner = var(phi_spawner, na.rm = TRUE),
    var_env     = var(phi_env, na.rm = TRUE),
    var_random  = var(phi_random, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    total_var = var_spawner + var_env + var_random,
    prop_spawner = var_spawner / total_var,
    prop_env     = var_env     / total_var,
    prop_random  = var_random  / total_var
  )

shapley.df_varx <- shapley.df_var|>select(pop,Spawner=prop_spawner,Environmental=prop_env,Random=prop_random)|>pivot_longer(-c(pop),values_to = "prop", names_to = "eff")

shapley.df_varx <- shapley.df_varx|>mutate(pop_name_short = str_trim(str_remove(pop, "\\b[Ll]ake\\b")))|>mutate(pop_name_short=factor(pop_name_short,levels=levelxs))
shapley.df_varx$eff <- factor(shapley.df_varx$eff,levels=c("Environmental","Random", "Spawner"))
ggplot()+geom_bar(data=shapley.df_varx,aes(x=pop_name_short,y=prop,fill=eff),position='fill',stat='identity')+
  scale_fill_manual(values=c("Environmental" = "dodgerblue", "Random" = "goldenrod2","Spawner"="salmon"))+
  geom_hline(yintercept = c(0.25, 0.5, 0.75), lty = 3, color = "grey90")+
  scale_y_continuous(breaks=c(0,0.5,1))+
  theme_minimal(base_size=14)+
  labs(title=NULL,x=NULL,y="Interannual Variability",fill="Effects:")+
  theme(legend.position='bottom')+
  coord_cartesian(expand = FALSE)+
  theme(axis.text.x = element_text(angle = -90, hjust = 0.1, vjust=0.1))
