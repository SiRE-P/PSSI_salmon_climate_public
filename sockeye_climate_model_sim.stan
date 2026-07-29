data{
  int N_years;
  int N_pop;
  int N_rivers;
  int N_domains;
  //int N_presmolt_obs;
  int N_smolt_obs;
  int N_spawner_obs;
  
  int N_smolt_age_obs;
  int N_return_age_obs;
  
  int N_fw_ages;
  int N_mar_ages;
  
  //pop IDs
  array[N_pop] int<lower=1, upper=N_rivers> river_id;
  array[N_pop] int<lower=1, upper=N_domains> domain_id;
  
  //abundance observations
  
  //fry oberved, use "presmolt" in the  model to differentiate to "fry" in the BH
  //array[N_presmolt_obs] int presmolt_obs;
  //array[N_presmolt_obs] int presmolt_year;
  //[N_presmolt_obs] int presmolt_pop;

  //smolts observed
  array[N_smolt_obs] int smolt_obs;
  array[N_smolt_obs] int smolt_year;
  array[N_smolt_obs] int smolt_pop;
  
  //spawners observed
  array[N_spawner_obs] int spawner_obs;
  array[N_spawner_obs] int spawner_year;
  array[N_spawner_obs] int spawner_pop;
    
  //age class proportion observations
  array[N_smolt_age_obs, N_fw_ages] real smolt_age_prop;
  array[N_smolt_age_obs] int smolt_age_n_eff;
  array[N_smolt_age_obs] int smolt_age_year;
  array[N_smolt_age_obs] int smolt_age_pop;
  
  array[N_return_age_obs, N_fw_ages * N_mar_ages] real return_age_prop;
  array[N_return_age_obs] int return_age_n_eff;
  array[N_return_age_obs] int return_age_year;
  array[N_return_age_obs] int return_age_pop;
  
  //initialization priors
  matrix[N_fw_ages,N_pop] s_priors_ln; // priors for smolts unique to stocks and age class
  matrix[N_fw_ages,N_pop] s_priors_sd; // priors for smolts unique to stocks and age class
  
  array[N_fw_ages, N_mar_ages, N_pop] real r_priors_ln;
  array[N_fw_ages, N_mar_ages, N_pop] real r_priors_sd; 
  
  //environmental covariates  
  array[N_years, N_pop] real FTwre;    // FT winter rearing
  array[N_years, N_pop] int<lower=0> FTwre_missid;
  int<lower=0> N_FTwre_miss; 

  //array[N_years, N_pop] real GCMTwre;   //GCM hindcast data to impute missing FT data
  //array[N_years, N_pop] int<lower=0> GCMTwre_missid;
  //int<lower=0> N_GCMTwre_miss;
  //vector[N_pop] GCMTwre_pop_mean; 
  //vector[N_pop] GCMTwre_pop_sd;

  array[N_years, N_pop] real FTsre;    // FT summer rearing
  array[N_years, N_pop] int<lower=0> FTsre_missid; 
  int<lower=0> N_FTsre_miss; 

  //array[N_years, N_pop] real GCMTsre;
  //array[N_years, N_pop] int<lower=0> GCMTsre_missid;
  //int<lower=0> N_GCMTsre_miss;
  //vector[N_pop] GCMTsre_pop_mean; 
  //vector[N_pop] GCMTsre_pop_sd;

  array[N_years, N_rivers] real TEMP;
  array[N_years, N_rivers] int<lower=0> TEMP_missid;
  int<lower=0> N_TEMP_miss;

  //array[N_years, N_rivers] real SST;
  //array[N_years, N_rivers] int<lower=0> SST_missid;
  //int<lower=0> N_SST_miss;
  
  array[N_years, N_rivers] real MLD;
  array[N_years, N_rivers] int<lower=0> MLD_missid;
  int<lower=0> N_MLD_miss;

  //array[N_years, N_rivers] real MLDbccm;
  //array[N_years, N_rivers] int<lower=0> MLDbccm_missid;
  //int<lower=0> N_MLDbccm_miss;
  
  vector[N_years] TEMPwoo; 
  vector<lower=0>[N_years] TEMPwoo_missid; 
  int<lower=0> N_TEMPwoo_miss;
  
  //vector[N_years] SSTwoo; 
  
  vector[N_years] TEMPsoo; 
  vector<lower=0>[N_years] TEMPsoo_missid; 
  int<lower=0> N_TEMPsoo_miss;
  
  //vector[N_years] SSTsoo; 
  
  array[N_years, N_pop] real FTum;
  array[N_years, N_pop] int<lower=0> FTum_missid;
  int<lower=0> N_FTum_miss;

  //array[N_years, N_pop] real GCMTum;
  //array[N_years, N_pop] int<lower=0> GCMTum_missid;
  //int<lower=0> N_GCMTum_miss;
  //vector[N_pop] GCMTum_pop_mean;
  //vector[N_pop] GCMTum_pop_sd;

  array[N_years, N_pop] real FDum;
  array[N_years, N_pop] int<lower=0> FDum_missid;
  int<lower=0> N_FDum_miss;

  //array[N_years, N_pop] real GCMDum;
  //array[N_years, N_pop] int<lower=0> GCMDum_missid;
  //int<lower=0> N_GCMDum_miss;
  //vector[N_pop] GCMDum_pop_mean;
  //vector[N_pop] GCMDum_pop_sd;

  // Priors for alpha (smolts)
  real smolt_alpha1_prior_mean;       // mean for age-1 alpha (natural scale)
  //real smolt_alpha2_prior_mean;       // mean for age-2 alpha (natural scale)
  //real smolt_alpha1_year_sd_prior;    // SD prior for age-1 annual deviations
  //real smolt_alpha2_year_sd_prior;    // SD prior for age-2 annual deviations
  
  // Priors for alpha (returns)
  real return_alpha_ref_prior_mean;   // mean for reference (fw=1, mar=2)
  
  // Priors for beta parameters
  real smolt_beta_prior_mean;
  real smolt_beta_prior_sd;
  real return_beta_prior_mean;
  real return_beta_prior_sd;
}

parameters{
  //overdispersion
  real log_spawner_kappa_mu;
  real<lower=1e-6> log_spawner_kappa_sigma;
  vector[N_pop] log_spawner_kappa_z;
  vector<lower=1>[N_pop] phi_spawners;
  
  real log_smolt_kappa_mu;
  real<lower=1e-6> log_smolt_kappa_sigma;
  vector[N_pop] log_smolt_kappa_z;
  vector<lower=1>[N_pop] phi_smolts;
  
  //initialization abundances
  array[N_fw_ages, N_pop] real<lower=1e-6> smo_init; 
  array[N_fw_ages, N_mar_ages, N_pop] real<lower=1e-6> ret_init;
  
  //Process parameters
  // Alpha
  //Fry
  real a_fry_mu_global_mu;              
  real<lower=0> a_fry_mu_global_sigma;
  vector[N_pop] a_fry_mu_global_z;
  
  //holdover
  vector[N_pop] holdover_logit_mu;
  vector<lower=0>[N_pop] holdover_logit_sigma;
  matrix[N_years, N_pop] holdover_logit_z;
  
  //age2 survival
  vector[N_pop] age2_survival_logit_mu;
  vector<lower=0>[N_pop] age2_survival_logit_sigma;
  matrix[N_years, N_pop] age2_survival_logit_z;
  
  // annual deviations for fry productivity (log scale)
  array[N_pop, N_years] real a_fry_log_year_z;
  vector<lower=0>[N_pop] a_fry_log_year_sigma;
  
  //Returns
  // Returns (MVN on logit scale)
  vector[N_fw_ages * N_mar_ages] a_returns_logit_mu;                 // mean by (fw × mar)
  cholesky_factor_corr[N_fw_ages * N_mar_ages] L_corr_returns;       // correlation among age combos
  vector<lower=0>[N_fw_ages * N_mar_ages] sigma_returns;             // SD per combo
  matrix[N_fw_ages * N_mar_ages, N_pop] a_returns_logit_raw;         // standardized draws per population
  
  // annual deviations for return productivity (logit scale)
  array[N_fw_ages, N_mar_ages, N_pop, N_years] real a_returns_logit_year_z;
  vector<lower=0>[N_pop] a_returns_logit_year_sigma; //reduced dimensions without fw_age and mar_age
  
  //Beta
  real b_fry_mu; // mean density dependent effect
  real<lower=1e-6> b_fry_sigma;  // variance across stocks
  vector[N_pop] b_fry_z;  // z score error for b_returns per population
  
  //real b_returns_mu; // mean density dependent effect
  //real<lower=1e-6> b_returns_sigma;  // variance across stocks
  //vector[N_pop] b_returns_z;  // z score error for b_returns per population
    
  //environmental imputation parameters
  //vector[N_pop] FTwre_mui;
  //cholesky_factor_corr[N_pop] L_corr_FTwre;
  //vector<lower=0>[N_pop] sigma_FTwre;
  vector[N_FTwre_miss] FTwre_latent;
  //vector[N_pop] FTwre_latent_a;
  //real FTwre_GCMTwre_b;
  //real<lower=0> FTwre_latent_sigma;
  //vector[N_GCMTwre_miss] GCMTwre_latent;
  //array[N_years,N_pop] real GCMTwre_latent_pop;

  //vector[N_pop] FTsre_mui;
  //cholesky_factor_corr[N_pop] L_corr_FTsre;
  //vector<lower=0>[N_pop] sigma_FTsre;
  vector[N_FTsre_miss] FTsre_latent;
  //vector[N_pop] FTsre_latent_a;
  //real FTsre_GCMTsre_b;
  //real<lower=0> FTsre_latent_sigma;
  //vector[N_GCMTsre_miss] GCMTsre_latent;
  //array[N_years,N_pop] real GCMTsre_latent_pop;
  
  vector[N_TEMP_miss] TEMP_latent;
  //vector[N_SST_miss] SST_latent;
  vector[N_MLD_miss] MLD_latent;
  //vector[N_MLDbccm_miss] MLDbccm_latent;

  //real TEMP_latent_a;
  //real TEMP_SST_b;
  //real<lower=0> TEMP_latent_sigma;

  //vector[N_rivers] MLD_latent_a;
  //vector[N_rivers] MLD_MLDbccm_b;
  //vector<lower=0>[N_rivers] MLD_latent_sigma;

  //real TEMPwoo_latent_a;
  //real TEMPwoo_SSTwoo_b;
  //real<lower=0> TEMPwoo_latent_sigma;
  vector[N_TEMPwoo_miss] TEMPwoo_latent;
  
  //real TEMPsoo_latent_a;
  //real TEMPsoo_SSTsoo_b;
  //real<lower=0> TEMPsoo_latent_sigma;
  vector[N_TEMPsoo_miss] TEMPsoo_latent;
  
  //vector[N_pop] FTum_mui;
  //cholesky_factor_corr[N_pop] L_corr_FTum;
  //vector<lower=0>[N_pop] sigma_FTum;
  vector[N_FTum_miss] FTum_latent;
  //real FTum_latent_a;
  //real FTum_GCMTum_b;
  //real<lower=0> FTum_latent_sigma;
  //vector[N_GCMTum_miss] GCMTum_latent;
  //array[N_years,N_pop] real GCMTum_latent_pop;

  //vector[N_pop] FDum_mui;
  //cholesky_factor_corr[N_pop] L_corr_FDum;
  //vector<lower=0>[N_pop] sigma_FDum;
  vector[N_FDum_miss] FDum_latent;
  //real FDum_latent_a;
  //real FDum_GCMDum_b;
  //real<lower=0> FDum_latent_sigma;
  //vector[N_GCMDum_miss] GCMDum_latent;
  //array[N_years,N_pop] real GCMDum_latent_pop;

  //env covariates
  //matrix[N_pop,3] b_FTwre; //Three winters
  //matrix[N_pop,2] b_FTsre; //Two summers
  //vector[N_pop] b_TEMPcs;
  //vector[N_pop] b_MLDcs;
  //matrix[N_pop,3] b_TEMPwoo; //Three winters
  //matrix[N_pop,2] b_TEMPsoo; //Two summers
  //vector[N_pop] b_FTum;
  //vector[N_pop] b_FDum;

  //added hyper mu for freshwater variables as well
  real b_FTwre_hyper_mu; 
  real<lower=0> b_FTwre_hyper_sigma;  
  matrix[N_domains,3] b_FTwre_hyper_z;
  //vector[3] b_FTwre_mu; //3 winters;
  //vector<lower=0>[3] b_FTwre_sigma;
  matrix<lower=0>[N_domains,3] b_FTwre_sigma;
  matrix[N_pop,3] b_FTwre_z; 
  
  real b_FTsre_hyper_mu; 
  real<lower=0> b_FTsre_hyper_sigma;  
  matrix[N_domains,2] b_FTsre_hyper_z;
  //vector[2] b_FTsre_mu; //2 summers;
  //vector<lower=0>[2] b_FTsre_sigma;
  matrix<lower=0>[N_domains,2] b_FTsre_sigma;
  matrix[N_pop,2] b_FTsre_z; 

  real b_FTum_hyper_mu; 
  real<lower=0> b_FTum_hyper_sigma;  
  vector[N_domains] b_FTum_hyper_z;
  //real b_FTum_mu; 
  vector<lower=0>[N_domains] b_FTum_sigma;  
  vector[N_pop] b_FTum_z;

  real b_FDum_hyper_mu; 
  real<lower=0> b_FDum_hyper_sigma;  
  vector[N_domains] b_FDum_hyper_z;
  //real b_FDum_mu; 
  vector<lower=0>[N_domains] b_FDum_sigma;  
  vector[N_pop] b_FDum_z;

  //hyper mu for marine variables
  real b_TEMPcs_hyper_mu; 
  real<lower=0> b_TEMPcs_hyper_sigma;  
  vector[N_domains] b_TEMPcs_hyper_z;
  vector<lower=0>[N_domains] b_TEMPcs_sigma;  
  vector[N_pop] b_TEMPcs_z;

  real b_MLDcs_hyper_mu; 
  real<lower=0> b_MLDcs_hyper_sigma;  
  vector[N_domains] b_MLDcs_hyper_z;
  vector<lower=0>[N_domains] b_MLDcs_sigma;  
  vector[N_pop] b_MLDcs_z;

  real b_TEMPwoo_hyper_mu; 
  real<lower=0> b_TEMPwoo_hyper_sigma;  
  matrix[N_domains,3] b_TEMPwoo_hyper_z;
  matrix<lower=0>[N_domains,3] b_TEMPwoo_sigma;  
  matrix[N_pop,3] b_TEMPwoo_z;

  real b_TEMPsoo_hyper_mu; 
  real<lower=0> b_TEMPsoo_hyper_sigma;  
  matrix[N_domains,2] b_TEMPsoo_hyper_z;
  matrix<lower=0>[N_domains,2] b_TEMPsoo_sigma;  
  matrix[N_pop,2] b_TEMPsoo_z; 

  //covariate links
  //vector[N_rivers] MLDstnd_TEMPstnd_latent_a;
  //vector[N_rivers] MLDstnd_TEMPstnd_b;
  //vector<lower=0>[N_rivers] MLDstnd_TEMPstnd_latent_sigma;
  vector[N_pop] FTum_FDum_latent_a;
  vector[N_pop] FTum_FDum_b;
  vector<lower=0>[N_pop] FTum_FDum_latent_sigma;
}

transformed parameters{
 
  //environmental imputation
  //Marine
  array[N_years, N_rivers] real TEMP_complete;
  //array[N_years, N_rivers] real SST_complete;
  array[N_years, N_rivers] real MLD_complete;
  //array[N_years, N_rivers] real MLDbccm_complete;
  vector[N_years] TEMPwoo_complete;
  vector[N_years] TEMPsoo_complete;

  {
    int miss_i_temp = 1;
    //int miss_i_sst = 1;
    int miss_i_mld = 1;
    //int miss_i_mldbccm = 1;
    int miss_i_woo = 1;
    int miss_i_soo = 1;
    
    for (y in 1:N_years) {
      for (r in 1:N_rivers) {
        if (TEMP_missid[y, r] == 1) {
          TEMP_complete[y, r] = TEMP_latent[miss_i_temp];
          miss_i_temp += 1;
        } else {
          TEMP_complete[y, r] = TEMP[y, r];
        }
        
        //if (SST_missid[y, r] == 1) {
        //  SST_complete[y, r] = SST_latent[miss_i_sst];
        //  miss_i_sst += 1;
        //} else {
        //  SST_complete[y, r] = SST[y, r];
        //}
  
        if (MLD_missid[y, r] == 1) {
          MLD_complete[y, r] = MLD_latent[miss_i_mld];
          miss_i_mld += 1;
        } else {
          MLD_complete[y, r] = MLD[y, r];
        }
  
        //if (MLDbccm_missid[y, r] == 1) {
        //  MLDbccm_complete[y, r] = MLDbccm_latent[miss_i_mldbccm];
        //  miss_i_mldbccm += 1;
        //} else {
        //  MLDbccm_complete[y, r] = MLDbccm[y, r];
        //}
      }
      
      if (TEMPwoo_missid[y] == 1) {
        TEMPwoo_complete[y] = TEMPwoo_latent[miss_i_woo];
        miss_i_woo += 1;
      } else {
        TEMPwoo_complete[y] = TEMPwoo[y];
      }
      
      if (TEMPsoo_missid[y] == 1) {
        TEMPsoo_complete[y] = TEMPsoo_latent[miss_i_soo];
        miss_i_soo += 1;
      } else {
        TEMPsoo_complete[y] = TEMPsoo[y];
      }
    }
  }
  
  //Freshwater
  array[N_years, N_pop] real FTwre_complete;
  array[N_years, N_pop] real FTsre_complete;
  array[N_years, N_pop] real FTum_complete;
  array[N_years, N_pop] real FDum_complete;
  //array[N_years, N_pop] real GCMTwre_complete;
  //array[N_years, N_pop] real GCMTsre_complete;
  //array[N_years, N_pop] real GCMTum_complete;
  //array[N_years, N_pop] real GCMDum_complete;
  
  {
    int miss_i_wre = 1;
    int miss_i_sre = 1;
    int miss_i_tum  = 1;
    int miss_i_dum = 1;
    //int miss_i_gwre = 1;
    //int miss_i_gsre = 1;
    //int miss_i_gtum  = 1;
    //int miss_i_gdum = 1;
    
    for (y in 1:N_years) {
      for (p in 1:N_pop) {
        // FTwre
        if (FTwre_missid[y, p] == 1) {
          FTwre_complete[y, p] = FTwre_latent[miss_i_wre];
          miss_i_wre += 1;
        } else {
          FTwre_complete[y, p] = FTwre[y, p];
        }
        
        // FTsre
        if (FTsre_missid[y, p] == 1) {
          FTsre_complete[y, p] = FTsre_latent[miss_i_sre];
          miss_i_sre += 1;
        } else {
          FTsre_complete[y, p] = FTsre[y, p];
        }
        
        // FTum
        if (FTum_missid[y, p] == 1) {
          FTum_complete[y, p] = FTum_latent[miss_i_tum];
          miss_i_tum += 1;
        } else {
          FTum_complete[y, p] = FTum[y, p];
        }
        
        // FDum
        if (FDum_missid[y, p] == 1) {
          FDum_complete[y, p] = FDum_latent[miss_i_dum];
          miss_i_dum += 1;
        } else {
          FDum_complete[y, p] = FDum[y, p];
        }

        //if (GCMTwre_missid[y, p] == 1) {
        //  //GCMTwre_complete[y, p] = GCMTwre_latent[miss_i_gwre];
        //  GCMTwre_complete[y, p] = GCMTwre_latent_pop[y,p];
        //  //miss_i_gwre += 1;
        //} else {
        //  GCMTwre_complete[y, p] = GCMTwre[y, p];
        //}

        //if (GCMTsre_missid[y, p] == 1) {
        //  //GCMTsre_complete[y, p] = GCMTsre_latent[miss_i_gsre];
        //  GCMTsre_complete[y, p] = GCMTsre_latent_pop[y,p];
        //  //miss_i_gsre += 1;
        //} else {
        //  GCMTsre_complete[y, p] = GCMTsre[y, p];
        //}

        //if (GCMTum_missid[y, p] == 1) {
        //  //GCMTum_complete[y, p] = GCMTum_latent[miss_i_gtum];
        //  GCMTum_complete[y, p] = GCMTum_latent_pop[y,p];
        //  //miss_i_gtum += 1;
        //} else {
        //  GCMTum_complete[y, p] = GCMTum[y, p];
        //}

        //if (GCMDum_missid[y, p] == 1) {
        //  //GCMDum_complete[y, p] = GCMDum_latent[miss_i_gdum];
        //  GCMDum_complete[y, p] = GCMDum_latent_pop[y,p];
        //  //miss_i_gdum += 1;
        //} else {
        //  GCMDum_complete[y, p] = GCMDum[y, p];
        //}
      }
    }
  }
  
  //non-centred spawner parameters
  //vector[N_pop] phi_spawners;
  //for (p in 1:N_pop) {
  //  phi_spawners[p]   = fmin(1e3, fmax(1e-6, 1 / exp(log_spawner_kappa_mu + log_spawner_kappa_sigma * log_spawner_kappa_z[p])));
  //}
  
  //non-centered smolt parameters
  //vector[N_pop] phi_smolts;
  //for (p in 1:N_pop) {
  //  phi_smolts[p]   = fmin(1e3, fmax(1e-6, 1 / exp(log_smolt_kappa_mu + log_smolt_kappa_sigma * log_smolt_kappa_z[p])));
  //}
  
  //process model parameters
  //Alpha
  matrix[N_pop,3] b_FTwre; //Three winters
  matrix[N_pop,2] b_FTsre; //Two summers
  vector[N_pop] b_TEMPcs;
  vector[N_pop] b_MLDcs;
  matrix[N_pop,3] b_TEMPwoo; //Three winters
  matrix[N_pop,2] b_TEMPsoo; //Two summers
  vector[N_pop] b_FTum;
  vector[N_pop] b_FDum;

  //smolts
  array[N_years, N_pop] real a_fry;
  matrix[N_years, N_pop] holdover;
  matrix[N_years, N_pop] age2_survival;
  vector[N_pop] a_fry_log_year_mu;
  a_fry_log_year_mu = a_fry_mu_global_mu + a_fry_mu_global_sigma * a_fry_mu_global_z;
    
  //added hierarchy for env variables by domain and hyper
  
  matrix[N_domains,3] b_FTwre_mu;
  matrix[N_domains,2] b_FTsre_mu;
  
  for(d in 1:N_domains){
    for(w in 1:3){
        b_FTwre_mu[d,w] = b_FTwre_hyper_mu+b_FTwre_hyper_sigma*b_FTwre_hyper_z[d,w];}
    for(s in 1:2){
        b_FTsre_mu[d,s] = b_FTsre_hyper_mu+b_FTsre_hyper_sigma*b_FTsre_hyper_z[d,s];}
  }
  
  for (p in 1:N_pop){
    //for(w in 1:3){
    //  b_FTwre[p,w] = b_FTwre_mu[w]+b_FTwre_sigma[w]*b_FTwre_z[p,w];}
    //for(s in 1:2){
    //  b_FTsre[p,s] = b_FTsre_mu[s]+b_FTsre_sigma[s]*b_FTsre_z[p,s];}    
    for(w in 1:3){
      b_FTwre[p,w] = b_FTwre_mu[domain_id[p],w]+b_FTwre_sigma[domain_id[p],w]*b_FTwre_z[p,w];}
    for(s in 1:2){
      b_FTsre[p,s] = b_FTsre_mu[domain_id[p],s]+b_FTsre_sigma[domain_id[p],s]*b_FTsre_z[p,s];}

    for (y in 1:N_years) {
      real log_fw = a_fry_log_year_mu[p] + a_fry_log_year_sigma[p] * a_fry_log_year_z[p, y];
      if(y - 2 >= 1){
        log_fw += b_FTwre[p,1] * FTwre_complete[y - 2, p]; //incubation winter
        log_fw += b_FTsre[p,1] * FTsre_complete[y - 1, p] + b_FTwre[p,2] * FTwre_complete[y - 1, p];
      }
      a_fry[y, p] = exp(log_fw);
      holdover[y,p] = inv_logit(holdover_logit_mu[p] + holdover_logit_sigma[p] * holdover_logit_z[y,p]);
      age2_survival[y,p] = inv_logit(age2_survival_logit_mu[p] + age2_survival_logit_sigma[p] * age2_survival_logit_z[y,p]
      + b_FTsre[p,2] * FTsre_complete[y, p] + b_FTwre[p,3] * FTwre_complete[y, p]);
    }
  }
  
  
  //returns
  matrix[N_fw_ages * N_mar_ages, N_pop] a_returns_logit;
  a_returns_logit = rep_matrix(a_returns_logit_mu, N_pop) + diag_pre_multiply(sigma_returns, L_corr_returns) * a_returns_logit_raw;

  vector[N_domains] b_TEMPcs_mu;
  vector[N_domains] b_MLDcs_mu;
  matrix[N_domains,3] b_TEMPwoo_mu;
  matrix[N_domains,2] b_TEMPsoo_mu;
  vector[N_domains] b_FTum_mu;
  vector[N_domains] b_FDum_mu;
  
  b_TEMPcs_mu = b_TEMPcs_hyper_mu+b_TEMPcs_hyper_sigma*b_TEMPcs_hyper_z;
  b_MLDcs_mu = b_MLDcs_hyper_mu+b_MLDcs_hyper_sigma*b_MLDcs_hyper_z;
  for(d in 1:N_domains){
    for(w in 1:3){
        b_TEMPwoo_mu[d,w] = b_TEMPwoo_hyper_mu+b_TEMPwoo_hyper_sigma*b_TEMPwoo_hyper_z[d,w];}
    for(s in 1:2){
        b_TEMPsoo_mu[d,s] = b_TEMPsoo_hyper_mu+b_TEMPsoo_hyper_sigma*b_TEMPsoo_hyper_z[d,s];}
  }
  b_FTum_mu = b_FTum_hyper_mu+b_FTum_hyper_sigma*b_FTum_hyper_z;
  b_FDum_mu = b_FDum_hyper_mu+b_FDum_hyper_sigma*b_FDum_hyper_z;

  //convert back to 4D array for process model
  array[N_years, N_fw_ages, N_mar_ages, N_pop] real a_returns;

  //added hierarchy for env variables by domain, and hyper mu for marine variables
  for (p in 1:N_pop) {
    b_TEMPcs[p] = b_TEMPcs_mu[domain_id[p]]+b_TEMPcs_sigma[domain_id[p]]*b_TEMPcs_z[p];
    b_MLDcs[p] = b_MLDcs_mu[domain_id[p]]+b_MLDcs_sigma[domain_id[p]]*b_MLDcs_z[p];
    b_TEMPsoo[p,1] = b_TEMPsoo_mu[domain_id[p],1]+b_TEMPsoo_sigma[domain_id[p],1]*b_TEMPsoo_z[p,1];
    b_TEMPsoo[p,2] = b_TEMPsoo_mu[domain_id[p],2]+b_TEMPsoo_sigma[domain_id[p],2]*b_TEMPsoo_z[p,2];
    b_TEMPwoo[p,1] = b_TEMPwoo_mu[domain_id[p],1]+b_TEMPwoo_sigma[domain_id[p],1]*b_TEMPwoo_z[p,1];
    b_TEMPwoo[p,2] = b_TEMPwoo_mu[domain_id[p],2]+b_TEMPwoo_sigma[domain_id[p],2]*b_TEMPwoo_z[p,2];
    b_TEMPwoo[p,3] = b_TEMPwoo_mu[domain_id[p],3]+b_TEMPwoo_sigma[domain_id[p],3]*b_TEMPwoo_z[p,3];
    b_FTum[p] = b_FTum_mu[domain_id[p]]+b_FTum_sigma[domain_id[p]]*b_FTum_z[p];
    b_FDum[p] = b_FDum_mu[domain_id[p]]+b_FDum_sigma[domain_id[p]]*b_FDum_z[p];
    for (y in 1:N_years) {
      for (fw_age in 1:N_fw_ages) {
        for (mar_age in 1:N_mar_ages) {
          int idx = (fw_age - 1) * N_mar_ages + mar_age;
          real log_mar = a_returns_logit[idx, p] +  a_returns_logit_year_sigma[p] * a_returns_logit_year_z[fw_age, mar_age, p, y] + 
          b_FTum[p] * FTum_complete[y , p] + b_FDum[p] * FDum_complete[y , p];
          if(y - mar_age >= 1){
            log_mar += b_TEMPcs[p] * TEMP_complete[y - mar_age, river_id[p]] + b_MLDcs[p] * MLD_complete[y - mar_age, river_id[p]] + b_TEMPwoo[p,1] * TEMPwoo_complete[y - (mar_age-1)];
            if(mar_age >= 2){
              for(m in 2:mar_age){
                log_mar += b_TEMPsoo[p,m-1] * TEMPsoo_complete[y - (mar_age - m + 1)] + b_TEMPwoo[p,m] * TEMPwoo_complete[y - (mar_age - m)];}
            }
          }
          a_returns[y, fw_age, mar_age, p] = inv_logit(log_mar);
        }
      }
    }
  }
  
  //Beta
  vector[N_pop] b_fry;
  //vector[N_pop] b_returns;
  b_fry = exp(b_fry_mu + b_fry_sigma * b_fry_z);
  //b_returns = exp(b_returns_mu + b_returns_sigma * b_returns_z);
  
  //process model
  array[N_years, N_pop] real fry;
  array[N_years, N_pop] real holdover_fry;
  //array[N_years, N_pop] real presmolts;
  
  array[N_years, N_fw_ages, N_pop] real smolts_age;
  array[N_years, N_pop] real smolts;
  
  array[N_years, N_fw_ages, N_mar_ages, N_pop] real  returns_age;
  array[N_years, N_pop] real returns;
  array[N_years, N_pop] real spawners;
  
  
  //BH
  for(p in 1:N_pop){
    for(y in 1:N_years){
      //initialize at 0
      smolts[y, p] = 0;
      returns[y, p] = 0;
      spawners[y, p] = 0;
      
      //spawner to smolt, modelling "fry" to represent smolt and presomlt counts
      int by = y - 2;
      if(by >= 1){
        fry[y,p] = (a_fry[y,p] * spawners[by,p]) / (1 + b_fry[p] * (spawners[by,p] / 1e6));
        smolts_age[y,1,p] = fry[y,p] * (1 - holdover[y,p]);
        holdover_fry[y,p]  = fry[y,p] * holdover[y,p];
        
        if (y > 1){
          smolts_age[y,2,p] = holdover_fry[y-1,p] * age2_survival[y-1,p];
        } else {
          smolts_age[y,2,p] = smo_init[2,p];
        }
        //calculate total presmolts and smolts
        //presmolts[y,p] = smolts_age[y,1,p] + holdover_fry[y,p] + smolts_age[y,2,p];
        smolts[y,p]    = smolts_age[y,1,p] + smolts_age[y,2,p];
      } else{
        fry[y,p] = smo_init[1,p];
        smolts_age[y,1,p] = smo_init[1,p];
        smolts_age[y,2,p] = smo_init[2,p];
        holdover_fry[y,p] = smo_init[1,p] * holdover[y,p];
        //presmolts[y,p]    = smo_init[1,p] + holdover_fry[y,p] + smo_init[2,p];
        smolts[y,p]       = smo_init[1,p] + smo_init[2,p];
      }
      
      //smolt to returns
      for(fw_age in 1:N_fw_ages){
        for(mar_age in 1:N_mar_ages){
          int sy = y - mar_age; //define smolt migration year
          if (sy >= 1) {
            returns_age[y, fw_age, mar_age, p] = fmin(1e9, fmax(1e-9,
            (a_returns[y, fw_age, mar_age, p] * smolts_age[sy, fw_age, p])));
          } else{
            returns_age[y, fw_age, mar_age, p] = ret_init[fw_age,mar_age,p];
          }
          //calculate total returns
          returns[y,p] += returns_age[y, fw_age, mar_age, p];
        }
      }
      spawners[y,p] = returns[y, p] * (1-0.2); //fishing proportion
    }
  }
}


model{
 
  //environmental data imputation
    
  //Impute missing TEMP data from SST
  TEMP_latent ~ normal(0, 1);
  //SST_latent  ~ normal(0,1);
  //TEMP_latent_a ~ normal(0, 0.1);
  //TEMP_SST_b ~ normal(1, 0.5);
  //TEMP_latent_sigma ~ lognormal(log(0.5),0.4); 
  //for (y in 1:N_years){
  //  for(r in 1:N_rivers){
  //    TEMP_complete[y,r] ~ normal(TEMP_latent_a + TEMP_SST_b * SST_complete[y,r], TEMP_latent_sigma);}}
  
  //Impute missing MLD data from MLDbccm
  MLD_latent ~ normal(0, 1);
  //MLDbccm_latent ~ normal(0,1);
  //MLD_latent_a ~ normal(0, 0.1);
  //MLD_MLDbccm_b ~ normal(1, 0.5);
  //MLD_latent_sigma ~ lognormal(log(0.5),0.4);
  //for (y in 1:N_years){
  //  for(r in 1:N_rivers){
  //    MLD_complete[y,r] ~ normal(MLD_latent_a[r] + MLD_MLDbccm_b[r] * MLDbccm_complete[y,r], MLD_latent_sigma[r]);}}
    
  //Imputate missing Open Ocean temp values from complete SST as a reference
  TEMPwoo_latent ~ normal(0, 1);
  //TEMPwoo_latent_a ~ normal(0, 0.1);
  //TEMPwoo_SSTwoo_b ~ normal(1, 0.5);
  //TEMPwoo_latent_sigma ~ exponential(1);
  //TEMPwoo_complete ~ normal(TEMPwoo_latent_a + TEMPwoo_SSTwoo_b * SSTwoo, TEMPwoo_latent_sigma);
    
  TEMPsoo_latent ~ normal(0, 1);
  //TEMPsoo_latent_a ~ normal(0, 0.1);
  //TEMPsoo_SSTsoo_b ~ normal(1, 0.5);
  //TEMPsoo_latent_sigma ~ exponential(1);
  //TEMPsoo_complete ~ normal(TEMPsoo_latent_a + TEMPsoo_SSTsoo_b * SSTsoo, TEMPsoo_latent_sigma);
    
  //L_corr_FTwre ~ lkj_corr_cholesky(2);
  //sigma_FTwre ~ exponential(1);
  //FTwre_mui ~ normal(0, 1);
  FTwre_latent ~ normal(0, 1);
  //GCMTwre_latent ~ normal(0,1);
  //for(p in 1:N_pop){GCMTwre_latent_pop[,p] ~ normal(GCMTwre_pop_mean[p],GCMTwre_pop_sd[p]);}
    
  // if (N_FTwre_miss > 0) {
    //   matrix[N_pop, N_pop] L_FTwre = diag_pre_multiply(sigma_FTwre, L_corr_FTwre);
    //   for (y in 1:N_years)
    //   target += multi_normal_cholesky_lpdf(to_vector(FTwre_complete[y, ]) | FTwre_mui, L_FTwre);
    // }

  //FTwre_latent_a ~ normal(0, 0.1);
  //FTwre_GCMTwre_b ~ normal(1, 0.5);
  //FTwre_latent_sigma ~ exponential(1);
  //for (y in 1:N_years){
  //  for(p in 1:N_pop){
  //    FTwre_complete[y,p] ~ normal(FTwre_latent_a[p] + FTwre_GCMTwre_b * GCMTwre_complete[y,p], FTwre_latent_sigma);}}
  

  //L_corr_FTsre ~ lkj_corr_cholesky(2);
  //sigma_FTsre ~ exponential(1);
  //FTsre_mui ~ normal(0, 1);
  FTsre_latent ~ normal(0, 1);
  //GCMTsre_latent ~ normal(0,1);
  //for(p in 1:N_pop){GCMTsre_latent_pop[,p] ~ normal(GCMTsre_pop_mean[p],GCMTsre_pop_sd[p]);}
      
  // if (N_FTsre_miss > 0) {
    //   matrix[N_pop, N_pop] L_FTsre = diag_pre_multiply(sigma_FTsre, L_corr_FTsre);
    //   for (y in 1:N_years)
    //   target += multi_normal_cholesky_lpdf(to_vector(FTsre_complete[y, ]) | FTsre_mui, L_FTsre);
    // }

  //FTsre_latent_a ~ normal(0, 0.1);
  //FTsre_GCMTsre_b ~ normal(1, 0.5);
  //FTsre_latent_sigma ~ exponential(1);
  //for (y in 1:N_years){
  //  for(p in 1:N_pop){
  //    FTsre_complete[y,p] ~ normal(FTsre_latent_a[p] + FTsre_GCMTsre_b * GCMTsre_complete[y,p], FTsre_latent_sigma);}}

  //L_corr_FTum ~ lkj_corr_cholesky(2);
  //sigma_FTum ~ exponential(1);
  //FTum_mui ~ normal(0, 1);
  FTum_latent ~ normal(0, 1);
  //GCMTum_latent ~ normal(0,1);
  //for(p in 1:N_pop){GCMTum_latent_pop[,p] ~ normal(GCMTum_pop_mean[p],GCMTum_pop_sd[p]);}
        
  //if (N_FTum_miss > 0) {
    //  matrix[N_pop, N_pop] L_FTum = diag_pre_multiply(sigma_FTum, L_corr_FTum);
    //  for (y in 1:N_years)
    //  target += multi_normal_cholesky_lpdf(to_vector(FTum_complete[y, ]) | FTum_mui, L_FTum);
    //}

  //FTum_latent_a ~ normal(0, 0.1);
  //FTum_GCMTum_b ~ normal(1, 0.5);
  //FTum_latent_sigma ~ exponential(1);
  //for (y in 1:N_years){
  //  for(p in 1:N_pop){
  //    FTum_complete[y,p] ~ normal(FTum_latent_a + FTum_GCMTum_b * GCMTum_complete[y,p], FTum_latent_sigma);}}

  //L_corr_FDum ~ lkj_corr_cholesky(2);
  //sigma_FDum ~ exponential(1);
  //FDum_mui ~ normal(0, 1);
  FDum_latent ~ normal(0, 1);
  //GCMDum_latent ~ normal(0,1);
  //for(p in 1:N_pop){GCMDum_latent_pop[,p] ~ normal(GCMDum_pop_mean[p],GCMDum_pop_sd[p]);}

  // if (N_FDum_miss > 0) {
    //   matrix[N_pop, N_pop] L_FDum = diag_pre_multiply(sigma_FDum, L_corr_FDum);
    //   for (y in 1:N_years)
    //   target += multi_normal_cholesky_lpdf(to_vector(FDum_complete[y, ]) | FDum_mui, L_FDum);
    // }

  //FDum_latent_a ~ normal(0, 0.1);
  //FDum_GCMDum_b ~ normal(1, 0.5);
  //FDum_latent_sigma ~ exponential(1);
  //for (y in 1:N_years){
  //  for(p in 1:N_pop){
  //    FDum_complete[y,p] ~ normal(FDum_latent_a + FDum_GCMDum_b * GCMDum_complete[y,p], FDum_latent_sigma);}}

  //Link coastal migration MLD and TEMP, and Freshwater Disc and Temp
  //MLDstnd_TEMPstnd_latent_a ~ normal(0, 0.1);
  //MLDstnd_TEMPstnd_b ~ normal(-0.5, 0.5);
  //MLDstnd_TEMPstnd_latent_sigma ~ lognormal(log(0.5),0.4);
  FTum_FDum_latent_a ~ normal(0, 0.1);
  FTum_FDum_b ~ normal(-0.5, 0.5);
  FTum_FDum_latent_sigma ~ lognormal(log(0.5),0.4);

  //Temp > MLD
  //for (y in 1:N_years){
  //  for (r in 1:N_rivers){MLDstnd[y,r] ~ normal(MLDstnd_TEMPstnd_latent_a[r] + MLDstnd_TEMPstnd_b[r] * TEMPstnd[y,r], MLDstnd_TEMPstnd_latent_sigma[r]);}
  //  for (p in 1:N_pop){FTum_complete[y,p] ~ normal(FTum_FDum_latent_a[p] + FTum_FDum_b[p] * FDum_complete[y,p], FTum_FDum_latent_sigma[p]);}}

  //MLD > Temp
  for (y in 1:N_years){
    //for (r in 1:N_rivers){TEMPstnd[y,r] ~ normal(MLDstnd_TEMPstnd_latent_a[r] + MLDstnd_TEMPstnd_b[r] * MLDstnd[y,r], MLDstnd_TEMPstnd_latent_sigma[r]);}
    for (p in 1:N_pop){FTum_complete[y,p] ~ normal(FTum_FDum_latent_a[p] + FTum_FDum_b[p] * FDum_complete[y,p], FTum_FDum_latent_sigma[p]);}}

  //env priors
  //b_TEMPcs ~ normal(0, 0.5);
  //b_MLDcs ~ normal(0, 0.5);
  //to_vector(b_TEMPwoo) ~ normal(0, 0.5);
  //to_vector(b_TEMPsoo) ~ normal(0, 0.5);
  //to_vector(b_FTwre) ~ normal(0, 0.5);
  //to_vector(b_FTsre) ~ normal(0, 0.5);
  //b_FTum ~ normal(0, 0.5);
  //b_FDum ~ normal(0, 0.5);

  b_FTwre_hyper_mu ~ normal(0,0.5);
  b_FTwre_hyper_sigma ~ lognormal(log(0.25), 0.2);
  to_vector(b_FTwre_hyper_z) ~ std_normal();
  //b_FTwre_mu ~ normal(0,0.5);
  //b_FTwre_sigma ~ lognormal(log(0.5), 0.2);
  to_vector(b_FTwre_sigma) ~ lognormal(log(0.5), 0.2);
  to_vector(b_FTwre_z) ~ std_normal();

  b_FTsre_hyper_mu ~ normal(0,0.5);
  b_FTsre_hyper_sigma ~ lognormal(log(0.25), 0.2);
  to_vector(b_FTsre_hyper_z) ~ std_normal();
  //b_FTsre_mu ~ normal(0,0.5);
  //b_FTsre_sigma ~ lognormal(log(0.5), 0.2);
  to_vector(b_FTsre_sigma) ~ lognormal(log(0.5), 0.2);
  to_vector(b_FTsre_z) ~ std_normal();

  b_TEMPcs_hyper_mu ~ normal(0,0.5);
  b_TEMPcs_hyper_sigma ~ lognormal(log(0.5), 0.2);
  b_TEMPcs_hyper_z ~ std_normal();
  b_TEMPcs_sigma ~ lognormal(log(0.5), 0.2);
  b_TEMPcs_z ~ std_normal();

  b_MLDcs_hyper_mu ~ normal(0,0.5);
  b_MLDcs_hyper_sigma ~ lognormal(log(0.5), 0.2);
  b_MLDcs_hyper_z ~ std_normal();
  b_MLDcs_sigma ~ lognormal(log(0.5), 0.2);
  b_MLDcs_z ~ std_normal();

  b_TEMPwoo_hyper_mu ~ normal(0,0.5);
  b_TEMPwoo_hyper_sigma ~ lognormal(log(0.5), 0.2);
  to_vector(b_TEMPwoo_hyper_z) ~ std_normal();
  to_vector(b_TEMPwoo_sigma) ~ lognormal(log(0.5), 0.2);
  to_vector(b_TEMPwoo_z) ~ std_normal();

  b_TEMPsoo_hyper_mu ~ normal(0,0.5);
  b_TEMPsoo_hyper_sigma ~ lognormal(log(0.5), 0.2);
  to_vector(b_TEMPsoo_hyper_z) ~ std_normal();
  to_vector(b_TEMPsoo_sigma) ~ lognormal(log(0.5), 0.2);
  to_vector(b_TEMPsoo_z) ~ std_normal();

  b_FTum_hyper_mu ~ normal(0,0.5);
  b_FTum_hyper_sigma ~ lognormal(log(0.5), 0.2);
  b_FTum_hyper_z ~ std_normal();
  //b_FTum_mu ~ normal(0,0.5);
  b_FTum_sigma ~ lognormal(log(0.5), 0.2);
  b_FTum_z ~ std_normal();

  b_FDum_hyper_mu ~ normal(0,0.5);
  b_FDum_hyper_sigma ~ lognormal(log(0.5), 0.2);
  b_FDum_hyper_z ~ std_normal();
  //b_FDum_mu ~ normal(0,0.5);
  b_FDum_sigma ~ lognormal(log(0.5), 0.2);
  b_FDum_z ~ std_normal();

  //initialization
  // init returns - before time series
  for(p in 1:N_pop){
    for(fw_age in 1:N_fw_ages){
      smo_init[fw_age,p] ~ lognormal(s_priors_ln[fw_age,p], s_priors_sd[fw_age,p]);
      for(mar_age in 1:N_mar_ages){
        ret_init[fw_age,mar_age,p] ~ lognormal(r_priors_ln[fw_age,mar_age,p], r_priors_sd[fw_age,mar_age,p]);
      }
    }
  }
  
  //process model priors
  //Alpha
  //fry
  a_fry_mu_global_mu ~ normal(log(smolt_alpha1_prior_mean), 1);
  a_fry_mu_global_sigma ~ lognormal(log(0.5), 0.4);
  a_fry_mu_global_z ~ std_normal();
  
  //holdover
  holdover_logit_mu ~ normal(logit(0.1), 1);
  holdover_logit_sigma ~ lognormal(log(0.5), 0.3);
  to_vector(holdover_logit_z) ~ std_normal();
  
  //age2 survival
  age2_survival_logit_mu ~ normal(logit(0.8), 0.2);
  age2_survival_logit_sigma ~ lognormal(log(0.5), 0.3);
  to_vector(age2_survival_logit_z) ~ std_normal();
  
  // annual deviations for fry productivity (log scale)
  a_fry_log_year_sigma ~ lognormal(log(0.5), 0.4);
  for(p in 1:N_pop)
  to_vector(a_fry_log_year_z[p,]) ~ std_normal();
  
  //returns
  L_corr_returns ~ lkj_corr_cholesky(3);
  //sigma_returns ~ lognormal(log(0.5), 0.4);
  sigma_returns ~ lognormal(log(2.5), 0.4);
  to_vector(a_returns_logit_raw) ~ std_normal();
  
  // Priors on mean productivity (logit scale)
  for (fw_age in 1:N_fw_ages) {
    for (mar_age in 1:N_mar_ages) {
      int idx = (fw_age - 1) * N_mar_ages + mar_age;
      if (fw_age == 1 && mar_age == 2) {
        // reference combo
        a_returns_logit_mu[idx] ~ normal(logit(return_alpha_ref_prior_mean), 1);
      } else {
        a_returns_logit_mu[idx] ~ normal(0, 0.5); // centered near ref on logit scale
      }
    }
  }
  
  a_returns_logit_year_sigma ~ lognormal(log(0.5), 0.2); //simplified dimensions without fw_age and mar_age
  
  for (fw_age in 1:N_fw_ages){
    for(mar_age in 1:N_mar_ages){
      for(p in 1:N_pop){
        //a_returns_logit_year_sigma[fw_age, mar_age, p] ~ lognormal(log(0.5), 0.4);
        to_vector(a_returns_logit_year_z[fw_age, mar_age, p,]) ~ std_normal();
      }
    }
  }
  
  //Beta
  b_fry_mu ~ normal(log(smolt_beta_prior_mean), smolt_beta_prior_sd);
  //b_fry_sigma ~ lognormal(log(0.5), 0.4);
  b_fry_sigma ~ lognormal(log(1), 0.4);
  b_fry_z ~ std_normal();
  
  //b_returns_mu ~ normal(log(return_beta_prior_mean), return_beta_prior_sd);
  //b_returns_sigma ~ lognormal(log(0.5), 0.4);
  //b_returns_z ~ std_normal();
  
  //observation model
  //spawner abundance
  //priors
  log_spawner_kappa_mu ~ normal(log(0.1), 0.3); 
  //log_spawner_kappa_mu ~ normal(log(0.05), 0.3); 
  //log_spawner_kappa_sigma ~ normal(0, 0.2);
  log_spawner_kappa_sigma ~ lognormal(log(0.5), 0.4);
  log_spawner_kappa_z ~ std_normal();
  phi_spawners ~ lognormal(log(50),0.3);
  
  for(i in 1:N_spawner_obs){
    spawner_obs[i] ~ neg_binomial_2(spawners[spawner_year[i], spawner_pop[i]], phi_spawners[spawner_pop[i]]);
  }
  
  //smolt abundance
  //priors
  log_smolt_kappa_mu ~ normal(log(0.05), 0.3); 
  //log_smolt_kappa_sigma ~ normal(0, 0.3);
  log_smolt_kappa_sigma ~ lognormal(log(0.5), 0.4);
  log_smolt_kappa_z ~ std_normal();
  phi_smolts ~ lognormal(log(50),0.3);

  //total presmolt likelihood
  //for(i in 1:N_presmolt_obs){
  //  presmolt_obs[i] ~ neg_binomial_2(presmolts[presmolt_year[i], presmolt_pop[i]], phi_smolts[presmolt_pop[i]]);
  //}

  //total smolts likelihood
  for(i in 1:N_smolt_obs){
    smolt_obs[i] ~ neg_binomial_2(smolts[smolt_year[i], smolt_pop[i]], phi_smolts[smolt_pop[i]]);
  }
  
  //age prop likelihood
  for (i in 1:N_smolt_age_obs) {
    int y = smolt_age_year[i];
    int p = smolt_age_pop[i];
    vector[N_fw_ages] logits;
    for (fw_age in 1:N_fw_ages){
      logits[fw_age] = log(smolts_age[y, fw_age, p] + 1e-9);
    }
    
    vector[N_fw_ages] pi = softmax(logits);
    target += smolt_age_n_eff[i] * dot_product(to_vector(smolt_age_prop[i]), log(pi));
  }
  
  //returns
  //age prop likelihood
  for (i in 1:N_return_age_obs) {
    int y = return_age_year[i];
    int p = return_age_pop[i];
    vector[N_fw_ages * N_mar_ages] logits;
    for(fw_age in 1:N_fw_ages){
      for (mar_age in 1:N_mar_ages){
        int idx = (fw_age - 1) * N_mar_ages + mar_age;
        logits[idx] = log(returns_age[y, fw_age, mar_age, p] + 1e-9);
      }
    }
    vector[N_fw_ages * N_mar_ages] pi = softmax(logits);
    target += return_age_n_eff[i] * dot_product(to_vector(return_age_prop[i]), log(pi));
  }
}
