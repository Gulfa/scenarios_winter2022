library(dplyr)
library(data.table)
library(metapop)
library(metapopnorge)
library(ggplot2)

run_current <- TRUE
run_variant <- TRUE


initial_S_dist = matrix(0, nrow=9, ncol=10)
initial_S_dist[,10] <- 0.7
initial_S_dist[,1] <- 0.3

param_file <- "parameter_files/parameters_vaccination.xlsx"

calc_seas <- function(date, amount, filename="parameter_files/norm_pred.csv"){
  dat <- fread(filename)
  dat$out <- 1 - dat$seas*amount
  dat <- rbind(dat, dat)
  i_start <- as.numeric(date - lubridate::make_date(year(date), 1, 1))
  selection <- dat[i_start:(i_start + 365), ]
  selection$out <- selection$out / selection[1, out]
  return(selection$out)
}


get_params <- function(initial_S_dist, waning_time=rep(20,11), seas=0, rr_inf= seq(0.95, 0, by=-0.1),
                       severity=1, include_new_variant=NULL){
    params <- read_param_file_simplified(param_file)
    L <- 6*30.5*2
    n_vac <- length(rr_inf)
    n_strain <- 1
    N <- 9
    seasonality <- calc_seas(as.Date("2022-10-01"), amount=seas)

    inc <- 20/0.00076/severity
    I_ini <- 0.6*inc*3
    P_ini <- 0.6*inc*2
    A_ini <- 0.4*inc*5
    Ea_ini <- 0.4*inc*2
    Es_ini <- 0.6*inc*2
    
    age_groups <- get_age_groups()

    S_ini <- round(age_groups*initial_S_dist)

    beta_strain <- rep(1, n_strain)
    rr_hosp <- rep(1, 10)
    rr_death <- rep(1, 10)
    
    if(!is.null(include_new_variant)){
      n_strain <- 2
      beta_strain <- c(1, include_new_variant$beta)
      rr_inf <- c(rr_inf, include_new_variant$rr_inf)
      rr_hosp <- c(rep(1, 10), rep(include_new_variant$severity, 10))
      rr_death <- c(rep(1, 10), rep(include_new_variant$severity, 10))
    }
    
    I_ini=array(round(I_ini/90), dim=c(N,n_vac,n_strain))
    Ea_ini=array(round(Ea_ini/90), dim=c(N,n_vac,n_strain))
    Es_ini=array(round(Es_ini/90), dim=c(N,n_vac,n_strain))
    A_ini=array(round(A_ini/90), dim=c(N,n_vac,n_strain))
    P_ini=array(round(P_ini/90), dim=c(N,n_vac,n_strain))
    
    if(n_strain==2){
      I_ini[,,2] <- round(I_ini[,,2]*include_new_variant$initial_frac)
      Ea_ini[,,2] <- round(Ea_ini[,,2]*include_new_variant$initial_frac)
      Es_ini[,,2] <- round(Es_ini[,,2]*include_new_variant$initial_frac)
      A_ini[,,2] <- round(A_ini[,,2]*include_new_variant$initial_frac)
      P_ini[,,2] <- round(P_ini[,,2]*include_new_variant$initial_frac)

    }
    
    vac_pars=list(rr_inf = rr_inf, 
                  rr_hosp = rr_hosp,
                  rr_death = rr_death,
                  rr_icu = rep(1,10),
                  rr_los_hosp = rep(1,10),
                  rr_inf_asymp = rep(1,10),
                  rr_trans = rep(1,10)
                  )

    T_waning <- matrix(waning_time, nrow=N, ncol=n_vac, byrow=T)
    params <- c(params,
                list(
                  N_steps=L,
                  n_vac=n_vac,
                  n_strain=n_strain,
                  dt=0.5,
                  T_waning=T_waning,
                  vaccinations=array(0,dim=c(L, N, n_vac)),
                                        #    import_vec=rep(0, L),
                                        #    import_age_prio=rep(1,dim=(N,n_vac,1)),
                  beta_day=matrix(seasonality[1:L], ncol=N, nrow=L),
                                        #    vac_time_full_effect=array(14.0, N),
                  beta_strain=beta_strain,
                  cross_protection=matrix(0, ncol=n_strain, nrow=n_strain),
                  n=9,
                  S_ini=S_ini,
                  import_vec=array(0, dim=c(L,N, n_vac, n_strain)),
                  I_ini=I_ini,
                  I_imp_ini=array(0, dim=c(N,n_vac,n_strain)),
                  Ea_ini=Ea_ini,
                  Es_ini=Es_ini,
                  A_ini=A_ini,
                  P_ini=P_ini,
                  H_ini=array(0, dim=c(N,n_vac,n_strain)),
                  ICU_H_ini=array(0, dim=c(N,n_vac,n_strain)),
                  ICU_R_ini=array(0, dim=c(N,n_vac,n_strain)),
                  ICU_P_ini=array(0, dim=c(N,n_vac,n_strain)),
                  B_D_ini=array(0, dim=c(N,n_vac,n_strain)),
                  B_D_H_ini=array(0, dim=c(N,n_vac,n_strain)),
                  B_D_ICU_ini=array(0, dim=c(N,n_vac,n_strain)),
                  R_ini=array(0, dim=c(N,n_vac,n_strain)),
                  D_ini=array(0, dim=c(N,n_vac,n_strain)),
                  tot_infected_ini=array(0, dim=c(N,n_vac,n_strain)),
                  tot_hosp_ini=array(0, dim=c(N,n_vac,n_strain)),
                  tot_resp_ini=array(0, dim=c(N,n_vac,n_strain)),
                  tot_vac_ini=array(0, dim=c(N, n_vac)),
                  beta_norm=age_groups,
                  reg_pop_long=age_groups,
                  N_regions=1,
                  waning_immunity_vax = array(1000, dim=c(N,n_vac,n_strain)),
                  waning_inf = 1000,
                  age_groups=N,
                  beta_mode=1,
                  include_waning=2
                )
                )
    basic_params <- fix_params(params, N, n_vac, n_strain, vac_pars)
    params <- fix_beta_mode_params(basic_params)
}


fit_beta <- function(params, n=300, n_parts =10, n_threads=10, n_best=10){

  data <- fread("https://raw.githubusercontent.com/folkehelseinstituttet/surveillance_data/master/covid19/data_covid19_hospital_by_time_latest.csv")
  hosp_inc <- data[location_code=="norge" & date >= as.Date("2022-10-01") & date < as.Date("2022-11-03")]

  filter <- get_filter_tot_hosp(hosp_inc$n_hospital_main_cause, params, n_particles =n_parts, n_threads=1)
  beta_1 <- fix_beta_large(params, params$S_ini, params$I_ini, 1, beta=params$beta_day[1,], use_eig=TRUE)
  basic_params <- copy(params)
  
  check_params <- function(p, return_history=FALSE){
    params$beta_day <- basic_params$beta_day*p*beta_1
    run <- filter$run(save_history = return_history, pars = params)
    if(return_history){
      return(filter$history(sample.int(filter$n_particles, 1)))
    }
    return(run)
  }
  p <- truncnorm::rtruncnorm(n, 0.5, 3, 1, 0.5)
  liks <- parallel::mclapply(1:n, function(i) {
    check_params(p[i])}, mc.cores=n_threads, mc.preschedule = FALSE
    )

  args <- doBy::which.maxn(unlist(liks), n_best)
  mean(p[args])
  param_sets <- list()
  for(i in 1:length(args)){
    param_sets[[i]] <- modifyList(params, list(beta_day=basic_params$beta_day*p[args[i]]*beta_1))
  }
  return(list(param_sets=param_sets,
              data=hosp_inc))
}

  
  




diff_time <- function(VE1, VE2, w, A=0.9){
   log(VE1/A) / w - log(VE2/A) / w
}

get_waning_times <- function(w, A){
  rev(c(diff_time(0.95, 0.85, w, A=A),
    diff_time(0.85, 0.75, w, A=A),
    diff_time(0.75, 0.65, w, A=A),
    diff_time(0.65, 0.55, w, A=A),
    diff_time(0.55, 0.45, w, A=A),
    diff_time(0.45, 0.35, w, A=A),
    diff_time(0.35, 0.25, w, A=A),
    diff_time(0.25, 0.15, w, A=A),
    diff_time(0.15, 0.05, w, A=A)))
}

i <- 0
scenarios <- list()
for(initial_conditions in c("low", "medium", "high")){
  for(waning in c(0, 0.005, 0.0075)){
    for (season in c(0, 0.2, 0.3)){
      for(severity in c(0.7, 1, 1.3)){
        scenarios[[length(scenarios) +1 ]] <- list(
          initial_conditions = glue::glue("parameter_files/{initial_conditions}_proportion_in_protection_compartments_by_age_group.csv"),
          initial_conditions_nice = initial_conditions,
          waning=waning,
          waning_nice=list("0"= "No waning", "0.005"="Medium Waning", "0.0075" = "Fast waning")[[as.character(waning)]],
          season=season,
          severity=severity,
          severity_nice=list("0.7"="Low severity", "1"="Medium severity", "1.3"="High severity")[[as.character(severity)]],
          season_nice=list("0" = "No seasonality", "0.2"="20% seasonality", "0.3"="30% seasonality")[[as.character(season)]],
          name=paste(i, initial_conditions, waning, season, severity))
        i <- i + 1
      }
    }
  }
}

variant_scenarios <- list()
for(initial_conditions in c("medium")){
  for(waning in c(0.005)){
    for (season in c(0.2)){
      for(severity in c(1)){
        for(beta_variant in c(1, 1.2)){
          for(var_rr in list(seq(0.95, 0, by=-0.1),
                             c(1, 1, 1, 0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35))){
            for(sev in c(1, 2)){
              for(initial_frac in c(0.05, 0.1)){
                variant_scenarios[[length(variant_scenarios) +1 ]] <- list(
                  initial_conditions = glue::glue("parameter_files/{initial_conditions}_proportion_in_protection_compartments_by_age_group.csv"),
                  initial_conditions_nice = initial_conditions,
                  waning=waning,
                  waning_nice=list("0"= "No waning", "0.005"="Medium Waning", "0.0075" = "Fast waning")[[as.character(waning)]],
                  season=season,
                  severity=severity,
                  severity_nice=list("0.7"="Low severity", "1"="Medium severity", "1.3"="High severity")[[as.character(severity)]],
                  season_nice=list("0" = "No seasonality", "0.2"="20% seasonality", "0.3"="30% seasonality")[[as.character(season)]],
                  include_new_variant=list(
                    beta=beta_variant,
                    rr_inf=var_rr,
                    severity=sev,
                    initial_frac=initial_frac,
                    name=paste(beta_variant, min(var_rr), sev)
                  ),
                  name=paste(i, initial_conditions, waning, season, severity))
                i <- i + 1
              }
            }
          }
        }
      }
    }
  }
}






update_severity <- function(params, severity, change_icu_prob=3.8, change_los_hosp=2.3, base_sev=10){
    params$hosp_prob <- params$hosp_prob*severity
    params$icu_prob <- params$icu_prob*min((1 + (change_icu_prob-1)*(severity-1)/(base_sev - 1)), change_icu_prob)
    params$length_hosp <- params$length_hosp*min((1 + (change_los_hosp-1)*(severity-1)/(base_sev - 1)), change_los_hosp)
    params$prob_death_hosp <- params$prob_death_icu <- params$prob_death_non_hosp <- params$prob_death_non_hosp*severity
    return(params)

}


run_scenario <- function(scenario, n=300,n_threads=10, n_best=10, n_parts=10, sims_per_best=5 ){
  print(scenario$name)
  initial_S <- fread(scenario$initial_conditions)
#  initial_S <- rbind(initial_S, data.frame(V1=100, age_group=0, protection_rounded=0.1, proportion=0))
                                        #  initial_S <- initial_S %>% arrange(age_group, protection_rounded)
  initial_S[, protection_rounded:=as.character(protection_rounded)]
  initial_S <- as.matrix(initial_S %>% tidyr::pivot_wider(id_cols=age_group, names_from=protection_rounded, values_from=proportion, values_fill=0) %>% select (-age_group))
  
  R <- 1.0
  waning_time <- c(1e10, get_waning_times(scenario$waning, 0.9))
  if(waning==0){
    waning_time <- rep(1e10, 10)
  }
  new_params <- get_params(initial_S, waning_time=waning_time, seas=scenario$season, severity=scenario$severity,
                           include_new_variant=scenario$include_new_variant)
  new_params <- update_severity(new_params, scenario$severity)
  ## params <- new_params
  ## beta_1 <- fix_beta_large(params, params$S_ini, params$I_ini, 1, beta=params$beta_day[1,], use_eig=TRUE)
  ## new_params$beta_day <- new_params$beta_day * beta_1
  ## res <- run_params(new_params, 120)
  
  param_sets <- fit_beta(new_params, n=n, n_parts=n_parts, n_threads=n_threads, n_best=n_best)
  r <- run_param_sets(param_sets$param_sets, L=120, sims_per_best,1,n_threads, silent=FALSE) %>% mutate(name=scenario$name,
                                                                                                   waning=scenario$waning_nice,
                                                                                                   seasonality=scenario$season_nice,
                                                                                                   severity=scenario$severity_nice,
                                                                                                   initial_conditions=scenario$initial_conditions_nice,
                                                                                                   new_variant= if(is.null(scenario$include_new_variant)) "No" else scenario$include_new_variant$variant_name
                                                                                                   )
  if(!is.null(scenario$include_new_variant)){
    r <- r %>% mutate(variant_beta=scenario$include_new_variant$beta,
                      variant_severity=scenario$include_new_variant$severity,
                      variant_initial_frac=scenario$include_new_variant$initial_frac,
                      variant_rr=min(scenario$include_new_variant$rr_inf)
                      )
  }else{
    r <- r %>% mutate(variant_beta=NA,
                      variant_severity=NA,
                      variant_rr=NA,
                      variant_initial_frac=NA
                      )
    
  }
  return(r)
}




if(run_current){
  all_results <- parallel::mclapply(scenarios, function(x) run_scenario(x, n_threads=10, n=400, n_best=10, sims_per_best=5), mc.cores=20, mc.preschedule = F)
  
  
  results <- rbindlist(all_results, fill=TRUE)
  results[, sim:=paste(sim, waning, seasonality, initial_conditions)]
  results[, date:=as.Date("2022-09-30") + time]
  
  
  saveRDS(results[, .(date, sim, initial_conditions,
                      waning, seasonality, severity, name, hosp_incidence, incidence, hosp,tot_resp, resp, D)], "min_results.RDS")
  
  saveRDS(results, "all_results.RDS")
}


if(run_variant){
  all_results <- parallel::mclapply(variant_scenarios, function(x) run_scenario(x, n_threads=10, n=400, n_best=10, sims_per_best=5), mc.cores=20, mc.preschedule = F)
  
  
  results <- rbindlist(all_results, fill=TRUE)
  results[, sim:=paste(sim, waning, seasonality, initial_conditions)]
  results[, date:=as.Date("2022-09-30") + time]

  saveRDS(results[, .(date, sim, initial_conditions,
                      waning, seasonality, severity, name, variant_beta, variant_initial_frac, variant_rr,variant_severity,
                      hosp_incidence, incidence, hosp,tot_resp, resp, D,
                      hosp_incidence_strain_1, hosp_incidence_strain_2,
                      incidence_strain_1, incidence_strain_2)], "min_variant_results.RDS")

  
  saveRDS(results, "variant_results.RDS")
}






