library(dplyr)
library(data.table)
library(ggplot2)

results <- readRDS("min_results.RDS")
results[, initial_conditions:=factor(initial_conditions, levels=c("low", "medium", "high"))]
results[, waning:=factor(waning, levels=c("No waning", "Medium Waning", "Fast waning"))]
results[, seasonality:=factor(seasonality, levels=c("No seasonality", "20% seasonality", "30% seasonality"))]


data <- fread("https://raw.githubusercontent.com/folkehelseinstituttet/surveillance_data/master/covid19/data_covid19_hospital_by_time_latest.csv")
hosp_inc <- data[location_code=="norge" & date >= as.Date("2022-10-01") & date < as.Date("2022-11-03")]


level <- 0.975

sum_r_comb <- results %>% group_by(date) %>% summarise(hosp_inc_l=quantile(hosp_incidence, 1 - level),
                                                       hosp_inc_h=quantile(hosp_incidence, level),
                                                       inc_l=quantile(incidence, 1 - level),
                                                       inc_h=quantile(incidence, level),
                                                       hosp_l=quantile(hosp, 1 - level),
                                                       hosp_h=quantile(hosp, level),
                                                       resp_l=quantile(resp, 1 - level),
                                                       resp_h=quantile(resp, level)
                                                       )

ggplot(sum_r_comb) + geom_ribbon(aes(x=date, ymin=hosp_inc_l, ymax=hosp_inc_h), alpha=0.9, color="#393C61") + geom_point(aes(x=as.Date(date), y=n_hospital_main_cause), data=hosp_inc, color="red", size=4) + xlab("Time")  +ylab("Hospital Incidence") + theme_bw() + theme(text = element_text(size=34)) 
ggsave("overall_hosp_inc.png", width=20, height = 10)

       
ggplot(sum_r_comb) + geom_ribbon(aes(x=date, ymin=inc_l, ymax=inc_h), alpha=0.9, color="#393C61")  + xlab("Time")  +ylab("Incidence") + theme_bw() + theme(text = element_text(size=30))
ggsave("overall_inc.png", width=20, height = 10)

ggplot(sum_r_comb %>% filter(date > as.Date("2022-10-25"))) + geom_ribbon(aes(x=date, ymin=hosp_l, ymax=hosp_h), alpha=0.9, color="#393C61") + xlab("Time")  +ylab("Hospital Prevalence") + theme_bw() + theme(text = element_text(size=30))
ggsave("overall_hosp.png", width=20, height = 10)

ggplot(sum_r_comb %>% filter(date > as.Date("2022-10-25"))) + geom_ribbon(aes(x=date, ymin=resp_l, ymax=resp_h), alpha=0.9, color="#393C61") + xlab("Time")  +ylab("Respiratory Prevalence") + theme_bw() + theme(text = element_text(size=30)) 
ggsave("overall_resp.png", width=20, height = 10)



sum_r <- results %>% group_by(initial_conditions, waning, seasonality, date, severity) %>% summarise(hosp_inc_l=quantile(hosp_incidence, 1 - level),
                                                                                          hosp_inc_h=quantile(hosp_incidence, level),
                                                                                          inc_l=quantile(incidence, 1 - level),
                                                                                          inc_h=quantile(incidence, level))
                                                                                          

ggplot(sum_r %>% filter(initial_conditions=="medium" & waning=="No waning" & severity=="Medium severity")) + geom_ribbon(aes(x=date, ymin=hosp_inc_l, ymax=hosp_inc_h), alpha=0.9) + geom_point(aes(x=as.Date(date), y=n_hospital_main_cause), data=hosp_inc, color="red", size=3) + xlab("Time")  +ylab("Hospital Incidence") + theme_bw() + theme(text = element_text(size=30)) + facet_wrap(.~seasonality)
ggsave("seasonality.png", width=20, height = 10)



ggplot(sum_r %>% filter(initial_conditions=="medium" & seasonality=="No seasonality" & severity=="Medium severity")) + geom_ribbon(aes(x=date, ymin=hosp_inc_l, ymax=hosp_inc_h), alpha=0.9) + geom_point(aes(x=as.Date(date), y=n_hospital_main_cause), data=hosp_inc, color="red", size=3) + facet_grid(.~waning) + xlab("Time")  +ylab("Hospital Incidence") + theme_bw() + theme(text = element_text(size=30))
ggsave("waning.png", width=20, height = 10)


ggplot(sum_r %>% filter(waning == "Medium Waning" & seasonality=="20% seasonality" & severity=="Medium severity")) + geom_ribbon(aes(x=date, ymin=hosp_inc_l, ymax=hosp_inc_h), alpha=0.9) + geom_point(aes(x=as.Date(date), y=n_hospital_main_cause), data=hosp_inc, color="red", size=3) + facet_grid(.~initial_conditions) + xlab("Time")  +ylab("Hospital Incidence") + theme_bw() + theme(text = element_text(size=30))
ggsave("initial_conds.png", width=20, height = 10)

sum_r$severity <- factor(sum_r$severity, levels=c("Low severity", "Medium severity", "High severity"))
ggplot(sum_r %>% filter(waning == "Medium Waning" & seasonality=="20% seasonality" & initial_conditions=="medium")) + geom_ribbon(aes(x=date, ymin=hosp_inc_l, ymax=hosp_inc_h), alpha=0.9) + geom_point(aes(x=as.Date(date), y=n_hospital_main_cause), data=hosp_inc, color="red", size=3) + facet_grid(.~severity) + xlab("Time")  +ylab("Hospital Incidence") + theme_bw() + theme(text = element_text(size=30))
ggsave("severity.png", width=20, height = 10)


ggplot(sum_r %>% filter(initial_conditions=="medium")) + geom_ribbon(aes(x=date, ymin=hosp_inc_l, ymax=hosp_inc_h, fill=severity), alpha=0.9) + geom_point(aes(x=as.Date(date), y=n_hospital_main_cause), data=hosp_inc, color="red") + facet_grid(seasonality~waning) + xlab("Time")  +ylab("Hospital Incidence") + theme_bw() + theme(text = element_text(size=30))
ggsave("scenarios.png")


ggplot(sum_r %>% filter(seasonality=="20% seasonality")) + geom_ribbon(aes(x=date, ymin=hosp_inc_l, ymax=hosp_inc_h, fill=severity), alpha=0.7) + geom_point(aes(x=as.Date(date), y=n_hospital_main_cause), data=hosp_inc, color="red") + facet_grid(waning~initial_conditions) + xlab("Time")  +ylab("Hospital Incidence") + theme_bw() + theme(text = element_text(size=30))
ggsave("scenarios_ic.png", width=20, height = 15)



library(data.table)
results <- readRDS("min_variant_results.RDS")
results[, date:=time+ as.Date("2022-09-30")]

results[, variant_severity_nice:="Same severity"]
results[variant_severity == 2, variant_severity_nice:="Double severity"]
results$variant_severity_nice <- factor(results$variant_severity_nice, levels=c("Same severity", "Double severity"))

results[, variant_if_nice:="5% 1st Oct"]
results[variant_initial_frac > 0.07 , variant_if_nice:="10% 1st Oct"]

results$variant_if_nice <- factor(results$variant_if_nice, levels=c("5% 1st Oct", "10% 1st Oct"))



sum_r <- results %>% group_by(date, variant_severity_nice, variant_beta, variant_rr,
                              variant_if_nice) %>% summarise(hosp_inc_l=quantile(hosp_incidence, 1 - level),
                                                                  hosp_inc_h=quantile(hosp_incidence, level),
                                                                  hosp_inc_l_strain_1=quantile(hosp_incidence_strain_1, 1 - level),
                                                                  hosp_inc_h_strain_1=quantile(hosp_incidence_strain_1, level),
                                                                  hosp_inc_l_strain_2=quantile(hosp_incidence_strain_2, 1 - level),
                                                                  hosp_inc_h_strain_2=quantile(hosp_incidence_strain_2, level),
                                                                  inc_l=quantile(incidence, 1 - level),
                                                                  inc_h=quantile(incidence, level),
                                                                  inc_l_strain_1=quantile(incidence_strain_1, 1 - level),
                                                                  inc_h_strain_1=quantile(incidence_strain_1, level),
                                                                  inc_l_strain_2=quantile(incidence_strain_2, 1 - level),
                                                                  inc_h_strain_2=quantile(incidence_strain_2, level))





variant_hosp <- cbind(sum_r %>% select(date,  variant_severity_nice,variant_if_nice, variant_beta, variant_rr, hosp_inc_l_strain_1, hosp_inc_l_strain_2) %>% tidyr::pivot_longer(starts_with("hosp")) %>% mutate(value_lower=value),
                      value_higher=sum_r %>% select(date,  variant_if_nice,variant_severity_nice, variant_beta, variant_rr, hosp_inc_h_strain_1, hosp_inc_h_strain_2) %>% tidyr::pivot_longer(starts_with("hosp")) %>% pull(value))

variant_hosp <- variant_hosp %>%mutate(name=recode(name, hosp_inc_l_strain_1="Current variant", hosp_inc_l_strain_2="New variant"))

ggplot(variant_hosp %>% filter(variant_beta == 1.2 & abs(variant_rr-0.05)< 0.01 )) + geom_ribbon(aes(x=date, ymin=value_lower, ymax=value_higher, fill=name)) + facet_grid(variant_severity_nice~variant_if_nice) + theme_bw() + theme(text = element_text(size=30))+labs(fill="Variant") + ylab("New hospitalisations") + xlab("Date") + ggtitle("More transmisible variant")
ggsave("variant_more_trans.png", width=22, heigh=12)

ggplot(variant_hosp %>% filter(variant_beta == 1.0 & abs(variant_rr-0.35)< 0.01 )) + geom_ribbon(aes(x=date, ymin=value_lower, ymax=value_higher, fill=name)) + facet_grid(variant_severity_nice~variant_if_nice) + theme_bw() + theme(text = element_text(size=30))+labs(fill="Variant") + ylab("New hospitalisations") + xlab("Date") + ggtitle("Imunity evading variant")
ggsave("variant_im_ev.png", width=22, heigh=12)

ggplot(variant_hosp %>% filter(variant_beta == 1.2 & abs(variant_rr-0.35)< 0.01 )) + geom_ribbon(aes(x=date, ymin=value_lower, ymax=value_higher, fill=name)) + facet_grid(variant_severity_nice~variant_if_nice) + theme_bw() + theme(text = element_text(size=30))+labs(fill="Variant") + ylab("New hospitalisations") + xlab("Date") + ggtitle("More transmisible and immunity evading variant")
ggsave("variant_more_trans_im_ev.png", width=22, heigh=12)




## ggplot(res) + geom_line(aes(x=time + as.Date("2022-09-30"), y=hosp_incidence, group=sim)) + geom_point(aes(x=as.Date(date), y=n_hospital_main_cause), data=param_sets$data, color="red") 

## ggplot(res) + geom_line(aes(x=time, y=incidence, group=sim))

## ggplot(res) + geom_line(aes(x=t, y=tot_infected/5.4e6, group=sim))





## colnames(res)

## ggplot(res %>% select(t, sim, paste0("N_vac_", 1:10, "")) %>% tidyr::pivot_longer(starts_with("N"))) +
##   geom_line(aes(x=t, y=value, color=name, group=paste(sim, name)))



## # New variants

## ggplot(res) + geom_line(aes(x=time, y=incidence, group=sim))

## ggplot(res) + geom_line(aes(x=t, y=tot_infected/5.4e6, group=sim))





## colnames(res)

## ggplot(res %>% select(t, sim, paste0("N_vac_", 1:10, "")) %>% tidyr::pivot_longer(starts_with("N"))) +
##   geom_line(aes(x=t, y=value, color=name, group=paste(sim, name)))



## # New variants
