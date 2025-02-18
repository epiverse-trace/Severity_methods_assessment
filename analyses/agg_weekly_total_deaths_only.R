
#Script where CFR is estimated from weekly aggregated data, comparing 3 methods:
# a) Allocation of cases to first day of the week
# b) Averaging incidence across days of the week
# c) Imputation method in EpiEstim
#Uses Marburg data to illustrate the example, corresponding to case study #6 in the Severiy applications,
#where we have access to weekly reports of cases and only total deaths
#*Also attempted to use this dataset for time varying cfr, but there's too few cases for this purpose, so using simulated data instead

#install.packages('EpiEstim', repos = c('https://mrc-ide.r-universe.dev', 'https://cloud.r-project.org'))
library(EpiEstim)
library(incidence2)
library(epiparameter)
library(tidyverse)

#Loading data and aggregating by week
Marburg_EqGuinea_linelist <- read_csv("Data/Marburg_EqGuinea_linelist.csv")
Marburg_EqGuinea_linelist$Onset_week <- as.Date(Marburg_EqGuinea_linelist$Onset_week)
Marburg_aggregated <- Marburg_EqGuinea_linelist %>% incidence(date_index = "Onset_week", interval = 7) %>% complete_dates()

plot(Marburg_aggregated)

#Bayesian confidence interval function
ci_bayes <- function(x, n, conf.level = 0.95) {
  # Use the beta distribution for Bayesian confidence intervals
  bin_out <- binom::binom.bayes(
    x = x,
    n = n,
    conf.level = conf.level,
    type = "highest"
  )

  # Format the output
  est <- c(bin_out$mean, bin_out$lower, bin_out$upper)
  ci <- paste0(round(100 * est[1]), "% (",
               round(100 * est[2]), "-",
               round(100 * est[3]), ")")
  return(ci)
}
#Setting a date mid-outbreak for real-time simulation
real_time <- "2023-03-09"

#Defining MVD onset to death delay
set.seed(1)
extract_param(type = "range", values = c(8, 2, 16), distribution = "gamma", samples = 77)
marburg_onset_death <- epiparameter(disease = "marburg",
                                    pathogen = NA,
                                    epi_name = "onset_to_death",
                                    create_prob_distribution(
                                      prob_distribution = "gamma",
                                      prob_distribution_params = c(shape = 2.095, scale = 4.513)))
plot(marburg_onset_death)

# Deaths up to real-time point: we don't really know when deaths happened and if using real-time point, deaths > cases
# therefore choosing to include deaths for cases that happened up until 21 days before real-time point (considering delay of onset-death)
# deaths_rt <- Marburg_EqGuinea_linelist %>% filter(Onset_week <= as.Date(real_time) - 21) %>% filter(Status == "dead") %>% summarise(deaths = length(Case_id))

# Actual number of deaths in real-time from situation report on 22/03/23
deaths_rt$deaths <- 27

### A) Allocation of cases to first day of the week ####

#Formatting data for CFR
Marburg_EqGuinea_linelist$Death_week <- Marburg_EqGuinea_linelist$Onset_week
MVD_cases_deaths <- incidence2::incidence(Marburg_EqGuinea_linelist, c("Onset_week","Death_week")) %>% complete_dates()
MVD_cases_deaths <- cfr::prepare_data(MVD_cases_deaths, cases_variable = "Onset_week", deaths_variable = "Death_week")

#Truncating data at real-time point
MVD_agg_rt <- MVD_cases_deaths[MVD_cases_deaths$date <= real_time,]

#Estimating known outcomes
known_outcomes_aggregated <- cfr::estimate_outcomes(MVD_agg_rt, delay_density = function(x) density(marburg_onset_death, x))
total_known_aggregated <- sum(known_outcomes_aggregated$estimated_outcomes)

#CFR of aggregated weekly data adjusted for delays:
aggregated_cfr <- ci_bayes(deaths_rt$deaths, total_known_aggregated)


#### B) Averaging method ####
#Estimating weekly incidence
weekly_incidence <- Marburg_EqGuinea_linelist %>%
  incidence2::incidence(date_index = "Onset_week", interval = 7) %>%
  mutate(Onset_week = str_extract(date_index, "^\\S+")) %>%
  mutate(Onset_week = as.Date(Onset_week))

#Create daily dates for each week and dividing cases evenly by day
averaged_cases <- weekly_incidence %>%
  rowwise() %>%
  mutate(cases = count / 7,
         date = list(seq(Onset_week, Onset_week + 6, by = "day"))) %>%
  unnest(date) %>%
  select(date, cases)

#Adding deaths column (which isn't actually used)
averaged_cases$deaths <- averaged_cases$cases

#Truncating data at real-time point
averaged_cases_rt <- averaged_cases[averaged_cases$date <= real_time,]

#Estimating known outcomes
known_outcomes_averaged <- cfr::estimate_outcomes(averaged_cases_rt, delay_density = function(x) density(marburg_onset_death, x))
total_known_averaged <- sum(known_outcomes_averaged$estimated_outcomes)

#CFR of aggregated weekly data adjusted for delays:
averaged_cfr <- ci_bayes(deaths_rt$deaths, total_known_averaged)

#### C) Imputation with EpiEstim ####

#Getting si for Marburg and formatting it to work with EpiEstim
Marburg_si <- epiparameter_db(disease = "Marburg", epi_name = "serial_interval", subset = is_parameterised)
Marburg_si <- discretise(Marburg_si)
wrap_si <- function(si) {
  domain <- seq(1L, to = si$prob_dist$qf(0.999), by = 1L)
  pmf <- si$prob_dist$d(domain)
  pmf[1] <- 0
  pmf <- pmf / sum(pmf)
  pmf
}
si_distr <- wrap_si(Marburg_si)

#Using estimate_R with dt=7 to reconstruct daily incidence
set.seed(1)
Rt_est <- estimate_R(incid = Marburg_aggregated$count,
                     dt = 7L,
                     dt_out = 7L,
                     recon_opt = "naive",
                     iter = 10L,
                     tol = 1e-6,
                     grid = list(precision = 0.001, min = -1, max = 1),
                     config = make_config(si_distr = si_distr),
                     method = "non_parametric_si")


plot(Rt_est, legend = FALSE)

#Reconstructed daily incidence data:
daily_incidence <- as.data.frame(Rt_est$I[-(1:6)])
names(daily_incidence) <- "cases"
#Length equal to no. of weeks in the og dataset (18)*7, i.e. 126 days
#The first 7 days have the same value for the incidence, presumably because, as it says on the vignette,
#the first week doesn't have previous data to reconstruct the incidence from
#"...the earliest the incidence reconstruction can start is at least the first day of the second aggregation window."

#Now using daily incidence to estimate CFR
#Creating data frame for {cfr} with daily incidence and the corresponding dates
Marburg_daily <- Marburg_EqGuinea_linelist %>% incidence(date_index = "Onset_week") %>% complete_dates() #Annoying that the length is 120 and the length of the incidence is 126! So removing first 5 days from incidence data
daily_incidence$date <- as.Date(Marburg_daily$date_index)
daily_incidence$deaths <- as.Date(Marburg_daily$date_index) #This column is needed but not actually used by estimate_outcomes

ggplot(daily_incidence, aes(x = date, y = cases)) +
  geom_point() + scale_y_continuous(breaks = seq(0, 5, by = 0.5), limits = c(0,5))

# Truncating data mid-outbreak
daily_incidence_rt <- daily_incidence[daily_incidence$date <= real_time,]

#Estimating known outcomes
known_outcomes_imputation <- cfr::estimate_outcomes(daily_incidence_rt, delay_density = function(x) density(marburg_onset_death, x))
total_known_imputation <- sum(known_outcomes_imputation$estimated_outcomes)

#The CFR adjusted for delays is:
reconstructed_cfr <- ci_bayes(deaths_rt$deaths, total_known_imputation)

# Total cases up to real-time point
total_cases_rt <- Marburg_EqGuinea_linelist %>% filter(Onset_week <= real_time) %>% summarise(length(Case_id))

#The CFR not adjusted for delays is:
naive_cfr <- ci_bayes(deaths_rt$deaths, total_cases_rt$`length(Case_id)`)

## True CFR at the end of the outbreak
total_MVD_deaths <- Marburg_EqGuinea_linelist %>% filter(Status == "dead") %>% summarise(total_deaths = length(Case_id))
total_cases <- Marburg_EqGuinea_linelist %>% summarise(length(Case_id))

outbreak_end_cfr <- ci_bayes(total_MVD_deaths$total_deaths, total_cases$`length(Case_id)`)

## Table comparing approaches: reconstructing daily incidence vs simply assigning cases to first day of the week

CFR_table <- data.table::data.table(
  Method = c("Reconstructed daily incidence",
             "Cases on first day of each week",
             "Average cases on each day",
             "Naive estimate in real-time",
             "CFR at the end of the outbreak"),
  `CFR (95% CI)` = c(reconstructed_cfr,
                            aggregated_cfr,
                            averaged_cfr,
                            naive_cfr,
                            outbreak_end_cfr)
)

write.table(CFR_table, quote = F, sep = ",", row.names = F)

# Conclusion: if you're only given weekly cases and total deaths, it is not worthy to reconstruct daily incidence, as we will use aggregated
# known outcomes by date X, using estimate_outcomes(). The method that simply aggregates all cases to the first day of the week offers the same results,
# which closely approximate the true CFR at the end of the outbreak (off by 1%).
# This is unless the cut-off point is in the middle of a week, especially close to the beginning of the week! In that case, if it's using the second approach, we will miss those cases and the CFR will be overestimated
# This is sometimes not noticeable when rounding up the number for the denominator
# In real life though, cases will be reported on a weekly basis?


#### Test: adding column with deaths to test time_varying CFR

# 1.  Adding column with dates of death (simulated as Marburg data didn't have this info)
o_d_params <- get_parameters(marburg_onset_death)
Marburg_EqGuinea_linelist$Death_week <- Marburg_EqGuinea_linelist$Onset_week +
                                                  rgamma(nrow(Marburg_EqGuinea_linelist), shape =o_d_params[1], scale = o_d_params[2])
Marburg_EqGuinea_linelist$Death_week <- as.Date(ifelse(Marburg_EqGuinea_linelist$Status == "dead", Marburg_EqGuinea_linelist$Death_week, NA))

# 2. Converting death data to incidence
weekly_deaths <- Marburg_EqGuinea_linelist %>% incidence(date_index = "Death_week", interval = 7) %>% complete_dates()
weekly_deaths$date_index <- str_extract(weekly_deaths$date_index, "^\\S+")
weekly_deaths$date_index <- as.Date(weekly_deaths$date_index)

# 3. Reconstructing also daily deaths with EpiEstim
set.seed(1)
Rt_est_deaths <- estimate_R(incid = weekly_deaths$count,
                            dt = 7L,
                            dt_out = 7L,
                            recon_opt = "naive",
                            iter = 10L,
                            tol = 1e-6,
                            grid = list(precision = 0.001, min = -1, max = 1),
                            config = make_config(si_distr = si_distr),
                            method = "non_parametric_si")

plot(Rt_est_deaths)

I_deaths <- as.data.frame(Rt_est_deaths$I[-c(1:6)])
names(I_deaths)= "deaths"

# 4. Creating dataframe with daily dates and deaths
complete_dates_death <- weekly_deaths %>% incidence("date_index", interval = 1) %>%
  complete_dates() %>%
  select(date = date_index) %>%
  mutate(date = as.Date(date)) # Extracting dates from weekly death data
daily_deaths <- bind_cols(complete_dates_death, I_deaths)

# 5. Merging with dataframe that contains daily incidence
daily_cases <- select(daily_incidence, date, cases)
daily_cases_deaths <- merge(daily_cases, daily_deaths, by = "date", all.x = T, all.y = T)
daily_cases_deaths <- replace(daily_cases_deaths, is.na(daily_cases_deaths), 0)

plot(daily_cases_deaths$date,daily_cases_deaths$cases,type="l", col = "blue")
lines(daily_cases_deaths$date,daily_cases_deaths$deaths, type = "l", col="red")

# 6. Using cfr to estimate time varying CFR
daily_cases_deaths[,2:3] <- round(daily_cases_deaths[,2:3]) # PROBLEM, cases are too low with Marburg data for this





