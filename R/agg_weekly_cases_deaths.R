
#Script where CFR is estimated from weekly aggregated data, comparing 3 methods:
# a) Allocation of cases to first day of the week
# b) Averaging incidence across days of the week
# c) Imputation method in EpiEstim
#Uses simulated data to illustrate the example, corresponding to case study #6 in the Severiy applications,
#where we have access to weekly reports of cases and deaths

library(simulist)
library(epiparameter)
library(cfr)
library(tidyverse)
library(incidence2)
library(EpiEstim)
library(gridExtra)

########### Option 1: Long delays and stable infection dynamics #########
# 1. Simulating linelist
# create contact distribution
contact_distribution <- epiparameter::epiparameter(
  disease = "COVID-19",
  epi_name = "contact distribution",
  prob_distribution = create_prob_distribution(
    prob_distribution = "pois",
    prob_distribution_params = c(mean = 2)
  )
)

# create infectious period
infectious_period <- epiparameter::epiparameter(
  disease = "COVID-19",
  epi_name = "infectious period",
  prob_distribution = create_prob_distribution(
    prob_distribution = "gamma",
    prob_distribution_params = c(shape = 1, scale = 1)
  )
)

# create onset to death
onset_death <- epiparameter_db(
  disease = "COVID-19",
  epi_name = "onset_to_death",
  single_epiparameter = TRUE)

# defining time-varying cfr
#config <- create_config(time_varying_death_risk = function(risk, time) risk * exp(-0.1 * time))

# simulating linelist
set.seed(10)
linelist <- sim_linelist(
  contact_distribution = contact_distribution,
  infectious_period = infectious_period,
  prob_infection = 0.5,
  onset_to_hosp = NULL,
  onset_to_death = onset_death,
  onset_to_recovery = NULL,
  non_hosp_death_risk = 0.2,
  outbreak_start_date = as.Date("2025-01-01"),
  outbreak_size = c(1000, 5000)

)

# 2. Converting to incidence and grouping by 7 days
# creating columns for deaths and recovery dates
linelist <- linelist %>%
  pivot_wider(
    names_from = outcome,
    values_from = date_outcome
  ) %>%
  rename(
    date_death = died,
  )

# converting to incidence
weekly_inc <- incidence(linelist,
                       date_index = c(cases = "date_onset", deaths = "date_death"),
                       interval = 7) %>%
  complete_dates()

# visualising cases, deaths and recoveries over time
plot(weekly_inc)

# 3. Method 1: CFR by assuming cases happen on day 1 of the week
# From weekly to daily incidence

weekly_cases <- weekly_inc %>%
  filter(!count == 0) %>%
  filter(count_variable == "cases") %>%
  mutate(date = str_extract(date_index, "^\\S+")) %>%
  mutate(cases = count) %>%
  select(date, cases)

daily_inc <- incidence(weekly_cases, date_index = "date", counts = "cases", interval = 1) %>%
  complete_dates()

weekly_deaths <- weekly_inc %>%
  filter(!count == 0) %>%
  filter(count_variable == "deaths") %>%
  mutate(date = str_extract(date_index, "^\\S+")) %>%
  mutate(deaths = count) %>%
  select(date, deaths)

daily_deaths <- incidence(weekly_deaths, date_index = "date", counts = "deaths", interval = 1) %>%
  complete_dates()

# merging both into a single dataset and formatting for CFR
daily_cases_deaths <- merge(daily_inc, daily_deaths, by = "date_index", all = TRUE)

df_transformed <- daily_cases_deaths %>%
  mutate(
    cases = ifelse(count_variable.x == "cases", count.x, 0),
    deaths = ifelse(count_variable.y == "deaths", count.y, 0)
  ) %>%
  mutate(
    cases = ifelse(is.na(cases), 0, cases),
    deaths = ifelse(is.na(deaths), 0, deaths)
  ) %>%
  select(date = date_index, cases, deaths) %>%
  mutate(date = as.Date(date))

#### Without uncertainty (using estimated outcomes and then dividing, w/o. having to round up cases)
outcomes_aggregated <- estimate_outcomes(df_transformed, delay_density = function(x) density(onset_death, x))
outcomes_aggregated$cfr <- outcomes_aggregated$deaths / outcomes_aggregated$estimated_outcomes

#### With uncertainty (using cfr_rolling and cfr_time_varying, first having to round up cases)
# Rolling cfr
cfr_aggregated_rolling <- cfr_rolling(df_transformed, delay_density = function(x) density(onset_death, x))

# Time-varying cfr
cfr_aggregated_tv <- cfr_time_varying(df_transformed, delay_density = function(x) density(onset_death, x))

#4. Method 2: Averaging incidence across days of the week
weekly_avg <- df_transformed %>%
  mutate(
    week = as.integer((date - min(date)) / 7),  # Assigns a week number
    cases_weekly = ifelse(cases > 0, cases / 7, NA),  # Only divide where cases exist
    deaths_weekly = ifelse(deaths > 0, deaths / 7, NA)
  ) %>%
  group_by(week) %>%
  mutate(
    cases = ifelse(any(!is.na(cases_weekly)), cases_weekly, 0),  # If no cases in the week, stay 0
    deaths = ifelse(any(!is.na(deaths_weekly)), deaths_weekly, 0)
  ) %>%
  ungroup() %>%
  select(-cases_weekly, -deaths_weekly, -week)

#### Without uncertainty (using estimated outcomes and then dividing, w/o. having to round up cases)
outcomes_averaged <- estimate_outcomes(weekly_avg, delay_density = function(x) density(onset_death, x))
outcomes_averaged$cfr <- outcomes_averaged$deaths / outcomes_averaged$estimated_outcomes

#### With uncertainty (using cfr_rolling and cfr_time_varying, first having to round up cases)
#Rolling cfr
weekly_avg[,2:3] <- round(weekly_avg[,2:3],0) # Need to round for cfr
cfr_averaged_rolling <- cfr_rolling(weekly_avg, delay_density = function(x) density(onset_death, x))

#Time-varying cfr
cfr_averaged_tv <- cfr_time_varying(weekly_avg, delay_density = function(x) density(onset_death, x))

#5. Method 3: using EpiEstim to reconstruct daily incidence

#Getting si for COVID and formatting it to work with EpiEstim
COVID_si <- epiparameter_db(disease = "COVID", epi_name = "serial_interval", subset = is_parameterised, single_epiparameter = TRUE)
COVID_si <- discretise(COVID_si)
wrap_si <- function(si) {
  domain <- seq(1L, to = si$prob_dist$qf(0.999), by = 1L)
  pmf <- si$prob_dist$d(domain)
  pmf[1] <- 0
  pmf <- pmf / sum(pmf)
  pmf
}
si_distr <- wrap_si(COVID_si)

#Using estimate_R with dt=7 to reconstruct daily incidence
Rt_est_incidence <- estimate_R(incid = weekly_cases$cases,
                               dt = 7L,
                               dt_out = 7L,
                               recon_opt = "naive",
                               iter = 10L,
                               tol = 1e-6,
                               grid = list(precision = 0.001, min = -1, max = 1),
                               config = make_config(si_distr = si_distr),
                               method = "non_parametric_si")


I_inc <- as.data.frame(Rt_est_incidence$I[-c(1:6)])
names(I_inc)="cases"

#Now reconstructing also daily deaths
Rt_est_deaths <- estimate_R(incid = weekly_deaths$deaths,
                            dt = 7L,
                            dt_out = 7L,
                            recon_opt = "naive",
                            iter = 10L,
                            tol = 1e-6,
                            grid = list(precision = 0.001, min = -1, max = 1),
                            config = make_config(si_distr = si_distr),
                            method = "non_parametric_si")


I_deaths <- as.data.frame(Rt_est_deaths$I[-c(1:6)])
names(I_deaths)= "deaths"

#Now using reconstructed daily incidence and deaths to estimate CFR
#First we create the data frame for {cfr} with daily incidence and the corresponding dates

daily_rec_incidence <- bind_cols(daily_inc$date_index, I_inc)
names(daily_rec_incidence)[1]= "date"
daily_rec_deaths <- bind_cols(daily_deaths$date_index, I_deaths)
names(daily_rec_deaths)[1]= "date"

daily_rec_data <- merge(daily_rec_incidence, daily_rec_deaths, by = "date", all.x = T, all.y = T)
daily_rec_data <- replace(daily_rec_data, is.na(daily_rec_data), 0)

plot(daily_rec_data$date,daily_rec_data$cases,type="l", col = "blue")
lines(daily_rec_data$date,daily_rec_data$deaths, type = "l", col="red")

#### Without uncertainty (using estimated outcomes and then dividing, w/o. having to round up cases)
outcomes_rec <- estimate_outcomes(daily_rec_data, delay_density = function(x) density(onset_death, x))
outcomes_rec$cfr <- outcomes_rec$deaths / outcomes_rec$estimated_outcomes
outcomes_rec$date <- as.Date(outcomes_rec$date)

#### With uncertainty (using cfr_rolling and cfr_time_varying, first having to round up cases)
#Now, to be able to use cfr for daily incidence, it's necessary to round the incidence and deaths to an integer
daily_rec_data[,2:3] <- round(daily_rec_data[,2:3])
daily_rec_data$date <- as.Date(daily_rec_data$date)

#Rolling cfr
cfr_reconstructed_rolling <- cfr_rolling(daily_rec_data, delay_density = function(x) density(onset_death, x), poisson_threshold = 101)

#Time-varying cfr
cfr_reconstructed_tv <- cfr_time_varying(daily_rec_data, delay_density = function(x) density(onset_death, x))

#6. Actual cfr throughout the outbreak
linelist_inc <- incidence(linelist, date_index = c(cases = "date_onset", deaths = "date_death"), interval = 1) %>%
  complete_dates()
data_for_cfr <- prepare_data(linelist_inc, cases_variable = "cases", deaths_variable = "deaths")
data_for_cfr$date <- as.Date(data_for_cfr$date)

### Without uncertainty
outcomes_linelist <- estimate_outcomes(data_for_cfr, delay_density = function(x) density(onset_death, x))
outcomes_linelist$cfr <- outcomes_linelist$deaths/outcomes_linelist$estimated_outcomes

### With uncertainty
cfr_linelist_rolling <- cfr_rolling(data_for_cfr, delay_density = function(x) density(onset_death, x))
cfr_linelist_tv <- cfr_time_varying(data_for_cfr, delay_density = function(x) density(onset_death, x))

#### Plots with uncertainty #####

plot_rolling_linelist <- ggplot(cfr_linelist_rolling) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_tv_linelist <- ggplot(cfr_linelist_tv) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

#### Plots rolling cfr
plot_reconstructed_rolling <- ggplot(cfr_reconstructed_rolling) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_aggregated_rolling <- ggplot(cfr_aggregated_rolling) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_averaged_rolling <- ggplot(cfr_averaged_rolling) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

grid.arrange(plot_rolling_linelist, plot_aggregated_rolling, plot_averaged_rolling, plot_reconstructed_rolling, nrow = 1,
             top = "1- From linelist; 2- From aggregated on day 1; 3- From averaged cases; 4- From reconstructed incidence")

### Plots time-varying cfr
plot_reconstructed_tv <- ggplot(cfr_reconstructed_tv) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_aggregated_tv <- ggplot(cfr_aggregated_tv) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_averaged_tv <- ggplot(cfr_averaged_tv) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

grid.arrange(plot_tv_linelist, plot_aggregated_tv, plot_averaged_tv, plot_reconstructed_tv, nrow = 1,
             top = "1- From linelist; 2- From aggregated on day 1; 3- From averaged cases; 4- From reconstructed incidence")


#### Plots without uncertainty #####

plot_linelist <- ggplot(outcomes_linelist) +
  geom_line(
    aes(x = date, y = cfr), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_reconstructed <- ggplot(outcomes_rec) +
  geom_line(
    aes(x = date, y = cfr), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_aggregated <- ggplot(outcomes_aggregated) +
  geom_line(
    aes(x = date, y = cfr), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  scale_y_continuous(limits = c(0, 1)) + # Values >1
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_averaged <- ggplot(outcomes_averaged) +
  geom_line(
    aes(x = date, y = cfr), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

grid.arrange(plot_linelist, plot_aggregated, plot_averaged, plot_reconstructed, nrow = 1,
             top = "1- From linelist; 2- From aggregated on day 1; 3- From averaged cases; 4- From reconstructed incidence")

######### Option 2: Short delay with varying infection dynamics ##############

# 1. Simulating linelist
# create contact distribution
contact_distribution <- epiparameter::epiparameter(
  disease = "COVID-19",
  epi_name = "contact distribution",
  prob_distribution = create_prob_distribution(
    prob_distribution = "pois",
    prob_distribution_params = c(mean = 2)
  )
)

# create infectious period
infectious_period <- epiparameter::epiparameter(
  disease = "COVID-19",
  epi_name = "infectious period",
  prob_distribution = create_prob_distribution(
    prob_distribution = "gamma",
    prob_distribution_params = c(shape = 1, scale = 1)
  )
)

# create onset to death
onset_death_short <- epiparameter_db(
  disease = "Ebola",
  epi_name = "onset_to_death",
  single_epiparameter = TRUE)

# defining time-varying cfr
#config <- create_config(time_varying_death_risk = function(risk, time) risk * exp(-0.1 * time))

# simulating linelist
set.seed(1)
linelist_s <- sim_linelist(
  contact_distribution = contact_distribution,
  infectious_period = infectious_period,
  prob_infection = 0.5,
  onset_to_hosp = NULL,
  onset_to_death = onset_death_short,
  onset_to_recovery = NULL,
  non_hosp_death_risk = 0.2,
  outbreak_start_date = as.Date("2025-01-01"),
  outbreak_size = c(1000, 5000)

)

# 2. Converting to incidence and grouping by 7 days
# creating columns for deaths and recovery dates
linelist_s <- linelist_s %>%
  pivot_wider(
    names_from = outcome,
    values_from = date_outcome
  ) %>%
  rename(
    date_death = died,
  )

# converting to incidence
weekly_inc_s <- incidence(linelist_s,
                        date_index = c(cases = "date_onset", deaths = "date_death"),
                        interval = 7) %>%
  complete_dates()

# visualising cases, deaths and recoveries over time
plot(weekly_inc_s)

# 3. Method 1: CFR by assuming cases happen on day 1 of the week
# From weekly to daily incidence

weekly_cases_s <- weekly_inc_s %>%
  filter(!count == 0) %>%
  filter(count_variable == "cases") %>%
  mutate(date = str_extract(date_index, "^\\S+")) %>%
  mutate(cases = count) %>%
  select(date, cases)

daily_inc_s <- incidence(weekly_cases_s, date_index = "date", counts = "cases", interval = 1) %>%
  complete_dates()

weekly_deaths_s <- weekly_inc_s %>%
  filter(!count == 0) %>%
  filter(count_variable == "deaths") %>%
  mutate(date = str_extract(date_index, "^\\S+")) %>%
  mutate(deaths = count) %>%
  select(date, deaths)

daily_deaths_s <- incidence(weekly_deaths_s, date_index = "date", counts = "deaths", interval = 1) %>%
  complete_dates()

# merging both into a single dataset and formatting for CFR
daily_cases_deaths_s <- merge(daily_inc_s, daily_deaths_s, by = "date_index", all = TRUE)

df_transformed_s <- daily_cases_deaths_s %>%
  mutate(
    cases = ifelse(count_variable.x == "cases", count.x, 0),
    deaths = ifelse(count_variable.y == "deaths", count.y, 0)
  ) %>%
  mutate(
    cases = ifelse(is.na(cases), 0, cases),
    deaths = ifelse(is.na(deaths), 0, deaths)
  ) %>%
  select(date = date_index, cases, deaths) %>%
  mutate(date = as.Date(date))

#### Without uncertainty (using estimated outcomes and then dividing, w/o. having to round up cases)
outcomes_aggregated_s <- estimate_outcomes(df_transformed_s, delay_density = function(x) density(onset_death_short, x))
outcomes_aggregated_s$cfr <- outcomes_aggregated_s$deaths / outcomes_aggregated_s$estimated_outcomes

#### With uncertainty (using cfr_rolling and cfr_time_varying, first having to round up cases)
# Rolling cfr
cfr_aggregated_rolling_s <- cfr_rolling(df_transformed_s, delay_density = function(x) density(onset_death_short, x))

# Time-varying cfr
cfr_aggregated_tv_s <- cfr_time_varying(df_transformed_s, delay_density = function(x) density(onset_death_short, x))

#4. Method 2: Averaging incidence across days of the week
weekly_avg_s <- df_transformed_s %>%
  mutate(
    week = as.integer((date - min(date)) / 7),  # Assigns a week number
    cases_weekly = ifelse(cases > 0, cases / 7, NA),  # Only divide where cases exist
    deaths_weekly = ifelse(deaths > 0, deaths / 7, NA)
  ) %>%
  group_by(week) %>%
  mutate(
    cases = ifelse(any(!is.na(cases_weekly)), cases_weekly, 0),  # If no cases in the week, stay 0
    deaths = ifelse(any(!is.na(deaths_weekly)), deaths_weekly, 0)
  ) %>%
  ungroup() %>%
  select(-cases_weekly, -deaths_weekly, -week)

#### Without uncertainty (using estimated outcomes and then dividing, w/o. having to round up cases)
outcomes_averaged_s <- estimate_outcomes(weekly_avg_s, delay_density = function(x) density(onset_death_short, x))
outcomes_averaged_s$cfr <- outcomes_averaged_s$deaths / outcomes_averaged_s$estimated_outcomes

#### With uncertainty (using cfr_rolling and cfr_time_varying, first having to round up cases)
#Rolling cfr
weekly_avg_s[,2:3] <- round(weekly_avg_s[,2:3],0) # Need to round for cfr
cfr_averaged_rolling_s <- cfr_rolling(weekly_avg_s, delay_density = function(x) density(onset_death_short, x))

#Time-varying cfr
cfr_averaged_tv_s <- cfr_time_varying(weekly_avg_s, delay_density = function(x) density(onset_death_short, x))

#5. Method 3: using EpiEstim to reconstruct daily incidence

#Getting si for COVID and formatting it to work with EpiEstim
COVID_si <- epiparameter_db(disease = "COVID", epi_name = "serial_interval", subset = is_parameterised, single_epiparameter = TRUE)
COVID_si <- discretise(COVID_si)
wrap_si <- function(si) {
  domain <- seq(1L, to = si$prob_dist$qf(0.999), by = 1L)
  pmf <- si$prob_dist$d(domain)
  pmf[1] <- 0
  pmf <- pmf / sum(pmf)
  pmf
}
si_distr <- wrap_si(COVID_si)

#Using estimate_R with dt=7 to reconstruct daily incidence
Rt_est_incidence_s <- estimate_R(incid = weekly_cases_s$cases,
                               dt = 7L,
                               dt_out = 7L,
                               recon_opt = "naive",
                               iter = 10L,
                               tol = 1e-6,
                               grid = list(precision = 0.001, min = -1, max = 1),
                               config = make_config(si_distr = si_distr),
                               method = "non_parametric_si")


I_inc_s <- as.data.frame(Rt_est_incidence_s$I[-c(1:6)])
names(I_inc_s)="cases"

#Now reconstructing also daily deaths
Rt_est_deaths_s <- estimate_R(incid = weekly_deaths_s$deaths,
                            dt = 7L,
                            dt_out = 7L,
                            recon_opt = "naive",
                            iter = 10L,
                            tol = 1e-6,
                            grid = list(precision = 0.001, min = -1, max = 1),
                            config = make_config(si_distr = si_distr),
                            method = "non_parametric_si")


I_deaths_s <- as.data.frame(Rt_est_deaths_s$I[-c(1:6)])
names(I_deaths_s)= "deaths"

#Now using reconstructed daily incidence and deaths to estimate CFR
#First we create the data frame for {cfr} with daily incidence and the corresponding dates

daily_rec_incidence_s <- bind_cols(daily_inc_s$date_index, I_inc_s)
names(daily_rec_incidence_s)[1]= "date"
daily_rec_deaths_s <- bind_cols(daily_deaths_s$date_index, I_deaths_s)
names(daily_rec_deaths_s)[1]= "date"

daily_rec_data_s <- merge(daily_rec_incidence_s, daily_rec_deaths_s, by = "date", all.x = T, all.y = T)
daily_rec_data_s <- replace(daily_rec_data_s, is.na(daily_rec_data_s), 0)

plot(daily_rec_data_s$date,daily_rec_data_s$cases,type="l", col = "blue")
lines(daily_rec_data_s$date,daily_rec_data_s$deaths, type = "l", col="red")

#### Without uncertainty (using estimated outcomes and then dividing, w/o. having to round up cases)
outcomes_rec_s <- estimate_outcomes(daily_rec_data_s, delay_density = function(x) density(onset_death_short, x))
outcomes_rec_s$cfr <- outcomes_rec_s$deaths / outcomes_rec_s$estimated_outcomes
outcomes_rec_s$date <- as.Date(outcomes_rec_s$date)

#### With uncertainty (using cfr_rolling and cfr_time_varying, first having to round up cases)
#Now, to be able to use cfr for daily incidence, it's necessary to round the incidence and deaths to an integer
daily_rec_data_s[,2:3] <- round(daily_rec_data_s[,2:3])
daily_rec_data_s$date <- as.Date(daily_rec_data_s$date)

#Rolling cfr
cfr_reconstructed_rolling_s <- cfr_rolling(daily_rec_data_s, delay_density = function(x) density(onset_death_short, x), poisson_threshold = 101)

#Time-varying cfr
cfr_reconstructed_tv_s <- cfr_time_varying(daily_rec_data_s, delay_density = function(x) density(onset_death_short, x))

#6. Actual cfr throughout the outbreak
linelist_inc_s <- incidence(linelist_s, date_index = c(cases = "date_onset", deaths = "date_death"), interval = 1) %>%
  complete_dates()
data_for_cfr_s <- prepare_data(linelist_inc_s, cases_variable = "cases", deaths_variable = "deaths")
data_for_cfr_s$date <- as.Date(data_for_cfr_s$date)

### Without uncertainty
outcomes_linelist_s <- estimate_outcomes(data_for_cfr_s, delay_density = function(x) density(onset_death_short, x))
outcomes_linelist_s$cfr <- outcomes_linelist_s$deaths/outcomes_linelist_s$estimated_outcomes

### With uncertainty
cfr_linelist_rolling_s <- cfr_rolling(data_for_cfr_s[11:62,], delay_density = function(x) density(onset_death_short, x), poisson_threshold = 150)
cfr_linelist_tv_s <- cfr_time_varying(data_for_cfr_s, delay_density = function(x) density(onset_death_short, x))

#### Plots with uncertainty #####

plot_rolling_linelist_s <- ggplot(cfr_linelist_rolling_s) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_tv_linelist_s <- ggplot(cfr_linelist_tv_s) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

#### Plots rolling cfr
plot_reconstructed_rolling_s <- ggplot(cfr_reconstructed_rolling_s) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_aggregated_rolling_s <- ggplot(cfr_aggregated_rolling_s) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_averaged_rolling_s <- ggplot(cfr_averaged_rolling_s) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

grid.arrange(plot_rolling_linelist_s, plot_aggregated_rolling_s, plot_averaged_rolling_s, plot_reconstructed_rolling_s, nrow = 1,
             top = "1- From linelist; 2- From aggregated on day 1; 3- From averaged cases; 4- From reconstructed incidence")

### Plots time-varying cfr
plot_reconstructed_tv_s <- ggplot(cfr_reconstructed_tv_s) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_aggregated_tv_s <- ggplot(cfr_aggregated_tv_s) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_averaged_tv_s <- ggplot(cfr_averaged_tv_s) +
  geom_ribbon(
    aes(x = date, ymin = severity_low, ymax = severity_high),
    alpha = 0.5, fill = "deepskyblue3") +
  geom_line(
    aes(x = date, y = severity_estimate), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

grid.arrange(plot_tv_linelist_s, plot_aggregated_tv_s, plot_averaged_tv_s, plot_reconstructed_tv_s, nrow = 1,
             top = "1- From linelist; 2- From aggregated on day 1; 3- From averaged cases; 4- From reconstructed incidence")


#### Plots without uncertainty #####

plot_linelist_s <- ggplot(outcomes_linelist_s) +
  geom_line(
    aes(x = date, y = cfr), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_reconstructed_s <- ggplot(outcomes_rec_s) +
  geom_line(
    aes(x = date, y = cfr), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_aggregated_s <- ggplot(outcomes_aggregated_s) +
  geom_line(
    aes(x = date, y = cfr), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  scale_y_continuous(limits = c(0, 1)) + # Values >1
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

plot_averaged_s <- ggplot(outcomes_averaged_s) +
  geom_line(
    aes(x = date, y = cfr), colour = "royalblue4"
  ) +
  scale_x_date(date_labels = "%b-%Y") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Date", y = "CFR"
  ) + theme_bw()

grid.arrange(plot_linelist, plot_aggregated, plot_averaged, plot_reconstructed, nrow = 1,
             top = "1- From linelist; 2- From aggregated on day 1; 3- From averaged cases; 4- From reconstructed incidence")
