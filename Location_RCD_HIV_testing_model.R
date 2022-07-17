### Code to model HIV transmission amongst a network of PWID who are different locations. the goal is to investigate the difference between location based testing strategies and rcd strategies. 

# Load libraries 
library(numbers)
library(statnet)
library(ergm.ego)
library(ergm)
library(wesanderson)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)

## set up network, time steps 
### time units are monthly 
n <-500 ## population size of the entire population
n_years <-  10
time_steps <- n_years*12 ## monthly time steps for 5 years 

### demographic information
## 18 - 65
age_distribution <- 0.1*c(rep(2,(25-18)), rep(22,(35-25)), rep(16,(45-35)), rep(60,(66-45)))
age_distribution <- age_distribution/sum(age_distribution)
age_vec = sample(18:65, size = n, prob = age_distribution, replace = TRUE)
age_group_vec <- sapply(age_vec, function(x) cut(x, breaks = c(0,30,50,100)))

## number of locations individuals frequent 
n_locations = 10
location_distribution <- rep(1/n_locations, n_locations) ## probability an individual will be assigned to a specific location - NEED TO UPDATE

library(igraph)
library(randnet)
pl_fit_alpha = 3.431536 ## data from bridgers paper 
mean_ego_network_dist = 3.556054
sim_model <- BlockModel.Gen(lambda = mean_ego_network_dist, n = n, beta = 0.6, K = n_locations, alpha = pl_fit_alpha, Pi = location_distribution)
## lambda = avg degree, n = number of nodes, beta = prob within block edge, k = number of blocks 
# A - adj matrix; g = community membership; P = probability matrix of network; theta = ndoe degree parameter
inj_network <- sim_model$A
group_membership <- sim_model$g ## WILL NEED TO COUPLE THIS WILL location distribution 

## group membership from the network is the vector of assigned locations 
location_vec <- group_membership
## note: make stochastic block model network of PWID - we could just change this to use the actual network later 

### mortality rates - both age specific and overdose 
mortality_rates <- c(rep(0,18), rep(1.2, 25-18), rep(2.5, 35-25), rep(3.5, 45-35), rep(6.7, 55-45), rep(12.2, 66-44)) ## per 1,000 person-years -- need to update depending on age_distribution 
### needs a scaling for HIV positive ??? 
### overdose rates  ???? is this the most updated value we could use 
fixed_od_rate <- 13/1000/12#28.3/100000 # per 100,000 standard population https://www.cdc.gov/nchs/products/databriefs/db428.htm


pop_matrix <- matrix(,n,time_steps) ## keep track of each persons age
### for monthly data
n_years = time_steps/12
pop_matrix <- t(sapply(age_vec, function(x) seq(x, x+n_years, length = time_steps)))
pop_matrix <- t(apply(pop_matrix, 1, function(x) floor(x)))

###

prob_sexual <- 0.13 ### probability if they are an injection partner then they are also a sexual partner 
sexual_network <- t(apply(inj_network, 1, function(x) ifelse(x == 1, rbinom(1, size = 1, prob=prob_sexual),0)))

### 1 if they are a tie -- either injection or sexual contact 

### HIV testing
hiv_test_prob = 0.01*(20.4/12) ## https://academic.oup.com/cid/article/70/6/1096/5506442#205397035 prop tested per month in the past year 
art_update <- 0.9 ## once infected, likelihood of ART inititation [should this vary by month?]
achieve_viral_supress <- 1 ## made up ALIVE 
art_duration <- 0.70 ## probability per month of staying on ART -- not currently using 

### hiv transmission
prob_hiv_per_share <- 0.63 ## CID paper 
prob_hiv_per_sexual <- mean(c(0.11,0.04)) ## based on sexual act -- need to update

reduce_hiv_trans_w_viral_supress <- 0.94
reduce_hiv_trans_w_out_viral_spress <- 0.17

base_syring_share_prob = 0.011

### modeling starting/stopping patterns - mostly to set up the options and make a preliminary plot 

stop_list <- data.frame(year_vec = seq(0,10,by = 0.5), 
                        early_line = c(1.0, 0.95, 0.7, 0.52, 0.4, 0.33, 0.25, 0.18, 0.12, 0.1, 0.09, rep(0.05, 10)), 
                        delay_line = c(1.0, 1.0, 0.97, 0.96, 0.95, 0.92, 0.88, 0.83, 0.77, 0.69, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.17, 0.12, 0.1, 0.1, 0.09), 
                        persistent_line = c(seq(1.0, 0.95, length = 15), seq(0.94,0.9,length=21-15)))

start <- list(beta = 2)
exp_model_early <- nls(early_line ~ exp(beta * year_vec), data = stop_list, start = start)
exp_model_delay <- nls(delay_line ~ exp(beta * year_vec), data = stop_list, start = start)
exp_model_persistent <- nls(persistent_line ~ exp(beta * year_vec), data = stop_list, start = start)

exp_model_early_predict <- predict(exp_model_early, list(year_vec = seq(0,n_years*2, length = time_steps)) )
exp_model_delay_predict <- predict(exp_model_delay, list(year_vec = seq(0,n_years*2, length = time_steps)) )
exp_model_persistent_predict <- predict(exp_model_persistent, list(year_vec = seq(0,n_years*2, length = time_steps)))

## figuring out cessation/relapse patterns per individual

make_stop_start_matrix <- function(n_years, time_steps, exp_model_early, exp_model_delay, exp_model_persistent, n, age_vec){
  ### 1 = early cessation; 2 = delayed cessation; 3 = random; 4 = persistent 
  units = time_steps/n_years
  test_years <- seq(0,n_years*50, length = n_years*50*units)
  exp_model_early_predict <- predict(exp_model_early, list(year_vec = test_years))
  exp_model_delay_predict <- predict(exp_model_delay, list(year_vec = test_years))
  exp_model_persistent_predict <- predict(exp_model_persistent, list(year_vec = test_years))
  exp_model_random_predict_prob <- 0.5
  
  prob_values <- matrix(,length(exp_model_early_predict), 4)
  prob_values[,1] = exp_model_early_predict
  prob_values[,2] = exp_model_delay_predict
  prob_values[,3] = 0.5
  prob_values[,4] = exp_model_persistent_predict
  
  prob_membership <- 0.01*c(18.7, (15.7+17.5), 16.2, 31.9) ## distribution of 4 class membership
  active_status_prob = matrix(,n,time_steps)
  membership_vec <- sample(1:4, n, replace = TRUE, prob = prob_membership)
  ### assumes everyone is a new injector? 
  #  start_spot <- rpois(n, lambda = 14)
  start_spot = 1
  active_status_prob = sapply(1:n, function(x) prob_values[start_spot:(start_spot+time_steps-1),membership_vec[x]])
  return(list(active_status_prob = t(active_status_prob), membership_vec = membership_vec))
}

### YOU MUST MAKE THESE MATRICES TO RUN THE MAIN CODE
stop_start_param <- make_stop_start_matrix(n_years, time_steps, exp_model_early, exp_model_delay, exp_model_persistent, n, age_vec)

active_status_prob = stop_start_param$active_status_prob
stop_start_membership = stop_start_param$membership_vec

### main code
## ratio testing: proportion of testing being done from each strategy 

run_simulation_loc <- function(risk_group, skew_matrix, syring_share_on, moud_on, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec, location_vec, loc_testing_binary, rcd_testing_binary, random_testing_binary, ratio_testing){
  moud_effect <- c(0.6, 0.6, 0.3, 0.1) ## how much should injecting events decrease by if using MOUD
  likelihood_on_moud <- c(0.3, 0.3, 0.4, 0.5) ## what is the probability that a person in each stop/start group is on MOUD when it is available - from ALIVE/BECKY 
  indiv_moud_on <- sapply(stop_start_membership, function(x) rbinom(1,size = 1, prob = likelihood_on_moud[x])) ## is the individual on MOUD
  
  ### the simulation results get stored here 
  mortality_matrix <- matrix(NA,n,time_steps)### mortality status: is the person alive or not? 1 = alive, 0 = deceased
  active_matrix <- matrix(NA,n,time_steps) ### active user: is the person currently an active user 1 = active, 0 = inactive
  hiv_status_matrix <- matrix(0,n,time_steps) ### hiv status: what is the individuals active status 1 = positive, 0 = negative 
  hiv_supp_matrix <- matrix(0,n,time_steps) ### hiv suppression: if infected, are they virially supressed? 1 = yes, 0 = no
  hiv_diag_matrix <- matrix(0,n,time_steps) ### hiv diag: has the individual been diagnoised with HIV or not  1 = yes, 0 = no
  
  ## where diagnoised -> 1 = location, 2 = random, 3 = rcd  
  diag_loc <- rep(NA, n) ## only for your initial diagnosis 
  
  inj_events_per_indiv <- matrix(NA,n,time_steps) ### number of injecting events per individual 
  mean_inj_events_per_partner <- matrix(NA,n,time_steps) ## mean number of injectinv events per individual
  
  ### start with everyone alive
  mortality_matrix[,1] = 1
  
  ### start with everyone being an active injector
  active_matrix[,1] = 1 
  
  ### random people to start as infected
  start_infected <- sample(1:n, 10)#ceiling(0.01*n))
  hiv_status_matrix[,1] = 0
  hiv_status_matrix[start_infected,1] = 1
  
  ## who infects whom -> 0 means no one infected you (you were a seed)
  wiw <- matrix(NA,n,n)
  for(bb in 1:length(start_infected)){
    wiw[start_infected[bb], start_infected[bb]] = 0
  }

  ### start with no one viral supp
  hiv_supp_matrix[start_infected,1] = 0
  ## simulate starting with time step 2
  for(ii in 2:time_steps){
    ## ii is the index for the time step 
    prev_step = ii-1 
    ## update HIV transmission 
    prev_hiv_status <- hiv_status_matrix[,prev_step]
    prev_infected <- which(prev_hiv_status == 1)
    prev_uninfected <- which(prev_hiv_status == 0)
    
    prev_hiv_supp <- hiv_supp_matrix[,prev_step]
    prev_suppressed <- which(prev_hiv_supp == 1)
    prev_unsuppress <- which(prev_hiv_supp == 0)
    if(length(prev_infected) > 0 ){ ## runs code by infected person 
      for(jj in 1:length(prev_infected)){
        ## jj is the index for the person 
        infected_index = prev_infected[jj]
        active_user <- ifelse(active_matrix[infected_index, prev_step] == 1, TRUE, FALSE) ## are an active user
        
        alive_person <- ifelse(mortality_matrix[infected_index, prev_step] == 1, TRUE, FALSE) ## person is alive 
        
        if(active_user == TRUE & alive_person == TRUE){
          early_infect <- ifelse(sum(hiv_status_matrix[infected_index,prev_step:ii], na.rm=T) == 1, 0, 0) ## new infection -- scale higher viral load  -- this is not correct and currently not being used 
          viral_supress = hiv_supp_matrix[jj,prev_step] ## were you previously virally supressed 
          
          viral_supp_scale <- 0.94 *viral_supress ## need to change this but it will scale viral supression 
          
          # reduce_hiv_trans_w_viral_supress <- 0.94
          # reduce_hiv_trans_w_out_viral_spress <- 0.17       
          
          sex_partners <- which(sexual_network[infected_index] == 1)
          
          ## infected injection partners
          inj_partners <- which(inj_network[infected_index,] == 1) ## index of injecting partners 
          inj_partners <- inj_partners[which(active_matrix[inj_partners,prev_step] == 1)] ## are their injecting partners active 
          ### change this to only look at HIV negative partners 
          
          n_inj_partners <- length(inj_partners)
          n_sex_partners <- length(sex_partners)
          
          inj_freq <- rpois(1, lambda_vec[risk_group[jj,ii]]) ## total number of injecting acts per month
          inj_events_per_indiv[jj,ii] = inj_freq
          if(n_sex_partners>0){
            sex_freq <- rpois(n_sex_partners, lambda = 5) ## this will need to be updated to account for heterosexual vs homosexual sexual events 
            
            prob_hiv_per_sex_w_suppress <- pmax(prob_hiv_per_sexual  - (viral_supp_scale*prob_hiv_per_sexual), 0)
            
            prob_hiv_trans_sex_partner <- pmin(prob_hiv_per_sex_w_suppress, 1)
            prob_hiv_trans_sex_partner <- pmin(early_infect*prob_hiv_trans_sex_partner + prob_hiv_trans_sex_partner, 1)
            prob_hiv_trans_sex_partner <- sapply(1:n_sex_partners, function(x) rbinom(1, size = sex_freq[x], prob = prob_hiv_trans_sex_partner))
            
            hiv_status_matrix[sex_partners,ii] = prob_hiv_trans_sex_partner
            
            wiw[sex_partners[which(prob_hiv_trans_sex_partner == 1)],infected_index] = 1
          }
          
          
          if(n_inj_partners > 0 ){
            ## number of injection events per individuals 
            skew_opt = skew_type[skew_matrix[jj,ii],]
            prob_inj_numb_by_partner <- rbeta(n_inj_partners, shape1 = skew_opt[1], shape2 = skew_opt[2])
            inj_acts_per_partner <- sample(1:n_inj_partners, size = inj_freq, replace = TRUE, prob = prob_inj_numb_by_partner/sum(prob_inj_numb_by_partner, na.rm=T)) 
            inj_acts_per_partner <- sapply(1:n_inj_partners, function(x) length(which(inj_acts_per_partner == x)))
            mean_inj_events_per_partner[jj,ii] = mean(inj_acts_per_partner, na.rm=T)
            
            syringe_service_effect <- syring_share_on[ii] ## will reduce the probability of sharing a syringe
            # moud_effect_share <- 0.1 * moud_on[ii] ## will reduce the probability of sharing per month -- will be 0 if moud is not available 
            
            moud_effect_inj <- moud_on[ii] * indiv_moud_on[jj] * moud_effect[stop_start_membership[jj]]  ## will reduce inj_act_per_month -- will be 0 if moud is not available
            
            inj_acts_per_partner <- floor(pmax((inj_acts_per_partner - moud_effect_inj*inj_acts_per_partner), 0)) ## set to be at least 0 injecting acts per partner
            prob_syringe_share <- pmax(base_syring_share_prob - (base_syring_share_prob*syringe_service_effect),0) ## reduce syringe sharing based on syringe services being available 
            share_syringe_events_per_partner <- sapply(inj_acts_per_partner, function(x) sum(rbinom(x, size = 1, prob = prob_syringe_share))) ## identify how many syringe shaing events there are per partner -- a single probability per person distributed based on the number of injecting events by partner 
            prob_hiv_per_share_w_suppress <- pmax(prob_hiv_per_share  - (viral_supp_scale*prob_hiv_per_share), 0) ## rescale the probability of HIV transmisison given viral suppression 
            prob_hiv_trans_inj_partner <- sapply(1:n_inj_partners, function(x) rbinom(1, size = share_syringe_events_per_partner[x], prob = prob_hiv_per_share_w_suppress)) ## using a binomial distribution, determine how many contacts become infected with HIV 
            hiv_status_matrix[inj_partners,ii] = prob_hiv_trans_inj_partner ## update the matrix on HIV status 
            wiw[sex_partners[which(inj_partners == 1)],infected_index] = 1
          }
        }
      }
    }
    hiv_status_matrix[prev_infected,ii] = 1 ### anyone who was previously infected stays infected 
    
    ### testing - you either will be randomly tested, RCD tested or location based testing
    
    
    baseline_testing <- rbinom(n, size = 1, prob = hiv_test_prob) ## right now there is no contact tracing increase in likelihood of being tested 
    new_diagnoised <- which(hiv_status_matrix[,ii] == 1 & baseline_testing == 1)
    hiv_diag_matrix[new_diagnoised,ii] = 1
    
    diagnoised_hiv <- which(hiv_diag_matrix[,ii] == 1) 
    n_diagnoised <- length(diagnoised_hiv)
    if(n_diagnoised > 0){ ## there are people who are diagnoised
      on_art_and_supress <- rbinom(n_diagnoised, size = 1, prob = art_update*achieve_viral_supress) 
      hiv_supp_matrix[diagnoised_hiv[which(on_art_and_supress == 1)],ii] = 1
      hiv_supp_matrix[diagnoised_hiv[which(on_art_and_supress == 0)], ii] = 0
    }
    
    ## update demography -- people stopping, dying, start
    pop_status <- pop_matrix[,prev_step] 
    prev_mortality <- mortality_matrix[,prev_step]
    mortality_vec <- sapply(pop_status, function(x) rbinom(1, size = 1, prob = 1-(fixed_od_rate+mortality_rates[x]/100/12))) ## this is wrong , if this is a 0 then the person has died 
    mortality_matrix[,ii] = mortality_vec
    mortality_matrix[which(prev_mortality == 0),ii] = 0 ## if you were previously dead you remain dead 
    active_status <- active_matrix[,prev_step]
    
    update_status = sapply(active_status_prob[,ii], function(x) rbinom(1, size = 1, prob = x))
    active_matrix[,ii] = update_status
  }
  return(list(mortality_matrix = mortality_matrix, active_matrix = active_matrix, hiv_status_matrix = hiv_status_matrix, hiv_supp_matrix = hiv_supp_matrix, hiv_diag_matrix = hiv_diag_matrix, inj_events_per_indiv = inj_events_per_indiv, mean_inj_events_per_partner = mean_inj_events_per_partner))
  
}

# median of 14 years of injecting at baseline (iqr = 7-20)

## make skewed distribution dependent upon degree 
## fit bridger data 

## stopping/starting behavior 

### individuals can have a risk category, right now it is set to be equal that can be changed that will determine their injecting behavior 
risk_cat <- sample(1:4, size = n, prob = c(0.25, 0.25, 0.25, 0.25), replace = TRUE) ## 1 = standard, 2 = higher risk, 3 = lower risk, 4 = homeless 
risk_group <- matrix(risk_cat,n,time_steps)

### how do we want injecting events distributed amongst contacts
skew_0 <- c(2,2) ## alpha, beta values for beta distribution 
skew_1 <- c(1,5)
skew_2 <- c(1,7)
skew_type <- matrix(,3,2)
skew_type[1,] = skew_0
skew_type[2,] = skew_1
skew_type[3,] = skew_2

## 
### make a plot to demonstrate distribution of injecting events/skew in contact+injecting events 
lambda_vec <- c(30, 60, 100, 150)
## 1 = lowest, 2 = moderate, 3 = high, 4 = highest

par(mfrow=c(1,1))
pal_test <- brewer.pal(9, 'Paired')
pal_test <- c(pal_test[2], pal_test[4], pal_test[8])
plot(NA, NA, xlim = c(0,1), ylim = range(lambda_vec))
points(c(0.91, 0.03, 0.03, 0.03), lambda_vec, col = pal_test, pch = 16)
lines(c(0.91, 0.03, 0.03, 0.03), lambda_vec, col = pal_test, pch = 16)
points(c(0.25, 0.25, 0.25, 0.25), lambda_vec, col = pal_test, pch = 16)
points(c(0.03, 0.03, 0.03, 0.91), lambda_vec, col = pal_test, pch = 16)

test_matrix <- matrix(,3,4)
test_matrix[1,] = c(0.25, 0.25, 0.25, 0.25)
test_matrix[2,] = c(0.91, 0.03, 0.03, 0.03)
test_matrix[3,] = c(0.03, 0.03, 0.03, 0.91)
colnames(test_matrix) <- lambda_vec
rownames(test_matrix) <- c('base', 'less', 'more')

barplot(test_matrix, beside = TRUE, col = pal_test, xlab = 'Avg inj events/month', ylab = 'Prob membership')
legend('top', legend = c('Base', 'Less inj', 'More inj'), col = pal_test, pch = 15, bty = 'n', ncol = 1)

### unorganized code to run some simulations 

## base scenario
lambda_vec <- c(30, 60, 100, 150)
base_risk_cat <- sample(1:4, size = n, prob = c(0.25, 0.25, 0.25, 0.25), replace = TRUE) 
base_risk_group <- matrix(base_risk_cat,n,time_steps)
base_skew_matrix<-matrix(1,n,time_steps)

### same skew, people are injecting less frequently 
less_freq_inj_risk_cat <- sample(1:4, size = n, prob = c(0.91, 0.03, 0.03, 0.03), replace = TRUE) 
less_freq_inj_risk_group <- matrix(less_freq_inj_risk_cat,n,time_steps)

### same skew, people are injecting more frequently
more_freq_inj_risk_cat <- sample(1:4, size = n, prob = c(0.03, 0.03, 0.03, 0.91), replace = TRUE) 
more_freq_inj_risk_group <- matrix(more_freq_inj_risk_cat,n,time_steps)

## even mix between 1 and 2
skew_matrix_2 <- matrix(sample(1:2,n, replace = TRUE),n,time_steps)
# 
# syring_share_on_2 <- c(rep(0.5, 2*12), rep(0, 12), rep(0.5, ((time_steps/12)-2-1)*12)) ## index if syring sharing services are available 
# moud_on_2 <- syring_share_on_2 ## index if MOUD services are available

syring_share_on <- rep(0.75, time_steps) #c(rep(0, time_steps/2), rep(0, time_steps/2)) ## index if syring sharing services are available 
moud_on <- rep(1, time_steps) #c(rep(0, time_steps/2), rep(0, time_steps/2)) ## index if MOUD services are available 1 if it is available, 0 if not 

syring_share_on_3m <- syring_share_on
syring_share_on_3m[25:(25+2)] = 0
moud_on_3m <- moud_on
moud_on_3m[25:(25+2)] = 0

syring_share_on_6m <- syring_share_on
syring_share_on_6m[25:(25+5)] = 0
moud_on_6m <- moud_on
moud_on_6m[25:(25+5)] = 0

syring_share_on_12m <- syring_share_on
syring_share_on_12m[25:(25+11)] = 0
moud_on_12m <- moud_on
moud_on_12m[25:(25+11)] = 0


n_sim = 25
### xx
total_numb_w_base <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  sim_base_results <- run_simulation(base_risk_group, base_skew_matrix, syring_share_on, moud_on, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_base[ii,] = colSums(sim_base_results$hiv_status_matrix)
}

### simulations with covid-19 disruptions lasting 3, 6 and 12 months - both MOUD and SSP
total_numb_w_covid_3m <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  sim_base_results <- run_simulation(base_risk_group, base_skew_matrix, syring_share_on_3m, moud_on_3m, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_covid_3m[ii,] = colSums(sim_base_results$hiv_status_matrix)
}

total_numb_w_covid_6m <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  sim_base_results <- run_simulation(base_risk_group, base_skew_matrix, syring_share_on_6m, moud_on_6m, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_covid_6m[ii,] = colSums(sim_base_results$hiv_status_matrix)
}

total_numb_w_covid_12m <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  sim_base_results <- run_simulation(base_risk_group, base_skew_matrix, syring_share_on_12m, moud_on_12m, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_covid_12m[ii,] = colSums(sim_base_results$hiv_status_matrix)
}

### simulations with only changes to SSP

total_numb_w_covid_12m_only_ssp <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  sim_base_results <- run_simulation(base_risk_group, base_skew_matrix, syring_share_on_12m, moud_on, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_covid_12m_only_ssp[ii,] = colSums(sim_base_results$hiv_status_matrix)
}

### simulation with only change to MOUD

total_numb_w_covid_12m_only_moud <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  sim_base_results <- run_simulation(base_risk_group, base_skew_matrix, syring_share_on, moud_on_12m, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_covid_12m_only_moud[ii,] = colSums(sim_base_results$hiv_status_matrix)
}

### figure looking at just COVID-19 disruptions
t1 = 2*12
t2 = 7*12-1
mean_numb_w_base <- colMeans(total_numb_w_base[,t1:t2])
mean_numb_w_covid_3m <- colMeans(total_numb_w_covid_3m[,t1:t2])
mean_numb_w_covid_6m <- colMeans(total_numb_w_covid_6m[,t1:t2])
mean_numb_w_covid_12m <- colMeans(total_numb_w_covid_12m[,t1:t2])
mean_numb_w_covid_12m_only_ssp <- colMeans(total_numb_w_covid_12m_only_ssp[,t1:t2])
mean_numb_w_covid_12m_only_moud <- colMeans(total_numb_w_covid_12m_only_moud[,t1:t2])

library(RColorBrewer)
pal_1 <- brewer.pal(9, 'Set1')
plot(NA, NA, xlim = c(1,60), ylim = c(-100,200), xlab = 'Months', ylab = '% Change over Baseline')
abline(h=0, lty = 2, col = 'red')
lines(100*(mean_numb_w_covid_3m - mean_numb_w_base)/mean_numb_w_base, type = 'b', col = pal_1[2], pch = 16)
lines(100*(mean_numb_w_covid_6m - mean_numb_w_base)/mean_numb_w_base, type = 'b', col = pal_1[3], pch = 16)
lines(100*(mean_numb_w_covid_12m - mean_numb_w_base)/mean_numb_w_base, type = 'b', col = pal_1[4], pch = 16)
lines(100*(mean_numb_w_covid_12m_only_ssp - mean_numb_w_covid_12m)/mean_numb_w_covid_12m, type = 'b', col = pal_1[4], pch = 15)
lines(100*(mean_numb_w_covid_12m_only_moud - mean_numb_w_covid_12m)/mean_numb_w_covid_12m, type = 'b', col = pal_1[4], pch = 14) ## something wrong with this one
legend('topright', bty = 'n', legend = c('3m', '6m', '12m', '12m+SSP', '12m+MOUD'), col = c(pal_1[2], pal_1[3], pal_1[4], pal_1[4], pal_1[4]), pch = c(16, 16, 16, 15, 14), ncol = 2)

### simulations with changes in distribution (skew) in injecting events 

total_numb_w_skew <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  sim_base_results <- run_simulation(base_risk_group, skew_matrix_2, syring_share_on, moud_on, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_skew[ii,] = colSums(sim_base_results$hiv_status_matrix)
}


total_numb_w_less_freq <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  sim_less_freq_inj_results <- run_simulation(less_freq_inj_risk_group, base_skew_matrix, syring_share_on, moud_on, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_less_freq[ii,] = colSums(sim_less_freq_inj_results$hiv_status_matrix)
}

total_numb_w_less_freq_skew2 <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  sim_less_freq_inj_results <- run_simulation(less_freq_inj_risk_group, skew_matrix_2, syring_share_on, moud_on, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_less_freq_skew2[ii,] = colSums(sim_less_freq_inj_results$hiv_status_matrix)
}

total_numb_w_more_freq <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  ## if 1 - no sharing is possible 
  sim_less_freq_inj_results_2 <- run_simulation(more_freq_inj_risk_group, base_skew_matrix, syring_share_on, moud_on, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_more_freq[ii,] = colSums(sim_less_freq_inj_results_2$hiv_status_matrix)
}

total_numb_w_more_freq_skew2 <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  ## if 1 - no sharing is possible 
  sim_less_freq_inj_results_2 <- run_simulation(more_freq_inj_risk_group, skew_matrix_2, syring_share_on, moud_on, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_more_freq_skew2[ii,] = colSums(sim_less_freq_inj_results_2$hiv_status_matrix)
}

### figure looking at just changes in distribution of injecting events 
t1 = 2*12
t2 = 7*12-1
mean_numb_w_skew <- colMeans(total_numb_w_skew[,t1:t2])
mean_numb_w_less_freq <- colMeans(total_numb_w_less_freq[,t1:t2])
mean_numb_w_less_freq_skew2 <- colMeans(total_numb_w_less_freq_skew2[,t1:t2])
mean_numb_w_more_freq  <- colMeans(total_numb_w_more_freq[,t1:t2])
mean_numb_w_more_freq_skew2 <- colMeans(total_numb_w_more_freq_skew2[,t1:t2])

# library(RColorBrewer)
# pal_1 <- brewer.pal(9, 'Set1')
plot(NA, NA, xlim = c(1,60), ylim = c(-100,200), xlab = 'Months', ylab = '% Change over Baseline')
abline(h=0, lty = 2, col = 'red')
lines(100*(mean_numb_w_skew - mean_numb_w_base)/mean_numb_w_base, type = 'b', col = pal_1[5], pch = 16)
lines(100*(mean_numb_w_less_freq - mean_numb_w_base)/mean_numb_w_base, type = 'b', col = pal_1[6], pch = 16)
lines(100*(mean_numb_w_less_freq_skew2 - mean_numb_w_base)/mean_numb_w_base, type = 'b', col = pal_1[7], pch = 16)
lines(100*(mean_numb_w_more_freq - mean_numb_w_base)/mean_numb_w_base, type = 'b', col = pal_1[8], pch = 16)
lines(100*(mean_numb_w_more_freq_skew2 - mean_numb_w_base)/mean_numb_w_base, type = 'b', col = pal_1[9], pch = 16)

legend('topright', bty = 'n', legend = c('skew', 'less-freq', 'less-freq + skew', 'more-freq', 'more-freq + skew'), col = c(pal_1[5], pal_1[6], pal_1[7], pal_1[8], pal_1[9]), pch = 16, ncol = 2) ## skew is in distribution of injecting events




plot(NA, NA, xlim = c(1,60), ylim = c(-100,200), xlab = 'Months', ylab = '% Change over Baseline')
abline(h=0, lty = 2, col = 'red')
lines(100*(mean_numb_w_skew - mean_numb_w_base)/mean_numb_w_base,  col = pal_1[5], pch = 16, lwd = 2)
lines(100*(mean_numb_w_less_freq - mean_numb_w_base)/mean_numb_w_base, col = pal_1[6], pch = 16, lwd = 2)
lines(100*(mean_numb_w_more_freq - mean_numb_w_base)/mean_numb_w_base, col = pal_1[8], pch = 16, lwd = 2)
lines(100*(mean_numb_w_covid_3m - mean_numb_w_base)/mean_numb_w_base,  col = pal_1[2], pch = 16, lwd = 2)
lines(100*(mean_numb_w_covid_6m - mean_numb_w_base)/mean_numb_w_base, col = pal_1[3], pch = 16, lwd = 2)
lines(100*(mean_numb_w_covid_12m - mean_numb_w_base)/mean_numb_w_base,  col = pal_1[4], pch = 16, lwd = 2)
legend('topright', bty = 'n', legend = c('3m', '6m', '12m', '12m+SSP', '12m+MOUD'), col = c(pal_1[2], pal_1[3], pal_1[4], pal_1[4], pal_1[4]), pch = c(16, 16, 16, 15, 14), ncol = 2)
legend('topleft', bty = 'n', legend = c('skew', 'less-freq', 'less-freq + skew', 'more-freq', 'more-freq + skew'), col = c(pal_1[5], pal_1[6], pal_1[7], pal_1[8], pal_1[9]), pch = 16, ncol = 2) ## skew is in distribution of injecting events


### simulations with changes in distribution (skew) in injecting events and covid-19 disruptions 

total_numb_w_skew <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  sim_base_results <- run_simulation(base_risk_group, skew_matrix_2, syring_share_on, moud_on, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_skew[ii,] = colSums(sim_base_results$hiv_status_matrix)
}

total_numb_w_covid <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  sim_base_results <- run_simulation(base_risk_group, base_skew_matrix, syring_share_on_2, moud_on_2, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_covid[ii,] = colSums(sim_base_results$hiv_status_matrix)
}

total_numb_w_less_freq <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  sim_less_freq_inj_results <- run_simulation(less_freq_inj_risk_group, base_skew_matrix, syring_share_on, moud_on, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_less_freq[ii,] = colSums(sim_less_freq_inj_results$hiv_status_matrix)
}


total_numb_w_less_freq_off <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  ## if 1 - no sharing is possible 
  sim_less_freq_inj_results_2 <- run_simulation(less_freq_inj_risk_group, base_skew_matrix, syring_share_on_2, moud_on_2, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_less_freq_off[ii,] = colSums(sim_less_freq_inj_results_2$hiv_status_matrix)
}

boxplot(total_numb_w_less_freq_off)

total_numb_w_less_freq_skew2 <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  sim_less_freq_inj_results <- run_simulation(less_freq_inj_risk_group, skew_matrix_2, syring_share_on, moud_on, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_less_freq_skew2[ii,] = colSums(sim_less_freq_inj_results$hiv_status_matrix)
}

total_numb_w_less_freq_off_skew2 <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  ## if 1 - no sharing is possible 
  sim_less_freq_inj_results_2 <- run_simulation(less_freq_inj_risk_group, skew_matrix_2, syring_share_on_2, moud_on_2, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_less_freq_off_skew2[ii,] = colSums(sim_less_freq_inj_results_2$hiv_status_matrix)
}

total_numb_w_more_freq <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  ## if 1 - no sharing is possible 
  sim_less_freq_inj_results_2 <- run_simulation(more_freq_inj_risk_group, base_skew_matrix, syring_share_on, moud_on, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_more_freq[ii,] = colSums(sim_less_freq_inj_results_2$hiv_status_matrix)
}


total_numb_w_more_freq_off <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  ## if 1 - no sharing is possible 
  sim_less_freq_inj_results_2 <- run_simulation(more_freq_inj_risk_group, base_skew_matrix, syring_share_on_2, moud_on_2, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_more_freq_off[ii,] = colSums(sim_less_freq_inj_results_2$hiv_status_matrix)
}

total_numb_w_more_freq_skew2 <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  ## if 1 - no sharing is possible 
  sim_less_freq_inj_results_2 <- run_simulation(more_freq_inj_risk_group, skew_matrix_2, syring_share_on, moud_on, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_more_freq_skew2[ii,] = colSums(sim_less_freq_inj_results_2$hiv_status_matrix)
}

total_numb_w_more_freq_off_skew2 <- matrix(,n_sim,time_steps,)
for(ii in 1:n_sim){
  print(ii)
  ## if 1 - no sharing is possible 
  sim_less_freq_inj_results_2 <- run_simulation(more_freq_inj_risk_group, skew_matrix_2, syring_share_on_2, moud_on_2, initiation_rate, cessation_rate, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec)
  total_numb_w_more_freq_off_skew2[ii,] = colSums(sim_less_freq_inj_results_2$hiv_status_matrix)
}

values <- data.frame(total = c(total_numb_w_base[,time_steps], total_numb_w_less_freq[,time_steps], total_numb_w_less_freq_off[,time_steps], total_numb_w_more_freq[,time_steps], total_numb_w_more_freq_off[,time_steps],total_numb_w_less_freq_skew2[,time_steps], total_numb_w_less_freq_off_skew2[,time_steps], total_numb_w_more_freq_skew2[,time_steps], total_numb_w_more_freq_off_skew2[,time_steps]), 
                     type = c(rep('base',n_sim), rep('less freq',n_sim), rep('less freq + covid',n_sim), rep('more freq', n_sim), rep('more freq + covid', n_sim), rep('less freq + skew', n_sim), rep('less freq + covid + skew', n_sim), rep('more freq + skew', n_sim), rep('more freq + covid + skew', n_sim)))

par(mfrow=c(1,3))
boxplot(total_numb_w_base, ylab = 'HIV cases', xlab = 'Months', col = pal_test[1], main = 'Base')
boxplot(total_numb_w_less_freq, ylab = 'HIV cases', xlab = 'Months', col = pal_test[2], main = 'Less Inj')
boxplot(total_numb_w_more_freq, ylab = 'HIV cases', xlab = 'Month', col = pal_test[3], main = 'More Inj')

par(mfrow=c(1,3))
boxplot(total_numb_w_base, ylab = 'HIV cases', xlab = 'Months', col = pal_test[1], main = 'Base')
boxplot(total_numb_w_covid, ylab = 'HIV cases', xlab = 'Months', main = 'COVID', col = 'red')
boxplot(total_numb_w_less_freq_off, ylab = 'HIV cases', xlab = 'Months', col = pal_test[2], main = 'Less Inj + COVID')


par(mfrow=c(1,3))
boxplot(total_numb_w_base, ylab = 'HIV cases', xlab = 'Months', col = pal_test[1], main = 'Base')
boxplot(total_numb_w_covid, ylab = 'HIV cases', xlab = 'Months', main = 'COVID', col = 'red')
boxplot(total_numb_w_skew, ylab = 'HIV cases', xlab = 'Months', main = 'Skew Dist', col = 'grey')

time_stop = 12*2
avg_values <- data.frame(total = c(colMeans(total_numb_w_base[,1:time_stop]), colMeans(total_numb_w_less_freq[,1:time_stop]), colMeans(total_numb_w_less_freq_off[,1:time_stop]), colMeans(total_numb_w_more_freq[,1:time_stop]), colMeans(total_numb_w_more_freq_off[,1:time_stop]), colMeans(total_numb_w_less_freq_skew2[,1:time_stop]), colMeans(total_numb_w_less_freq_off_skew2[,1:time_stop]), colMeans(total_numb_w_more_freq_skew2[,1:time_stop]), colMeans(total_numb_w_more_freq_off_skew2[,1:time_stop]), colMeans(total_numb_w_covid[,1:time_stop])), 
                         type = c(rep('base',time_stop), rep('less freq',time_stop), rep('less freq + covid',time_stop), rep('more freq', time_stop), rep('more freq + covid', time_stop), rep('less freq + skew', time_stop), rep('less freq + covid + skew', time_stop), rep('more freq + skew', time_stop), rep('more freq + covid + skew', time_stop), rep('covid', time_stop)))

base_values <- filter(avg_values, type == 'base')
base_values_vec <- rep(base_values$total, times = length(unique(avg_values$type)))
avg_values$norm <- 100*(avg_values$total - base_values_vec)/base_values_vec


