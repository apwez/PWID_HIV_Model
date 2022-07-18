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
n_years <-  5
time_steps <- n_years*12 ## monthly time steps for 5 years 

### demographic information
## 18 - 65
## from data - median = 31; IQR = 26-37 --> sd = 1.36*IQR
library(truncnorm)
max_age = 50; min_age = 18
age_vec <- rtruncnorm(n, a = min_age, b = max_age, mean = 31, sd = 11/1.35)
age_group_vec <- sapply(age_vec, function(x) cut(x, breaks = c(0,25,35,45,55)))

## number of locations individuals frequent 
n_locations = 7 
location_distribution <- rep(1/n_locations, n_locations) ## probability an individual will be assigned to a specific location - NEED TO UPDATE

### should we replace this with the actual network? 
library(igraph)
library(randnet)
pl_fit_alpha = 3.431536 ## data from bridgers paper 
mean_ego_network_dist = 3.556054
sim_model <- BlockModel.Gen(lambda = mean_ego_network_dist, n = n, beta = 0.6, K = n_locations, alpha = pl_fit_alpha, Pi = location_distribution)
## lambda = avg degree, n = number of nodes, beta = prob within block edge, k = number of blocks 
# A - adj matrix; g = community membership; P = probability matrix of network; theta = ndoe degree parameter
inj_network <- sim_model$A
group_membership <- sim_model$g 

## group membership from the network is the vector of assigned locations 
location_vec <- group_membership
## note: make stochastic block model network of PWID - we could just change this to use the actual network later 

### mortality rates - both age specific and overdose 
fixed_od_rate = 1.2/100/12 ## monthly OD rate 
mortality_rates <- fixed_od_rate ### ignoring other types of mortality right now 

pop_matrix <- matrix(,n,time_steps) ## keep track of each persons age - need to figure out birth rate 
### for monthly data
n_years = time_steps/12
pop_matrix <- t(sapply(age_vec, function(x) seq(x, x+n_years, length = time_steps)))
pop_matrix <- t(apply(pop_matrix, 1, function(x) floor(x)))

###

prob_sexual <- 0.001## from Brown paper this is set to 0.13 -- making it 0.01 since its all inj transmission 
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

### more complicated code that includes risk categories

### main code
## risk_group: 4 categories of risk behavior 
## skew_matrix: how do you distribute your injecting events
## syring_share_on: binary if syringe services are available
## moud_on: binary if moud resources are available
## inj_network: network of injecting partners
## sexual_network: network of sexual contacts
## active_status_prob: probability that you are an active user based on risk category
## stop-start_membership: cessation/relapse patterns
## lambda_vec: mean injecting events/month
## location_vec: locations associated with each individual
## testing_loc: which locations have testing available
## rcd_prob: probability you will be a seed for a reactive case detection event 
run_simulation_loc <- function(risk_group, skew_matrix, syring_share_on, moud_on, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec, location_vec, testing_loc, rcd_prob){
  
  testing_loc = seq(1,length(testing_loc))[which(testing_loc == 1)] ## make a list of the ids of locations with testing 
  
  moud_effect <- c(0.6, 0.6, 0.3, 0.1) ## how much should injecting events decrease by if using MOUD
  likelihood_on_moud <- c(0.3, 0.3, 0.4, 0.5) ## what is the probability that a person in each stop/start group is on MOUD when it is available - from ALIVE/BECKY 
  
  indiv_moud_on <- sapply(stop_start_membership, function(x) rbinom(1,size = 1, prob = likelihood_on_moud[x])) ## is the individual on MOUD
  
  ### the simulation results get stored here 
  mortality_matrix <- matrix(NA,n,time_steps)### mortality status: is the person alive or not? 1 = alive, 0 = deceased
  active_matrix <- matrix(NA,n,time_steps) ### active user: is the person currently an active user 1 = active, 0 = inactive
  hiv_status_matrix <- matrix(0,n,time_steps) ### hiv status: what is the individuals active status 1 = positive, 0 = negative 
  hiv_supp_matrix <- matrix(0,n,time_steps) ### hiv suppression: if infected, are they virially supressed? 1 = yes, 0 = no
  hiv_diag_matrix <- matrix(0,n,time_steps) ### hiv diag: has the individual been diagnoised with HIV or not  1 = yes, 0 = no
  
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
  print(start_infected)
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
            wiw[inj_partners[which(prob_hiv_trans_inj_partner == 1)],infected_index] = 1
          }
        }
      }
    }
    hiv_status_matrix[prev_infected,ii] = 1 ### anyone who was previously infected stays infected 
    
    
    
    ### testing - you either will be randomly tested, RCD tested or location based testing
    baseline_testing <- rbinom(n, size = 1, prob = hiv_test_prob) ## right now there is no contact tracing increase in likelihood of being tested 
    rcd_tested = rbinom(n, size = 1, prob = rcd_prob)
    loc_tested = sapply(1:n, function(x) ifelse(is.element(location_vec[x], testing_loc), 1,0))#rbinom(1,size=1,prob = 0.8), 0))
    diag_results <- sapply(1:n, function(x) ifelse((baseline_testing[x] + rcd_tested[x] + loc_tested[x]) >=1, sample(c(1, 2, 3), size = 1, prob = c(baseline_testing[x], rcd_tested[x], loc_tested[x])), 0)) ## 1 = random; 2 = rcd; 3 = location; 0 = not tested 
    
    new_diagnoised <- which(hiv_status_matrix[,ii] == 1 & diag_results > 0 & hiv_diag_matrix[,(ii-1)] == 0)
    
    hiv_diag_matrix[new_diagnoised,ii] = diag_results[new_diagnoised]
    
    ## anyone who was previously diagnosed stays previously diagnosed 
    prev_diag_indiv <- which(hiv_diag_matrix[,(ii-1)] > 0)
    prev_diag_route <- hiv_diag_matrix[prev_diag_indiv,(ii-1)]
    hiv_diag_matrix[prev_diag_indiv,ii] = prev_diag_route
    
    diagnoised_hiv <- which(hiv_diag_matrix[,ii] > 0) 
    n_diagnoised <- length(diagnoised_hiv)
    if(n_diagnoised > 0){ ## there are people who are diagnosed
      on_art_and_supress <- rbinom(n_diagnoised, size = 1, prob = art_update*achieve_viral_supress) 
      hiv_supp_matrix[diagnoised_hiv[which(on_art_and_supress == 1)],ii] = 1
      hiv_supp_matrix[diagnoised_hiv[which(on_art_and_supress == 0)], ii] = 0
    }
    
    ## update demography -- people stopping, dying, start
    pop_status <- pop_matrix[,prev_step] 
    prev_mortality <- mortality_matrix[,prev_step]
    mortality_vec <- sapply(pop_status, function(x) rbinom(1, size = 1, prob = 1-(fixed_od_rate))) ## this is wrong , if this is a 0 then the person has died 
    mortality_matrix[,ii] = mortality_vec
    mortality_matrix[which(prev_mortality == 0),ii] = 0 ## if you were previously dead you remain dead 
    active_status <- active_matrix[,prev_step]
    
    update_status = sapply(active_status_prob[,ii], function(x) rbinom(1, size = 1, prob = x))
    active_matrix[,ii] = update_status
  }
  return(list(mortality_matrix = mortality_matrix, active_matrix = active_matrix, hiv_status_matrix = hiv_status_matrix, hiv_supp_matrix = hiv_supp_matrix, hiv_diag_matrix = hiv_diag_matrix, inj_events_per_indiv = inj_events_per_indiv, mean_inj_events_per_partner = mean_inj_events_per_partner, wiw = wiw))
}


### YOU MUST MAKE THESE MATRICES TO RUN THE MAIN CODE
stop_start_param <- make_stop_start_matrix(n_years, time_steps, exp_model_early, exp_model_delay, exp_model_persistent, n, age_vec)

active_status_prob = stop_start_param$active_status_prob
stop_start_membership = stop_start_param$membership_vec

skew_0 <- c(2,2) ## alpha, beta values for beta distribution 
skew_1 <- c(1,5)
skew_2 <- c(1,7)
skew_type <- matrix(,3,2)
skew_type[1,] = skew_0
skew_type[2,] = skew_1
skew_type[3,] = skew_2

lambda_vec <- c(30, 60, 100, 150)
base_risk_cat <- sample(1:4, size = n, prob = c(0.25, 0.25, 0.25, 0.25), replace = TRUE) 
base_risk_group <- matrix(base_risk_cat,n,time_steps)
base_skew_matrix<-matrix(1,n,time_steps)

syring_share_on <- rep(0.75, time_steps) #c(rep(0, time_steps/2), rep(0, time_steps/2)) ## index if syring sharing services are available 
moud_on <- rep(1, time_steps) #c(rep(0, time_steps/2), rep(0, time_steps/2)) ## index if MOUD services are available 1 if it is available, 0 if not 

rcd_prob = 0.1

n_sim <- 20
testing_loc = rep(1, n_locations)
total_infected_all_loc_testing <- rep(NA, n_sim)
for(ii in 1:n_sim){
  sim_base_results <- run_simulation_loc(base_risk_group, base_skew_matrix, syring_share_on, moud_on, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec, location_vec, testing_loc, rcd_prob)
  total_infected_all_loc_testing[ii] = sum(sim_base_results$hiv_status_matrix[,60], na.rm=T)
}


testing_loc = rep(0, n_locations)
testing_loc[c(1,2,3)] = 1
total_infected_all_loc_testing_2 <- rep(NA, n_sim)
for(ii in 1:n_sim){
  sim_base_results_2 <- run_simulation_loc(base_risk_group, base_skew_matrix, syring_share_on, moud_on, inj_network, sexual_network, active_status_prob, stop_start_membership, lambda_vec, location_vec, testing_loc, rcd_prob)
  total_infected_all_loc_testing_2[ii] = sum(sim_base_results_2$hiv_status_matrix[,60], na.rm=T)
}

par(mfrow=c(1,1))
boxplot(total_infected_all_loc_testing, total_infected_all_loc_testing_2, names = c('all loc', '3 loc'))

image(sim_base_results_2$wiw)


### now we will simulate what would have happened with additional testing strategies 

simulate_other_testing <- function(location_vec, hiv_status_matrix, hiv_diag_matrix, wiw, testing_loc, rcd_prob){
  ### location_vec: vector of locations assigned to an individual
  ### hiv_status_matrix: col -> time points; rows -> individuals - HIV status per time point; 0 = uninfected; 1 = infected;
  ### hiv_diag_matrix: col -> time points; rows -> individuals - diagnosis status; 0 = undiagnosed; 1 = undiagnosed
  ### wiw: who infected whom> col -> infector; rows -> person infected; 1 = they infected each other; 0 = initial seed infection 
  ### testing_loc: vector with length equal to number of locations. If 1, they are a place that can have location based testing
  ### rcd_ratio: probability of a rcd testing strategy for positive as an initial event 
  
### This is not very realistic but the way it will currently work is: it will add in additional testing to see if you are tested at a location or a rcd seed
  time_steps <- ncol(hiv_status_matrix)
  testing_loc = seq(1,length(testing_loc))[which(testing_loc == 1)] ## make a list of the ids of locations with testing 
  update_hiv_diag_matrix <- hiv_diag_matrix
  update_hiv_status_matrix <- hiv_status_matrix
  for(ii in 1:time_steps){
    current_infected = which(update_hiv_status_matrix[,ii] == 1)
    n_infected = length(current_infected)
    current_diag = hiv_diag_matrix[current_infected,ii]
    updated_diag = current_diag
    rcd_tested = rbinom(n_infected, size = 1, prob = rcd_prob)
    loc_tested = sapply(1:n_infected, function(x) ifelse(is.element(location_vec[x], testing_loc), rbinom(1,size=1,prob = 0.8), 0))
    diag_results <- sapply(1:n_infected, function(x) ifelse((current_diag[x] + rcd_tested[x] + loc_tested[x]) >=1, sample(c(1, 2, 3), size = 1, prob = c(current_diag[x], rcd_tested[x], loc_tested[x])), 0)) ## 1 = random; 2 = rcd; 3 = location; 0 = not tested 
    update_hiv_diag_matrix[current_infected,ii] = diag_results
  }
  return(update_hiv_diag_matrix)
}

find_infection_chain <- function(wiw, seed, time_point_diag, update_hiv_status_matrix){
  ### seed is the person who is now diagnosed so they won't be infected their contacts for subsequent time steps
  infected_indiv = which(wiw[,seed] == 1) ## who was infected by the seed - need to double check if this is the column or the row that you need to check 
  infection_times_for_infected_indiv <- sapply(infected_indiv, function(x) min(which(update_hiv_status_matrix[x,]>0)))
  change_infection_indiv <- infected_indiv[which(infection_times_for_infected_indiv > time_point_diag)] ## this individual was infected after you were diagnoised and will not be uninfected
  update_hiv_status_matrix[change_infection_indiv]
  
}

update_hiv_diag_matrix <- simulate_other_testing(location_vec, sim_base_results$hiv_status_matrix, sim_base_results$hiv_diag_matrix, sim_base_results$wiw, testing_loc, rcd_prob)



