# libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tidyverse)
library(readxl)

# get and read data
clams <- read_csv("data/clams_clean.csv")
covariate_data <- read_excel("data/covariate_data.xlsx")
clams <- clams |> left_join(covariate_data, by = "location")
clams

# remove negative values from Coopers Falls, Lake Clear, Mirror Lake
clams <- clams |> 
  filter(growth > 0) |>
  mutate(x = case_when(
    body == "L" ~ 1,
    body == "R" ~ 0)) 
  |> group_by(name) 

# vb function

vb_function <- function(t, max_l, k, b, beta, x) {
  max_l * (1 - exp(-(k + beta * x) * (t - b)))}

# data for union falls pond
# clams_unionfallspond <- clams |>
#  filter(location == "UnionFallsPond") |>
#  group_by(name)

# data for rivers
#rivers <- clams |> 
#  filter(body == "R") |>
#  group_by(name)

# data for lakes
# lakes <- clams |> 
#  filter(body == "L") |>
#  group_by(name)

# MCMC Algorithm

# step 2 - compute prior log

compute_prior_log <- function(max_l, k, b, sigma2, beta) {
  dnorm(max_l, mean = 0, sd = 1000, log = TRUE) +
    dnorm(k, mean = 0, sd = 1000, log = TRUE) +
    dnorm(b, mean = 0, sd = 1000, log = TRUE) +
    log(1 / sigma2) +
    dnorm(beta, mean = 0, sd = 1000, log = TRUE)
}

# step 4 - compute log likelihood (make function)

compute_likelihood_log <- function(max_l, k, b, sigma2, beta) {
  dnorm(rivers$Length,
        mean =  vb_function(rivers$growth_period, max_l, k, b, beta, rivers$x),
        sd = sqrt(sigma2), log = TRUE) |>
    sum()
}

# step 5 - sample call

compute_likelihood_log(max_l = 30, k = 0.3, b = -0.8, sigma2 = 5, beta = 0)
compute_likelihood_log(max_l = 40, k = 0.35, b = -0.5, sigma2 = 30, beta = 5)
compute_likelihood_log(max_l = 50, k = 0.4, b = -0.1, sigma2 = 50, beta = 10)


# step 6 - compute posterior log

compute_posterior_log <- function(max_l, k, b, sigma2, beta) {
  compute_prior_log(max_l, k, b, sigma2, beta) + 
    compute_likelihood_log(max_l, k, b, sigma2, beta)
}

# step 7 - sample call

compute_posterior_log(max_l = 30, k = 0.3, b = -0.8, sigma2 = 5, beta = 0)
compute_posterior_log(max_l = 40, k = 0.35, b = -0.5, sigma2 = 30, beta = 5)
compute_posterior_log(max_l = 50, k = 0.4, b = -0.1, sigma2 = 50, beta = 10)


## define a proposal/jumping distribution to be N(0, 0.1)

## set-up the chain
niter <- 10000

# set up three empty vectors, one for each parameter
max_l_store <- rep(NA, 10000)
b_store <- rep(NA, 10000)
k_store <- rep(NA, 10000)
sigma2_store <- rep(NA, 10000)
beta_store <- rep(NA, 10000)
  
  
# define a starting value
max_l_store[1] <- 1.2
k_store[1] <- 0.3183
b_store[1] <- -0.8749
sigma2_store[1] <- 30
beta_store[1] <- 0
  
  
for (i in 2:niter) {
    
    ## start with max_l
    max_l_current <- max_l_store[i - 1]
    
    proposal_sd <- 3
    
    max_l_proposal <- rnorm(1, mean = max_l_current, sd = proposal_sd)
    
    alpha_l <- min(1, exp(compute_posterior_log(max_l_proposal, k_store[i-1], b_store[i-1], sigma2_store[i-1], beta_store[i-1]) - compute_posterior_log(max_l_current, k_store[i-1], b_store[i-1], sigma2_store[i-1], beta_store[i-1])))
    
    max_l_store[i] <- sample(c(max_l_proposal, max_l_current), size = 1,
                             prob = c(alpha_l, 1 - alpha_l))
    
    
   ## start K
    
    k_current <- k_store[i - 1]
    
    proposal_sd <- 0.05
    
    k_proposal <- rnorm(1, mean = k_current, sd = proposal_sd)
    
    
    alpha_k <- min(1, exp(compute_posterior_log(max_l_store[i], k_proposal, b_store[i-1], sigma2_store[i-1], beta_store[i-1]) - compute_posterior_log(max_l_store[i], k_current, b_store[i-1], sigma2_store[i-1], beta_store[i-1])))
    
    k_store[i] <- sample(c(k_proposal, k_current), size = 1,
                         prob = c(alpha_k, 1 - alpha_k))
    
    
    # start B
    
    b_current <- b_store[i - 1]
    
    proposal_sd <- 0.1
    
    b_proposal <- rnorm(1, mean = b_current, sd = proposal_sd)
    
    alpha_b <- min(1, exp(compute_posterior_log(max_l_store[i], k_store[i], b_proposal, sigma2_store[i-1], beta_store[i-1]) - compute_posterior_log(max_l_store[i], k_store[i], b_current, sigma2_store[i-1], beta_store[i-1])))
    
    b_store[i] <- sample(c(b_proposal, b_current), size = 1,
                         prob = c(alpha_b, 1 - alpha_b))
    
    # start sigma2
    
    sigma2_current <- sigma2_store[i-1]
    sigma2_current_log <- sigma2_current |> log()
    
    proposal_sd <- 0.2
    
    sigma2_proposal_log <- rnorm(1, mean = sigma2_current_log, sd = proposal_sd)
    sigma2_proposal <- exp(sigma2_proposal_log)
    
    alpha_sigma2 <- min(1, exp(compute_posterior_log(max_l_store[i], k_store[i], b_store[i], sigma2_proposal, beta_store[i-1]) + log(sigma2_proposal) - compute_posterior_log(max_l_store[i], k_store[i], b_store[i], sigma2_current, beta_store[i-1]) - log(sigma2_current)))
    
    sigma2_store[i] <- sample(c(sigma2_proposal, sigma2_current), size = 1,
                              prob = c(alpha_sigma2, 1 - alpha_sigma2))
    
    # start beta
    beta_current <- beta_store[i - 1]
    
    proposal_sd <- 0.1
    
    beta_proposal <- rnorm(1, mean = beta_current, sd = proposal_sd)
    
    alpha_beta <- min(1, exp(compute_posterior_log(max_l_store[i], k_store[i], b_store[i], sigma2_store[i], beta_proposal) - compute_posterior_log(max_l_store[i], k_store[i], b_store[i], sigma2_store[i], beta_current)))
    
    beta_store[i] <- sample(c(beta_proposal, beta_current), size = 1,
                            prob = c(alpha_beta, 1 - alpha_beta))
    
    ## updating one at a time
}
  
  
plot_df <- tibble(iter = 1:niter, max_l = max_l_store, k = k_store, b = b_store,
                  sigma2_store, beta_store) |> filter(iter > 2000)

# max_l plots
ggplot(data = plot_df, aes(x = iter, y = max_l)) +
  geom_line()

plot_df |>
  summarise(mean_l = mean(max_l),
            sd = sd(max_l))

ggplot(data = plot_df, aes(x = max_l)) +
  geom_histogram(colour = "skyblue4", fill = "skyblue1", bins = 15)

# k plots
ggplot(data = plot_df, aes(x = iter, y = k)) +
  geom_line()

plot_df |>
  summarise(mean_k = mean(k),
            sd = sd(k))

ggplot(data = plot_df, aes(x = k)) +
  geom_histogram(colour = "skyblue4", fill = "skyblue1", bins = 15)

# b plots
ggplot(data = plot_df, aes(x = iter, y = b)) +
  geom_line()

plot_df |>
  summarise(mean_b = mean(b),
            sd = sd(b))

ggplot(data = plot_df, aes(x = b)) +
  geom_histogram(colour = "skyblue4", fill = "skyblue1", bins = 15)

# sigma2 plots
ggplot(data = plot_df, aes(x = iter, y = sigma2_store)) +
  geom_line()

plot_df |>
  summarise(mean_sigma2 = mean(sigma2_store),
            sd = sd(sigma2_store))

ggplot(data = plot_df, aes(x = sigma2_store)) +
  geom_histogram(colour = "skyblue4", fill = "skyblue1", bins = 15)

# beta plots
ggplot(data = plot_df, aes(x = iter, y = beta_store)) +
  geom_line()

plot_df |>
  summarise(mean_sigma2 = mean(beta_store),
            sd = sd(beta_store))

ggplot(data = plot_df, aes(x = beta_store)) +
  geom_histogram(colour = "skyblue4", fill = "skyblue1", bins = 15)







