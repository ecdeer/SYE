## linear regression model: Y_i = \beta_0 + \beta_1 + \epsilon_i,
## where \epsilon_i ~iid N(0, \sigma^2)

library(tidyverse)
library(palmerpenguins)
penguins_comp <- penguins |> filter(!is.na(bill_length_mm) & !is.na(bill_depth_mm))
lm(bill_depth_mm ~ bill_length_mm, data = penguins_comp) |>
  summary()

## for simpliciity,
## suppose that we just fix \beta_0 at it's estimator from lm(), 20.885 and we fix
## sigma^2 at it's estimator from lm(), 1.922^2
## then we can use Bayesian to obtain a posterior distribution for \beta_1


#### PRIOR FOR \beta_1: diffuse normal: N(0, 1000)
#### helpful to work on log scale throughout: otherwise the numbers get really really tiny
compute_prior_log <- function(beta1) {
  dnorm(beta1, mean = 0, sd = 1000, log = TRUE)
}

compute_prior_log(beta1 = 0)
compute_prior_log(beta1 = 1)


#### LIKELIHOOD/DATA: assume Y_1, ... Y_n ~ iid N(\beta_0 + beta_1 X_i, \sigma^2),
#### where \beta_0 and \sigma^2 are assumed to be known.
beta0 <- 20.885
sigma2 <- 1.922^2

compute_likelihood_log <- function(beta1) {

  dnorm(penguins_comp$bill_depth_mm, mean = beta0 + beta1 * penguins_comp$bill_length_mm,
        sd = sqrt(sigma2), log = TRUE) |>
    sum() ## sum log likelihoods instead of multiplying like you would do with non-logged likelihoods
  ## can think back to the math-stat part of probability- the mean just differs depending on
  ## the bill length instead of just being constant.
}
compute_likelihood_log(beta1 = 0)
compute_likelihood_log(beta1 = -1)
compute_likelihood_log(beta1 = 1)
## higher likelihoods are "better"


## combining prior and data/likelihood
## again, summing not multiplying because we add things on the log scale
compute_posterior_log <- function(beta1) {
  compute_prior_log(beta1) + compute_likelihood_log(beta1)
}

compute_posterior_log(beta1 = 0)
compute_posterior_log(beta1 = -1)
compute_posterior_log(beta1 = 1)


## define a proposal/jumping distribution to be N(0, 0.1)


## set-up the chain
niter <- 10000
beta1_store <- rep(NA, 10000) ## empty vector

beta1_store[1] <- 1.2 ## define a starting value

for (i in 2:niter) {
  beta1_current <- beta1_store[i - 1]

  ## propose a new value with proposal/jumping distribution
  ## here, my proposal/jumping distribution is N(current_beta1, 0.1)
  ## this distribution must be symmetric around the current beta1 (or else would need adjustment)
  proposal_sd <- 0.01
  ## proposal_sd can't be too large or too small, as discussed in the reading.

  beta1_proposal <- rnorm(1, mean = beta1_current, sd = proposal_sd)

  ## compute acceptance probability
  alpha <- min(1, exp(compute_posterior_log(beta1_proposal) - compute_posterior_log(beta1_current)))

  ## maybe update beta1 with acceptance probability alpha
  beta1_store[i] <- sample(c(beta1_proposal, beta1_current), size = 1,
                           prob = c(alpha, 1 - alpha))
}

## plots
## filter removes "start-up"
plot_df <- tibble(iter = 1:niter, beta1 = beta1_store) |>
  filter(iter > 2000)

## make sure that this seems "stable" (isn't trailing off in an increasing
## or decreasing direction)
ggplot(data = plot_df, aes(x = iter, y = beta1)) +
  geom_line()

plot_df |>
  summarise(mean = mean(beta1),
            sd = sd(beta1))
## check these with estimator for beta1 and its SE from lm: they should be close!

## shape of posterior
ggplot(data = plot_df, aes(x = beta1)) +
  geom_histogram(colour = "skyblue4", fill = "skyblue1", bins = 15)
