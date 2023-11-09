# Monte Carlo Simulation Study for Estimators and CI Performance
Create a Monte Carlo simulation study in R to investigate properties of estimators and confidence intervals (CIs).

- Monte Carlo Simulation 
- Estimators 
- Confidence intervals
- Parametric and Non-parametric Bootstrapping

[Report for this project](./proj1 v3.Rmd)

Tuesday, September 27, 2022, 11:59 PM

Chien-Lan Hsueh (chienlan.hsueh at gmail.com)


## Overview and Project Goal
The project involves creating a Monte Carlo simulation study in R (generating data in R to investigate properties of estimators/CIs) and the creation of a report.

We’ll simulate data from different distributions and look at the performance of confidence intervals for capturing a parameter. Here: 1/(mean of the distribution).
Specifically, we’ll consider obtain a random sample from an $\exp{(λ)}$ distribution and attempt to capture $1/E(Y ) = \lambda$ and investigate the following performance aspects of our confidence intervals:

- coverage rate (the proportion of time they contain the true value)
- proportion of intervals that miss by being below (all values less than the truth)
- proportion of intervals that miss by being above (all values greater than the truth)
- average length of the interval
- standard error of the average length of the interval

## Part 1 - Compare Different MC Methods
Assume that we have an iid sample from an exp(λ) population and we have an interest in creating
confidence intervals for one divided by the mean of the distribution, here λ. We’ll compare six different
methods for constructing a confidence interval:

- The exact interval using the estimator 1/Y
- The large-sample normality based interval using the CLT and the Delta method on that same estimator
- The raw percentile non-parametric bootstrap interval using that estimator
- The raw percentile parametric bootstrap interval using that estimator
- The reflected percentile non-parametric bootstrap interval using that estimator
- The bootstrap t-interval using a non-parametric bootstrap using that estimator

## Part 2 - Data Simulation
Compare the above intervals’ performance when truly sampling from an exp(λ) distribution
for sample sizes of n = 10, 30, 100, and 500 and λ values of 0.5, 1, and 5. (All combinations, so 12 total
combinations.) You should simulate at least 1000 samples for each combination.

## Part 3 - Performance
Compare the above intervals’ performances when we are wrong about the data generating
process (to see how robust our methods are to misspecification). That is, generate data from a gamma(2, λ)
distribution under the different n and λ specifications above. 
In this case, you want to check if 1/E(Y ) = λ/2 is in the above intervals rather than λ itself. 