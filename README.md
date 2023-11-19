# Monte Carlo Simulation Study for Estimators and CI Performance

Create a Monte Carlo simulation study in R to investigate properties of estimators and confidence intervals (CIs).

[Project report](https://rpubs.com/clh2021/1113711)

Key features:

- Monte Carlo Simulation 
- Estimators and Confidence Intervals
- Large-sample normality based intervals using CLT and Delta method 
- Parametric and Non-parametric Bootstrapping
- R code benchmarks

R packages used:

- `here`: enables easy file referencing and builds file paths in a OS-independent way
- `stats`: loads this before loading `tidyverse` to avoid masking some `tidyverse` functions
- `tidyverse`: includes collections of useful packages like `dplyr` (data manipulation), `tidyr` (tidying data),  `ggplots` (creating graphs), etc.
- `scales`: formats and labels scales nicely for better visualization

## Project Report

[Project report](https://rpubs.com/clh2021/1113711) ([Github Markdown](./proj1_v3.md))([R Markdown](./proj1_v3.Rmd))

The analysis results with all theoretical backgrounds and math derivations are included.

Author: Chien-Lan Hsueh (chienlan.hsueh at gmail.com)

## Overview and Project Goal
The project involves creating a Monte Carlo simulation study in R (generating data in R to investigate properties of estimators and CIs) and the creation of a report.

We’ll simulate data from different distributions and look at the performance of confidence intervals for capturing a parameter. Here: 1/(mean of the distribution).
Specifically, we’ll consider obtain a random sample from an $\exp{(\lambda)}$ distribution and attempt to capture $1/E(Y) = \lambda$ and investigate the following performance aspects of our confidence intervals:

- coverage rate (the proportion of time they contain the true value)
- proportion of intervals that miss by being below (all values less than the truth)
- proportion of intervals that miss by being above (all values greater than the truth)
- average length of the interval
- standard error of the average length of the interval

## Part 1 - Compare Different MC Methods
Assume that we have an iid sample from an $\exp{(\lambda)}$ population and we have an interest in creating
confidence intervals for one divided by the mean of the distribution, here $\lambda$. We’ll compare six different methods for constructing a confidence interval:

1. The exact interval using the estimator $1/Y$
1. The large-sample normality based interval using the CLT and the Delta method on that same estimator
1. The raw percentile non-parametric bootstrap interval using that estimator
1. The raw percentile parametric bootstrap interval using that estimator
1. The reflected percentile non-parametric bootstrap interval using that estimator
1. The bootstrap t-interval using a non-parametric bootstrap using that estimator

## Part 2 - Data Simulation
Compare the above intervals’ performance when truly sampling from an $\exp{(\lambda)}$ distribution
for sample sizes of $n =$ 10, 30, 100 and 500 and $\lambda$ values of 0.5, 1, and 5. (All combinations, so 12 total combinations.) Simulate at least 1000 samples for each combination.

## Part 3 - Performance Comparison
Compare the above intervals’ performances when we are wrong about the data generating
process (to see how robust our methods are to misspecification). That is, generate data from a $\Gamma(2,  \lambda)$ distribution under the different $n$ and $\lambda$ specifications above. 
In this case, you want to check if $1/E(Y) = \lambda/2$ is in the above intervals rather than $\lambda$ itself. 
