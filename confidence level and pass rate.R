# Load and compare the simulation results by each method
# Compare the pass rates in a summary table

library(gt)
library(gtExtras)
library(tidyverse)
library(here)

# set confidence level
CL <- 0.95

# methods to compare
methods = c("Exact", "Delta_Normality", "Raw_Nonparametric", "Raw_Parametric", "Reflected", "Bootstrape_t")

# helper function to read the simulation results
read_all <- function(path, shape){
  f_files <- list.files(path, pattern = "^\\w.*\\.rds$")
  
  df <- tibble()
  for (f in f_files) {
    print(f)
    df_tmp <- readRDS(here(path, f))
    df <- bind_rows(df, df_tmp)
  }
  
  df <- mutate(df, shape = shape) %>% arrange(rowID)
  
  return(df)
}

# read simulation results 
df_1 <- read_all(here("method_3-4_CLH", "shape2"), 2) 
df_2 <- read_all(here("method_3-4_CLH", "shape4"), 4)
df_3 <- read_all(here("method_3-4_CLH", "shape05"), 0.5)

# combine the results 
df_all <- bind_rows(df_1, df_2, df_3) %>% 
  mutate(
    method = fct_relevel(method, methods),
    correct_dist = if_else(correct_dist, "Correct Dist", "Incorrect Dist")
  )


# filter for coverage_rate >= 0.975
gt1 <- df_all %>% 
  mutate(
    QA = coverage_rate >= 0.975
  ) %>% 
  group_by(shape, method, correct_dist) %>% 
  summarise(pass = sum(QA)) %>% 
  mutate(pass = paste0(pass, "/12")) %>% 
  pivot_wider(
    names_from = shape,
    values_from = pass,
    names_prefix = "shape_"
  ) %>% 
  gt() %>% 
  tab_header(
    title = "coverage_rate >= 0.975",
    subtitle = "(number of cases)"
  ) %>% 
  cols_align("left")


# filter for coverage_rate >= CL
gt2 <- df_all %>% 
  mutate(
    QA = coverage_rate >= CL
  ) %>% 
  group_by(shape, method, correct_dist) %>% 
  summarise(pass = sum(QA)) %>% 
  mutate(pass = paste0(pass, "/12")) %>% 
  pivot_wider(
    names_from = shape,
    values_from = pass,
    names_prefix = "shape_"
  ) %>% 
  gt() %>% 
  tab_header(
    title = "coverage_rate >= 0.95",
    subtitle = "(number of cases)"
  ) %>% 
  cols_align("left")


# filter for coverage_rate >= 0.90
gt3 <- df_all %>% 
  mutate(
    QA = coverage_rate >= 0.9
  ) %>% 
  group_by(shape, method, correct_dist) %>% 
  summarise(pass = sum(QA)) %>% 
  mutate(pass = paste0(pass, "/12")) %>% 
  pivot_wider(
    names_from = shape,
    values_from = pass,
    names_prefix = "shape_"
  ) %>% 
  gt() %>% 
  tab_header(
    title = "coverage_rate >= 0.9",
    subtitle = "(number of cases)"
  ) %>% 
  cols_align("left")

# comparison table
gt_two_column_layout(list(gt1, gt2))
gt3




