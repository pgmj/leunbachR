library(tidyverse)
set.seed(87365)

# import code
files <- list.files(path = "R")
for (i in files) {
  source(paste0("R/",i))
}

# data import --------------------------------------------------------------

d3a <- read_delim("FHSdata/3a_Dataset 3 - AbilityBias = 0.3, DifficultyBias = 0, n_item = 10, n = 500.csv") %>% 
  as.matrix()
d3b <- read_delim("FHSdata/3b_Dataset 3 - AbilityBias = 0.3, DifficultyBias = 0, n_item = 10, n = 1500.csv") %>% 
  as.matrix()
d3c <- read_delim("FHSdata/3c_Dataset 3 - AbilityBias = 0.3, DifficultyBias = 0, n_item = 20, n = 500.csv") %>% 
  as.matrix()
d3d <- read_delim("FHSdata/3d_Dataset 3 - AbilityBias = 0.3, DifficultyBias = 0, n_item = 20, n = 1500.csv") %>% 
  as.matrix()


# Fit IPF -----------------------------------------------------------------


## d3a ---------------------------------------------------------------------

fit3a <- leunbach_ipf(d3a, verbose = TRUE)
summary(fit3a)
print(fit3a)
calculate_statistics(fit3a)

#Expected Gamma =  -0.050  s.e. =  0.0395
#Observed Gamma =  -0.046    p  =  0.5405 (one-sided p.value)

## d3b ---------------------------------------------------------------------

fit3b <- leunbach_ipf(d3b, verbose = TRUE)
summary(fit3b)
print(fit3b)
calculate_statistics(fit3b)

#Expected Gamma =  -0.038  s.e. =  0.0230
#Observed Gamma =  -0.033    p  =  0.5927 (one-sided p.value)

## d3c ---------------------------------------------------------------------

fit3c <- leunbach_ipf(d3c, verbose = TRUE)
summary(fit3c)
print(fit3c)
calculate_statistics(fit3c)

#Expected Gamma =  -0.053  s.e. =  0.0346
#Observed Gamma =  -0.058    p  =  0.4468 (one-sided p.value)

df3c <- c(0L, 0L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 17L, 16L, 15L, 14L, 13L, 12L, 11L, 10L, 9L, 8L, 7L, 6L, 5L, 4L, 3L, 2L, 1L, 0L)
sum(df3c) - length(df3c)

## d3d ---------------------------------------------------------------------

fit3d <- leunbach_ipf(d3d, verbose = TRUE)
summary(fit3d)
print(fit3d)
calculate_statistics(fit3d)


# Orbits ------------------------------------------------------------------

## d3a ---------------------------------------------------------------------

orb3a <- analyze_orbits(fit3a, alpha = 0.05)
summary(orb3a)
print(orb3a)
# View orbit distribution for a specific total score
get_orbit(orb3a, total_score = 5)

# Plot orbit distributions
plot(orb3a, type = "orbits")
plot(orb3a, type = "cumulative")
plot(orb3a, type = "significant")

## d3b ---------------------------------------------------------------------

orb3b <- analyze_orbits(fit3b, alpha = 0.05)
summary(orb3b)
print(orb3b)

## d3c ---------------------------------------------------------------------

orb3c <- analyze_orbits(fit3c, alpha = 0.05)
summary(orb3c)
print(orb3c)

## d3d ---------------------------------------------------------------------

orb3d <- analyze_orbits(fit3d, alpha = 0.05)
summary(orb3d)
print(orb3d)

### 
### NOTE: check df, and add gamma test
###


# Direct equating ---------------------------------------------------------

## d3a ---------------------------------------------------------------------
leunbach_equate(fit3a, verbose = T)
leunbach_equate(fit3a, method = "newton")

## d3b ---------------------------------------------------------------------
leunbach_equate(fit3b)
leunbach_equate(fit3b, method = "newton")
leunbach_equate(fit3b, direction = "2to1")

## d3c ---------------------------------------------------------------------
leunbach_equate(fit3c)
leunbach_equate(fit3c, method = "newton")

## d3d ---------------------------------------------------------------------
leunbach_equate(fit3d)
leunbach_equate(fit3d, method = "newton")


# Bootstrap DE ------------------------------------------------------------

## d3a ---------------------------------------------------------------------
boot3a <- leunbach_bootstrap(fit3a, n_cores = 12, verbose = TRUE)
print(boot3a)
summary(boot3a)

# Significance of the likelihood ratio test:  p = 0.905
# Significance of gamma coefficient:          p = 0.689

## d3b ---------------------------------------------------------------------
boot3b <- leunbach_bootstrap(fit3b, n_cores = 12, verbose = TRUE)
print(boot3b)
summary(boot3b)

## d3c ---------------------------------------------------------------------
boot3c <- leunbach_bootstrap(fit3c, n_cores = 12, verbose = TRUE)
print(boot3c)
summary(boot3c)

## d3d ---------------------------------------------------------------------
boot3d <- leunbach_bootstrap(fit3d, n_cores = 12, verbose = TRUE)
print(boot3d)
summary(boot3d)


# Indirect equating -------------------------------------------------------

## d1 ---------------------------------------------------------------------
d1 <- read_delim("FHSdata/indirekt/Dataset 1 - AbilityBias = 0.3, DifficultyBias = 0, n_item = 10, n = 1000.csv") %>% 
  as.matrix()
# Fit models
fit_ab <- leunbach_ipf(d1[,c(1,2)], verbose = TRUE)
fit_bc <- leunbach_ipf(d1[,c(2,3)], verbose = TRUE)

# Direct equating (for comparison)
eq_ab <- leunbach_equate(fit_ab, direction = "1to2")  # A → B
eq_bc <- leunbach_equate(fit_bc, direction = "1to2")  # B → C

# Indirect equating:  A → B → C
indirect1 <- leunbach_indirect_equate(fit_ab, fit_bc,
                                      direction_ab = "1to2",
                                      direction_bc = "1to2",
                                      verbose = TRUE)
print(indirect1)


### Bootstrap ---------------------------------------------------------------

boot_indirect1 <- leunbach_indirect_bootstrap(fit_ab, fit_bc,
                                              direction_ab = "1to2",
                                              direction_bc = "1to2",
                                              nsim = 100,
                                              verbose = TRUE, n_cores = 12)
print(boot_indirect1)
summary(boot_indirect1)

# Plots
plot(boot_indirect1, type = "equating")
plot(boot_indirect1, type = "see")

# Get clean table
indirect_table <- get_indirect_equating_table(boot_indirect1)
indirect_table

ind1see <- c(0.46, 0.55, 0.4, 0.25, 0.15, 0.04, 0.04, 0.27, 0.51, 0.51, 0)
mean(ind1see[])
ind1seeDIGRAM <- c(0.51, 0.56, 0.39, 0.28, 0.13, 0.08, 0.05, 0.25, 0.51, 0.51, 0.27)
mean(ind1seeDIGRAM)

## d2 ---------------------------------------------------------------------
d2 <- read_delim("FHSdata/indirekt/Dataset 2 - AbilityBias = 0.3, DifficultyBias = 0, n_item = 20, n = 1000.csv") %>% 
  as.matrix()
# Fit models
fit_ab <- leunbach_ipf(d2[,c(1,2)], verbose = TRUE)
fit_bc <- leunbach_ipf(d2[,c(2,3)], verbose = TRUE)

# Direct equating (for comparison)
eq_ab <- leunbach_equate(fit_ab, direction = "1to2")  # A → B
eq_bc <- leunbach_equate(fit_bc, direction = "1to2")  # B → C

# Indirect equating:  A → B → C
indirect2 <- leunbach_indirect_equate(fit_ab, fit_bc,
                                      direction_ab = "1to2",
                                      direction_bc = "1to2",
                                      verbose = TRUE)
print(indirect2)

### Bootstrap ---------------------------------------------------------------

boot_indirect2 <- leunbach_indirect_bootstrap(fit_ab, fit_bc,
                                              direction_ab = "1to2",
                                              direction_bc = "1to2",
                                              nsim = 500,
                                              verbose = TRUE, n_cores = 12)
print(boot_indirect2)
summary(boot_indirect2)

# Plots
plot(boot_indirect2, type = "equating")
plot(boot_indirect2, type = "see")

# Get clean table
indirect_table2 <- get_indirect_equating_table(boot_indirect2)
indirect_table2

## d3 ---------------------------------------------------------------------
d3 <- read_delim("FHSdata/indirekt/Dataset 3 - AbilityBias = 0.3, DifficultyBias = 0, n_item = 10, n = 3000.csv") %>% 
  as.matrix()
# Fit models
fit_ab <- leunbach_ipf(d3[,c(1,2)], verbose = TRUE)
fit_bc <- leunbach_ipf(d3[,c(2,3)], verbose = TRUE)

# Direct equating (for comparison)
eq_ab <- leunbach_equate(fit_ab, direction = "1to2")  # A → B
eq_bc <- leunbach_equate(fit_bc, direction = "1to2")  # B → C

# Indirect equating:  A → B → C
indirect3 <- leunbach_indirect_equate(fit_ab, fit_bc,
                                      direction_ab = "1to2",
                                      direction_bc = "1to2",
                                      verbose = TRUE)
print(indirect3)


### Bootstrap ---------------------------------------------------------------

boot_indirect3 <- leunbach_indirect_bootstrap(fit_ab, fit_bc,
                                              direction_ab = "1to2",
                                              direction_bc = "1to2",
                                              nsim = 500,
                                              verbose = TRUE, n_cores = 12)
print(boot_indirect3)
summary(boot_indirect3)

# Plots
plot(boot_indirect3, type = "equating")
plot(boot_indirect3, type = "see")

# Get clean table
indirect_table3 <- get_indirect_equating_table(boot_indirect3)
indirect_table3

## d4 ---------------------------------------------------------------------
d4 <- read_delim("FHSdata/indirekt/Dataset 4 - AbilityBias = 0.3, DifficultyBias = 0, n_item = 20, n = 3000.csv") %>% 
  as.matrix()
# Fit models
fit_ab <- leunbach_ipf(d4[,c(1,2)], verbose = TRUE)
fit_bc <- leunbach_ipf(d4[,c(2,3)], verbose = TRUE)

# Direct equating (for comparison)
eq_ab <- leunbach_equate(fit_ab, direction = "1to2")  # A → B
eq_bc <- leunbach_equate(fit_bc, direction = "1to2")  # B → C

# Indirect equating:  A → B → C
indirect4 <- leunbach_indirect_equate(fit_ab, fit_bc,
                                      direction_ab = "1to2",
                                      direction_bc = "1to2",
                                      verbose = TRUE)
print(indirect4)

#Could you check what eq_ab$equating_table shows for source score 1? 
#Specifically, what is the expected anchor score (column 3) for source = 1?
#Also, what does eq_bc$equating_table show for anchor score 1? What is the expected target score?

### Bootstrap ---------------------------------------------------------------

boot_indirect4 <- leunbach_indirect_bootstrap(fit_ab, fit_bc,
                                              direction_ab = "1to2",
                                              direction_bc = "1to2",
                                              nsim = 200,
                                              verbose = TRUE, n_cores = 12)
print(boot_indirect4)
summary(boot_indirect4)

# Plots
plot(boot_indirect4, type = "equating")
plot(boot_indirect4, type = "see")

# Get clean table
indirect_table4 <- get_indirect_equating_table(boot_indirect4)
indirect_table4




# OLD code ----------------------------------------------------------------


fit3b <- leunbach_ipf(d3b, verbose = TRUE)
summary(fit3b)


fit3d <- leunbach_ipf(d3d, verbose = TRUE,max_iter = 10000,tol = 1e-10)
summary(fit3d)
print(fit3d)
# Analyze orbits
orb3d <- analyze_orbits(fit3d, alpha = 0.05, verbose = TRUE)

# View orbit distribution for a specific total score
get_orbit(orb3d, total_score = 10)
# Plot orbit distributions
plot(orb3d, type = "orbits")
plot(orb3d, type = "cumulative")
plot(orb3d, type = "significant")

# Get summary with significant cells
summary(orb)

eq_1to2 <- leunbach_equate(fit3d, direction = "1to2", verbose = TRUE)


boot_v6_d3d <- leunbach_bootstrap(fit3d, nsim = 1000, verbose = TRUE, seed = 234)
print(boot_v6_d3d)

# bootstrap ---------------------------------------------------------------


lboot <- leunbach_bootstrap(fit, verbose = TRUE)
# View results
print(lboot)
summary(lboot)

# Plot bootstrap distribution of LR
plot(lboot, type = "lr")

# Plot parameters with confidence intervals
plot(lboot, type = "gamma")
plot(lboot, type = "delta")
plot(lboot, type = "sigma")

# Get standard errors for equating
equating_se <- get_equating_se(lboot)
print(equating_se)



# equating ----------------------------------------------------------------

leunbach_equate(fit, verbose = T)



# equating SEE ------------------------------------------------------------

# Fit the model
fit <- leunbach_ipf(d3a, verbose = TRUE)

# Run bootstrap with SEE
boot3 <- leunbach_bootstrap(fit, nsim = 1000, verbose = TRUE, seed = 234)
boot4 <- leunbach_bootstrap(fit, nsim = 1000, verbose = TRUE, seed = 234, see_method = "rmse")
boot5 <- leunbach_bootstrap(fit, nsim = 1000, verbose = TRUE, seed = 234, see_method = "sd")
boot_v4 <- leunbach_bootstrap(fit, nsim = 1000, verbose = TRUE, seed = 234)
boot_v5 <- leunbach_bootstrap(fit, nsim = 1000, verbose = TRUE, seed = 234, see_type = "rounded")
boot_v6 <- leunbach_bootstrap(fit, nsim = 1000, verbose = TRUE, seed = 234)

# View full results with SEE tables
print(boot3)
print(boot_v4)
print(boot_v5)
print(boot_v6)
# Summary
summary(boot)

# Plot SEE
plot(boot, type = "see_1to2")
plot(boot, type = "see_2to1")

# Plot equating functions
plot(boot, type = "equating")

# Just get equating table without bootstrap
eq_1to2 <- leunbach_equate(fit, direction = "1to2", verbose = TRUE)
eq_2to1 <- leunbach_equate(fit, direction = "2to1", verbose = TRUE)


boot_v6_d3d <- leunbach_bootstrap(fit3d, nsim = 1000, verbose = TRUE, seed = 234)
print(boot_v6_d3d)


# testing d3b -------------------------------------------------------------

d3b <- read_delim("FHSdata/3b_Dataset 3 - AbilityBias = 0.3, DifficultyBias = 0, n_item = 10, n = 1500.csv") %>% 
  as.matrix()

fit3b <- leunbach_ipf(d3b, verbose = TRUE)
summary(fit3b)
# orbits
orb3b <- analyze_orbits(fit3b)
orb3b
summary(orb3b)
# Plots
plot(orb3b, type = "orbits")  # P-value distribution
plot(orb3b, type = "cumulative")    # P-values by total score
plot(orb3b, type = "orbit", total_score = 10)  # Specific orbit
plot(orb3b, type = "significant")    # P-values by total score


# Get misfitting persons
misfits <- get_misfitting_persons(orbits)
# equating
leunbach_equate(fit3b, direction = "1to2", verbose = TRUE)

# bootstrap
boot3b <- leunbach_bootstrap(fit3b, nsim = 1000, verbose = TRUE, seed = 234, n_cores = 12)
boot3b

# testing d3b -------------------------------------------------------------



fit3c <- leunbach_ipf(d3c, verbose = TRUE)
fit3c
summary(fit3c)
# orbits
orb3c <- analyze_orbits(fit3c)
orb3c
summary(orb3c)
# Plots
plot(orb3b, type = "orbits")  # P-value distribution
plot(orb3b, type = "cumulative")    # P-values by total score
plot(orb3b, type = "orbit", total_score = 10)  # Specific orbit
plot(orb3b, type = "significant")    # P-values by total score

# equating
leunbach_equate(fit3c, direction = "1to2", verbose = TRUE)

# bootstrap
boot3c <- leunbach_bootstrap(fit3c, nsim = 1000, verbose = TRUE, seed = 234, n_cores = 12)
boot3b


# orbits ------------------------------------------------------------------

# Plots
plot(orbits, type = "histogram")  # P-value distribution
plot(orbits, type = "scatter")    # P-values by total score
plot(orbits, type = "orbit", total_score = 10)  # Specific orbit

# Get misfitting persons
misfits <- get_misfitting_persons(orbits)



# kequate -----------------------------------------------------------------


# number of moments chosen based on AIC
glm_a5 <- glm(count~I(total) + I(total^2) + I(total^3) + I(total^4) + I(total^5), 
              family = "poisson", data = rx, x = TRUE)
glm_b3 <- glm(count~I(total) + I(total^2) + I(total^3), 
              family = "poisson", data = ry, x = TRUE)
glm_a_ns <- glm(count ~ splines::bs(total, df = 3), 
               family = "poisson", data = rx, x = TRUE)
glm_a_ns_nb <- MASS::glm.nb(count ~ splines::bs(total, df = 3), 
              data = rx, x = TRUE, control = glm.control(maxit = 1000))
AIC(glm_a5,glm_a_ns,glm_a_ns_nb)
library(kequate)
# Equivalent groups design: In an EG design, two groups that have been randomly selected from a common population are given separate but parallel tests.
eg_eq <- kequate(design = "EG", # assuming that these are two independent groups from the same sample
        r = glm_a_ns, s = glm_b3,
        x = c(0:10), y = c(0:10))

eg_eq@equating
lboot2 <- leunbach_bootstrap(lfit, n_cores = 4, nsim = 100, see_type = "expected")


data.frame(score = c(0:10,0:10),
           see = c(lboot2[["see_1to2"]],eg_eq@equating$SEEYx),
           model = c(rep("Leunbach",11),rep("glm",11))
) %>% ggplot(aes(x=score,y=see, color = model)) +
  geom_point(size = 3) +
  geom_line()

#After having obtained satisfactory univariate models, the bivariate model needs to be estimated. 
# The bivariate model should include the univariate moments of the respective univariate log-linear models, along with cross-moments between the two tests X and A. We consider for inclusion in the bivariate log-linear model the cross-moments up to the highest moments contained in the univariate log-linear models.


rxx <- rx %>% 
  as.data.frame() %>% 
  rename(count_a = count) 
ryy <- ry %>% 
  as.data.frame() %>% 
  rename(count_b = count) 
# rxy <- rxx %>% 
#   add_column(count_b = ryy$count_b)

rxy <- data.frame(
  count = c(rxx$count_a,ryy$count_b),
  test_a = c(0:10,rep(NA,11)),
  test_b = c(rep(NA,11),0:10)
)

glm(count~I(test_a) + I(test_a^2) + I(test_a^3) + I(test_a^4) + I(test_a^5) + I(test_b) + I(test_b^2) + I(test_b^3) + I(test_a):I(test_b) + I(test_a^2):I(test_b^2), data = rxy, family = "poisson", x = TRUE)



## Indirect equating

#This is work in progress.


data1 <- read.delim("data/data1_item.csv", sep = ",") %>% 
  select(a01:a10,b01:b10,c01:c10)
data1_sum <- read.delim("data/data1_item.csv", sep = ",") %>% 
  select(a_sum,b_sum,c_sum)

# equate test a and b
lfit_ab <- leunbach_ipf(data1_sum[,c(1,2)])
# equate test b and c
lfit_bc <- leunbach_ipf(data1_sum[,c(2,3)])

indirect1 <- leunbach_indirect_equate(lfit_ab, lfit_bc,
                                      direction_ab = "1to2",
                                      direction_bc = "1to2")

# subset equating table for test a -> c
id1table <- indirect1[["equating_table"]]
id1table




