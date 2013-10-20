Importance of modelling variance in population growth rates in accurately predicting future population size.
========================================================

This code is based of Morris & Doak's MATLAB code (CITATION, pg. 29-30).  This code does multiple simulations of discrete exponential gorwht trajectories with a set of observed population growth rates (lambda values).


```r
########## Simulation Parameters

lambdas <- c(1, 1.98, 1.02, 0.92, 0.53)  # different lambda values to use for simulation.  There can be as many as you want.  For this simulation, each will be drawn with an equal probability.

start_pop_size <- 29  # starting population size
tmax <- 100  # length of simulations
numreps <- 100  # number of replicate trajectories to simulate
```


Below I have written a function to conduct population simulations

```r

count_PVA_clapper_simulation <- function(lambdas = lambdas, start_pop_size = 29, 
    tmax = 100, numreps = 100) {
    lambdas_simulation <- sample(x = lambdas, size = tmax, replace = TRUE)  # generates a random sample of lambdas over the entire simulation, sampling with replacement
    
    simulated_lambdas_matrix <- matrix(data = 0, ncol = numreps, nrow = tmax - 
        1)
    for (i in 1:numreps) {
        simulated_lambdas_matrix[, i] <- sample(x = lambdas, size = (tmax - 
            1), replace = TRUE)
    }
    
    pop_sizes_simulation <- matrix(data = c(rep(x = start_pop_size, times = numreps), 
        rep(x = 0, times = numreps * (tmax - 1))), ncol = numreps, byrow = TRUE)
    
    for (i in 1:numreps) {
        for (j in 1:(tmax - 1)) {
            pop_sizes_simulation[j + 1, i] <- simulated_lambdas_matrix[j, i] * 
                pop_sizes_simulation[j, i]
        }
    }
    
    pop_sizes_simulation[pop_sizes_simulation < 1] <- 0  # identified all instances where population size dipped below 1 individual and labeled it as extinct. 
    
    min_pop_size_each_rep_year5 <- vector()
    for (i in 1:numreps) {
        min_pop_size_each_rep_year5[i] <- min(pop_sizes_simulation[1:5, i])
    }
    pop_sizes_simulation[5, ][min_pop_size_each_rep_year5 == 0] <- 0
    pop_sizes_simulation[5, ][pop_sizes_simulation[5, ] > 200] <- 200
    year_5_pop_sizes <- as.data.frame(pop_sizes_simulation[5, ])
    year_5_hist <- ggplot(year_5_pop_sizes, aes(x = pop_sizes_simulation[5, 
        ])) + geom_histogram() + xlab("")
    
    min_pop_size_each_rep_year20 <- vector()
    for (i in 1:numreps) {
        min_pop_size_each_rep_year20[i] <- min(pop_sizes_simulation[1:20, i])
    }
    pop_sizes_simulation[20, ][min_pop_size_each_rep_year20 == 0] <- 0
    pop_sizes_simulation[20, ][pop_sizes_simulation[20, ] > 200] <- 200
    year_20_pop_sizes <- as.data.frame(pop_sizes_simulation[20, ])
    year_20_hist <- ggplot(year_20_pop_sizes, aes(x = pop_sizes_simulation[20, 
        ])) + geom_histogram() + xlab("")
    
    min_pop_size_each_rep_year100 <- vector()
    for (i in 1:numreps) {
        min_pop_size_each_rep_year100[i] <- min(pop_sizes_simulation[1:100, 
            i])
    }
    pop_sizes_simulation[100, ][min_pop_size_each_rep_year100 == 0] <- 0
    pop_sizes_simulation[100, ][pop_sizes_simulation[100, ] > 200] <- 200
    year_100_pop_sizes <- as.data.frame(pop_sizes_simulation[100, ])
    year_100_hist <- ggplot(year_100_pop_sizes, aes(x = pop_sizes_simulation[100, 
        ])) + geom_histogram() + xlab("Population size")
    
    source("~/Documents/miscellaneous_R/multiplot.R")
    multiplot(year_5_hist, year_20_hist, year_100_hist, cols = 1)
    
    # Geometric mean
    mean_geom <- lambdas[1]^0.2 * lambdas[2]^0.2 * lambdas[3]^0.2 * lambdas[4]^0.2 * 
        lambdas[5]^0.2  # predicts an almost stable population (lambda = 1)
    return(mean_geom)
}

```




```r
count_PVA_lambda_variance_simulation <- function(mean_lambda, sd_lambda, start_size = 1000, 
    tiempo = 50, reps = 1000) {
    sim_lambdas_matrix <- matrix(data = rnorm(n = reps * (tiempo - 1), mean = mean_lambda, 
        sd = sd_lambda), ncol = reps, nrow = tiempo - 1)
    
    pop_sizes_sim <- matrix(data = c(rep(x = start_size, times = reps), rep(x = 0, 
        times = reps * (tiempo - 1))), ncol = reps, byrow = TRUE)
    
    for (i in 1:reps) {
        for (j in 1:(tiempo - 1)) {
            pop_sizes_sim[j + 1, i] <- sim_lambdas_matrix[j, i] * pop_sizes_sim[j, 
                i]
        }
    }
    
    min_pop_size_each_rep_year50 <- vector()
    for (i in 1:reps) {
        min_pop_size_each_rep_year50[i] <- min(pop_sizes_sim[1:50, i])
    }
    pop_sizes_sim[50, ][min_pop_size_each_rep_year50 < 1] <- 0
    # pop_sizes_sim[50, ][pop_sizes_sim[50, ] < 1] <- 0
    pop_sizes_sim[50, ][pop_sizes_sim[50, ] > 10000] <- 10000
    
    hist(pop_sizes_sim[50, ], main = "Distribution of Population Sizes")
}
```


