############################Proposal Cauchy sidtribution ################################
  
# To draw sample from N(0,1)

MH_draw <- function(n, initial_value, tuning_parameter){  # n : number of samples to be selected
                                     # initial_value : initial value of a possible sample from N(0,1)
count <- 0  # Number of samples selected
i <- 1
accept_rate <- 0


final_sample <- numeric(length = length(n)) #storing the samples

# Setting initial value of a possible sample from N(0,1)
theta_ini <- initial_value

while(count != n){

# Considering Cauchy proposal
 theta_current <- rcauchy(1, theta_ini, tuning_parameter)

# Computing acceptance probability
 acceptance_prob <- (dnorm(theta_current)*dcauchy(theta_ini, theta_current, tuning_parameter))/(dnorm(theta_ini)*dcauchy(theta_current, theta_ini, tuning_parameter))

# Decision whether to accept the proposed sample value or not
 if(log(runif(1)) < log(acceptance_prob)){
   
   final_sample[i] <- theta_current
   
   theta_ini <- theta_current
   accept_rate <- accept_rate + 1
 }else{
   
   final_sample[i] <- theta_ini
 
 }
 i <- i + 1
 count <- count + 1

}

 return(list(final_sample, accept_rate/n))

}


set.seed(1234)
# initial value 0.001


samples <- MH_draw(1000, 0.001, 1)
# Histogram of 10000 samples drawn from N(0,1) by MH Algorithm
hist(samples[[1]], main = "Histogram of 10000 Samples Drawn from N(0,1) Using  Metropolis-Hastings",
     probability = TRUE, xlab = "Sample Values", ylab = "Frequency Density",
     xlim = c(-4,4), ylim = c(0,0.5), col = "yellow")

plot.ts(samples[[1]])
acf(samples[[1]], lag.max = 50)

samples[[2]]  # acceptance probability 0.518  

###############################################################################33
#setting tuning parameter
tuning_parameter <- seq(0.01,10,0.01)
acc_prob <- array(0)


for(i in 1:length(tuning_parameter)){
  samples <- MH_draw(1000, 0.001, tuning_parameter[i])
  acc_prob[i] <- samples[[2]]
}
acc_prob[which.max(acc_prob)]
tuning_parameter[which.max(acc_prob)]

samples <- MH_draw(1000, 0.001, tuning_parameter[which.max(acc_prob)])
# Histogram of 10000 samples drawn from N(0,1) by MH Algorithm
hist(samples[[1]], main = "Histogram of 10000 Samples Drawn from N(0,1) Using  Metropolis-Hastings",
     probability = TRUE, xlab = "Sample Values", ylab = "Frequency Density",
     xlim = c(-4,4), ylim = c(0,0.5), col = "yellow")

plot.ts(samples[[1]])
acf(samples[[1]], lag.max = 50)
samples[[2]]
###################################################################################




###############################################################################33
#setting tuning parameter
tuning_parameter <- seq(0.01,10,0.01)
acc_prob <- array(0)


for(i in 1:length(tuning_parameter)){
  samples <- MH_draw(1000, 10, tuning_parameter[i])
  acc_prob[i] <- samples[[2]]
}
acc_prob[which.max(acc_prob)]
tuning_parameter[which.max(acc_prob)]

samples <- MH_draw(1000, 10, tuning_parameter[which.max(acc_prob)])
# Histogram of 10000 samples drawn from N(0,1) by MH Algorithm
hist(samples[[1]], main = "Histogram of 10000 Samples Drawn from N(0,1) Using  Metropolis-Hastings",
     probability = TRUE, xlab = "Sample Values", ylab = "Frequency Density",
     xlim = c(-4,4), ylim = c(0,0.5), col = "yellow")

plot.ts(samples[[1]])
acf(samples[[1]], lag.max = 50)
samples[[2]]
###################################################################################


###############################################################################33
#setting tuning parameter
tuning_parameter <- seq(0.01,10,0.01)
acc_prob <- array(0)


for(i in 1:length(tuning_parameter)){
  samples <- MH_draw(1000, -15, tuning_parameter[i])
  acc_prob[i] <- samples[[2]]
}
acc_prob[which.max(acc_prob)]
tuning_parameter[which.max(acc_prob)]

samples <- MH_draw(1000, -15, tuning_parameter[which.max(acc_prob)])
# Histogram of 10000 samples drawn from N(0,1) by MH Algorithm
hist(samples[[1]], main = "Histogram of 10000 Samples Drawn from N(0,1) Using  Metropolis-Hastings",
     probability = TRUE, xlab = "Sample Values", ylab = "Frequency Density",
      ylim = c(0,0.5), col = "yellow")

plot.ts(samples[[1]])
acf(samples[[1]], lag.max = 50)
samples[[2]]
###################################################################################























###################### Proposal Dexp sidtribution ############################33
library(nimble)
# To draw sample from N(0,1)

MH_draw <- function(n, initial_value, tuning_parameter){  # n : number of samples to be selected
  # initial_value : initial value of a possible sample from N(0,1)
  count <- 0  # Number of samples selected
  i <- 1
  accept_rate <- 0
  
  
  final_sample <- numeric(length = length(n)) #storing the samples
  
  # Setting initial value of a possible sample from N(0,1)
  theta_ini <- initial_value
  
  while(count != n){
    
    # Considering Cauchy proposal
    theta_current <- rdexp(1, location = theta_ini, scale = tuning_parameter)
    
    # Computing acceptance probability
    acceptance_prob <- (dnorm(theta_current)*ddexp(theta_ini,location = theta_current, scale = tuning_parameter))/(dnorm(theta_ini)*ddexp(theta_current, location = theta_ini, scale = tuning_parameter))
    
    # Decision whether to accept the proposed sample value or not
    if(log(runif(1)) < log(acceptance_prob)){
      
      final_sample[i] <- theta_current
      
      theta_ini <- theta_current
      accept_rate <- accept_rate + 1
    }else{
      
      final_sample[i] <- theta_ini
      
    }
    i <- i + 1
    count <- count + 1
    
  }
  
  return(list(final_sample, accept_rate/n))
  
}


set.seed(1234)
# initial value 0.001
samples <- MH_draw(1000,0.001, 1)
# Histogram of 10000 samples drawn from N(0,1) by MH Algorithm
hist(samples[[1]], main = "Histogram of 10000 Samples Drawn from N(0,1) Using  Metropolis-Hastings",
     probability = TRUE, xlab = "Sample Values", ylab = "Frequency Density",
     xlim = c(-4,4), ylim = c(0,0.5), col = "yellow")

plot.ts(samples[[1]])
acf(samples[[1]], lag.max = 50)

samples[[2]]  # acceptance probability 0.672  


###############################################################################33
#setting tuning parameter
tuning_parameter <- seq(0.01,10,0.01)
acc_prob <- array(0)


for(i in 1:length(tuning_parameter)){
  samples <- MH_draw(1000, 0.001, tuning_parameter[i])
  acc_prob[i] <- samples[[2]]
}
acc_prob[which.max(acc_prob)]
tuning_parameter[which.max(acc_prob)]

samples <- MH_draw(1000, 0.001, tuning_parameter[which.max(acc_prob)])
# Histogram of 10000 samples drawn from N(0,1) by MH Algorithm
hist(samples[[1]], main = "Histogram of 10000 Samples Drawn from N(0,1) Using  Metropolis-Hastings",
     probability = TRUE, xlab = "Sample Values", ylab = "Frequency Density",
     col = "yellow")

plot.ts(samples[[1]])
acf(samples[[1]], lag.max = 50)
samples[[2]]
###################################################################################






















###############################################################################33
#setting tuning parameter
tuning_parameter <- seq(0.01,10,0.01)
acc_prob <- array(0)


for(i in 1:length(tuning_parameter)){
  samples <- MH_draw(1000, 10, tuning_parameter[i])
  acc_prob[i] <- samples[[2]]
}
acc_prob[which.max(acc_prob)]
tuning_parameter[which.max(acc_prob)]

samples <- MH_draw(1000, 10, tuning_parameter[which.max(acc_prob)])
# Histogram of 10000 samples drawn from N(0,1) by MH Algorithm
hist(samples[[1]], main = "Histogram of 10000 Samples Drawn from N(0,1) Using  Metropolis-Hastings",
     probability = TRUE, xlab = "Sample Values", ylab = "Frequency Density",
     col = "yellow")

plot.ts(samples[[1]])
acf(samples[[1]], lag.max = 50)
samples[[2]]
###################################################################################

###############################################################################33
#setting tuning parameter
tuning_parameter <- seq(0.01,10,0.01)
acc_prob <- array(0)


for(i in 1:length(tuning_parameter)){
  samples <- MH_draw(1000, -15, tuning_parameter[i])
  acc_prob[i] <- samples[[2]]
}
acc_prob[which.max(acc_prob)]
tuning_parameter[which.max(acc_prob)]

samples <- MH_draw(1000, -15, tuning_parameter[which.max(acc_prob)])
# Histogram of 10000 samples drawn from N(0,1) by MH Algorithm
hist(samples[[1]], main = "Histogram of 10000 Samples Drawn from N(0,1) Using  Metropolis-Hastings",
     probability = TRUE, xlab = "Sample Values", ylab = "Frequency Density",
     col = "yellow")

plot.ts(samples[[1]])
acf(samples[[1]], lag.max = 50)
samples[[2]]
###################################################################################








