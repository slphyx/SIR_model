#
# Slphyx@Shift-Enter 
# 

library(rstan)

set.seed(12345)

STcode <-"
functions {
  real[] SIR(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real Susceptibles = y[1];
      real Infected = y[2];
      real Recovered = y[3];
      
      real beta = theta[1];
      real gamma = theta[2];
      
      real dSusceptibles = -beta*Susceptibles*Infected;
      real dInfected =  beta*Susceptibles*Infected-gamma*Infected;
      real dRecovered =  gamma*Infected;
      
      return {dSusceptibles, dInfected, dRecovered};
  }
}
data {
  int<lower = 1> n_obs;   // number of days observed
  real   y[n_obs];    // data, total number of infected individuals
  real   t0;                // initial time point 
  real   ts[n_obs];         // time points observed
  real  y_init[3];  //initial value of each compartment at t0
}

transformed data {
  real x_r[0];
  int x_i[0];
  int n_states = 3;
}

parameters {
  real<lower = 0> beta;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
}

transformed parameters{
  real<lower = 0> theta[2];
  theta[1] = beta;
  theta[2] = gamma;

  // ODE solutions
  real   y_hat[n_obs, n_states]; 
  
  real infected[n_obs]; 
  
  // Solve SIR model
  y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);
  
  // Extract I(t) 
  for (i in 1:n_obs) {
    infected[i] = y_hat[i, 2];
  }

}

model {
  
  //prior distributions
  beta ~ uniform(0.0001, 0.1);
  gamma ~ uniform(0.0001, 0.5);
  sigma ~ uniform(0.001,100);

  y ~ normal(infected,sigma)T[0,];
}

generated quantities {
  real y_mod[n_obs]; 
  
  // sampling the observed data from the posterior 
  for (i in 1:n_obs) {
    y_mod[i] = normal_rng(infected[i], sigma);  
  }
}

";


# this data was generated from solving the SIR model with beta=0.001 and 
# gamma=0.1. The I(t) compartment was added with the noises from  ~N(0,sd = 8)  
obs <- data.frame(times = 1:60, infected = c(6.44098, 13.6343, 0, 14.5879, 0, 14.1143, 22.4727, 22.7972, 35.5923, 
               55.7897, 70.4286, 86.1984, 114.689, 131.202, 171.406, 189.654, 
               222.74, 217.166, 227.565, 228.92, 246.056, 216.92, 211.376, 193.15, 
               186.952, 184.023, 146.325, 141.609, 126.664, 121.139, 102.171, 
               107.414, 106.47, 80.2677, 90.7381, 72.1889, 60.5055, 60.8584, 
               39.8837, 55.4039, 50.5164, 32.5718, 32.4639, 24.2174, 38.2163, 
               40.1582, 25.0519, 31.0153, 28.4255, 8.69976, 22.4933, 29.2557, 
               5.19587, 8.64702, 8.24683, 0, 19.7272, 18.3432, 9.25908, 9.47222))

plot(obs,xlab="time",ylab="infected")

stan_data <- list(
  n_obs = length(obs$times),   #number of observed data points
  t0 = 0,   # initial time
  ts = obs$times,  # output times for comparable with the data
  y = obs$infected,  #  data; number of infected cases
  y_init = c(499,1,0)  # initial value at time = t0
)

# initial function for giving the initial values of the parameters
init_fun <- function(...) list(beta=0.0005,gamma=0.2)

# run the mcmc
fit <- stan(model_code = STcode, data = stan_data, iter = 10000, chains = 2,init = init_fun)

# show the trace plot and the values
print(fit,pars = c("beta","gamma","sigma"),digits=5)
stan_trace(fit, pars = c("beta","gamma","sigma"))
stan_hist(fit,pars = c("beta","gamma","sigma"))


# plot the fitting results with the data
y_mod_samples <- extract(fit)$y_mod
y_mod_mean <- apply(y_mod_samples, 2, mean)
y_mod_ci <- apply(y_mod_samples, 2, quantile, probs = c(0.025, 0.975))

plot(obs,xlab="time",ylab="infected")
lines(x=1:60,y=y_mod_mean, col='red')
lines(x=1:60,y=y_mod_ci[1,])
lines(x=1:60,y=y_mod_ci[2,])

polygon(c(1:60, 60:1), c(y_mod_ci[1,], rev(y_mod_ci[2,])), 
        col=rgb(0.8, 0.8, 0.8, 0.5), border=NA)
