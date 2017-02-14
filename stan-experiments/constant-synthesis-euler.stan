data {
  int<lower=0> numData; // number data points
  //  real<lower=0> y[numData]; // measured values
  real<lower=0> y[numData]; // measured values
  real<lower=0> sigma[numData]; // s.e. of measured values
  real time[numData]; //time of measurements
}
parameters {
  
  //real<lower = 0> mean_value;
  real<lower = 0> init; 
  real<lower = 0> synthesis;
  real<lower = 0> decay;
}

transformed parameters {
  real<lower = 0> true_value[numData];
  true_value[1] = init;
  for(j in 2:numData) {
    real derivative = synthesis - decay * true_value[j-1];
    true_value[j] = true_value[j - 1] + (time[j] - time[j - 1]) * derivative;
  } 
}

model {

    target += normal_lpdf(y | true_value, sigma);
}

generated quantities {
  vector[numData] log_lik;
  for (n in 1:numData) log_lik[n] = normal_lpdf(y[n] | true_value[n], sigma[n]);
}
