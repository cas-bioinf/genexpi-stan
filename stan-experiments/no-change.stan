data {
  int<lower=0> numData; // number data points
  //  real<lower=0> y[numData]; // measured values
  real y[numData]; // measured values
  real<lower=0> sigma[numData]; // s.e. of measured values
}
parameters {
  
  //real<lower = 0> mean_value;
  real true_value; 
  real<lower=0> absoluteError;
  real<lower=0> relativeError;
  real<lower=0> errorSpread;
}

model {
  //target += gamma_lpdf(true_value | 2, 2);
  real error = true_value * relativeError + absoluteError;
  target += cauchy_lpdf(sigma | error, errorSpread);
  target += normal_lpdf(y | true_value, error);
  /*
    target += normal_lpdf(y[1] | change[1], sigma[1]);
    for(j in 2:numData) {
      target += normal_lpdf(y[j] | y[j-1] + change[j], sigma[j]);
      target += normal_lpdf(y[j] | y[j-1] + change[j], sigma[j]);
    }*/
      //target += normal_lpdf(y | theta, sigma);
      
}

generated quantities {
  vector[numData] log_lik;
  for (n in 1:numData) log_lik[n] = normal_lpdf(y[n] | true_value, sigma[n]);
  
}
