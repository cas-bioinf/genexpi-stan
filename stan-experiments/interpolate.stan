data {
  int<lower=0> numData; // number data points
//  real<lower=0> y[numData]; // measured values
  real y[numData]; // measured values
  real<lower=0> sigma[numData]; // s.e. of measured values
}
parameters {

  //real<lower = 0> mean_value;
  real mean_change; 
  real<lower=0> change_dev;
  real<lower=0> true_value[numData];
}

transformed parameters {
  //real<lower=0> true_value[numData];
  real change[numData - 1];
  for (j in 2:numData)
    change[j - 1] = true_value[j] - true_value[j-1];
}

model {
  //TODO using y ~ normal(0,1) could be faster (but then I do not get correct log likelihood from the samples)
  
  //target += gamma_lpdf(true_value | 2, 2);
  target += cauchy_lpdf(change_dev | 0, 1);
  target += normal_lpdf(change | mean_change, change_dev);
  for(j in 1:numData) {
    target += normal_lpdf(y | true_value, sigma);
  }
  /*
  target += normal_lpdf(y[1] | change[1], sigma[1]);
  for(j in 2:numData) {
    target += normal_lpdf(y[j] | y[j-1] + change[j], sigma[j]);
    target += normal_lpdf(y[j] | y[j-1] + change[j], sigma[j]);
  }*/
  //target += normal_lpdf(y | theta, sigma);
  
}
