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


model {
    target += normal_lpdf(init | y[1], sigma[1]);
    for(j in 2:numData){
      real derivative = synthesis - decay * y[j-1];
      target += normal_lpdf(y[j] | y[j-1] + derivative * (time[j] - time[j-1]), sigma[j]);
    }
}
