/*
The computation in this model seems to be more problematic than the -derivative and -euler variant. 
In particular, synthesis and decay are happy to blow out of proportions
*/

data {
  int<lower=0> numData; // number data points
  //  real<lower=0> y[numData]; // measured values
  real y[numData]; // measured values
  real<lower=0> sigma[numData]; // s.e. of measured values
  real time[numData]; //time of measurements
}
parameters {
  
  //real<lower = 0> mean_value;
  real init; 
  real<lower = 0> synthesis;
  real<lower = 0> decay;
}

transformed parameters {
  real true_value[numData];
  
  {

    real derivativeAtZero = synthesis - (decay * init);
    real signDerivativeAtZero = (derivativeAtZero > 0) ? 1 : -1; //sign of the derivative
    real constantFactor = fabs(derivativeAtZero) / (-decay);

    real ratio = synthesis / decay;
    //print("Synthesis: ", synthesis, " Decay: ", decay, " Derivative at 0:", derivativeAtZero, " constantFactor: ", constantFactor, "\n");
    
    for(j in 1:numData)
    {
      true_value[j] = signDerivativeAtZero * (constantFactor * exp(-decay * time[j])) + ratio;
      //print("Value at ", time[j], " = ", true_value[j]);
    }
  }
}

model {
  target += normal_lpdf(y | true_value, sigma);
}

generated quantities {
  vector[numData] log_lik;
  for (n in 1:numData) log_lik[n] = normal_lpdf(y[n] | true_value[n], sigma[n]);
  
}
