data {
  int N; //Number of observations
  int<lower=0, upper=10> Y[N]; //The answers
  int<lower=1,upper=2> group[N]; //1 - pre treatment, 2 - post treatment
}

parameters {
  //Mean satisfaction per group
  real<lower=0, upper=1> a[2];

  //True satisfaction of each person
  real<lower=0, upper=1> b[N];

  //Thresholds
  positive_ordered[10] c;

  //'b' gives the same amount of information about 'a' as 'tau' coin tosses
  //with probability 'a'
  real<lower=0> tau;
}

model {
  for(i in 1:N) {
    real expected = a[group[i]];
    b[i] ~ beta(tau * expected, tau * (1 - expected));
    Y[i] ~ ordered_logistic(b[i], c);
  }

  //Priors
  //We assume both a to be more likely close to 0.5, but with a lot of leeway
  a ~ beta(2,2);
  //We assume that 'b' give roughly between (e^0) = 1 and (e^2) = 7.5
  //coin tosses worth of information about 'a'
  tau ~ lognormal(1,0.5);
}
