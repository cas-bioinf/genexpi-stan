generateConstantSynthesis <- function(sco, time, repeats, init, synthesis, decay, relativeNoise = 0.2, minNoise = 0.5)
{
  derivativeAtZero = synthesis - (decay * init);
  signDerivativeAtZero = sign(derivativeAtZero); 
  constantFactor = abs(derivativeAtZero) / (-decay);
  
  ratio = synthesis / decay;

  numValues = length(time) * repeats;
  
  times = numeric(numValues);
  replicates = numeric(numValues);
  values = numeric(numValues);
  
  index = 1;
  for(t in time)
  {
    true_value = signDerivativeAtZero * (constantFactor * exp(-decay * t))  + ratio;
    noise = max(c(minNoise, (true_value * relativeNoise)));
    for(r in 1:repeats) {
      times[index] = t;
      replicates[index] = r;
      values[index] = rnorm(1, true_value, noise);
      if(values[index] < 0){
        values[index] = runif(1,0.0001, true_value + noise);
      }
      index = index + 1;
    }
  }

  result = data.frame("sco"= rep.int(sco, numValues),"spot" = rep.int(1, numValues), "time" = times, "replicate" = replicates, "value" = values);
  return (result) ;
}