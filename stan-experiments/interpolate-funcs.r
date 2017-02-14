library(rstan)
library(loo)
rstan_options(auto_write = TRUE);
options(mc.cores = parallel::detectCores());

plotProfileFit <- function(fit, data, targetSco, numSamples = 10) {
  true_value = extract(fit,'true_value')$true_value;
  
  samplesToPlot = true_value[sample(1:(dim(true_value)[1]),numSamples),];
  
  #Add the raw profile
  plotData = t(rbind(subset(data, sco ==targetSco)[,"val"], samplesToPlot));
  defaultWidth = 1;
  lineWidths = rep.int(1, numSamples + 1);
  lineWidths[1] = defaultWidth * 2;
  matplot(plotData, lwd = lineWidths,  type="l") 
}

plotNoChangeFit <- function(fit, data, targetSco) {
  relevant_data = subset(data, sco==targetSco);
  minTime = min(relevant_data$time)
  maxTime = max(relevant_data$time)
  
  smmry = summary(fit)$summary;
  val_mean = smmry["true_value","mean"];
  val_25 = smmry["true_value","25%"];
  val_75 = smmry["true_value","75%"];
  val_2_5 = smmry["true_value","2.5%"];
  val_97_5 = smmry["true_value","97.5%"];
  
  ribbon1 = data.frame("x" = c(minTime,maxTime),"ymin" = c(val_2_5,val_2_5), "ymax" = c(val_97_5, val_97_5))
  ribbon2 = data.frame("x" = c(minTime,maxTime),"ymin" = c(val_25,val_25), "ymax" = c(val_75, val_75))
  ribbonAes = aes(x = x, ymin = ymin, ymax = ymax)
  
  ggplot(relevant_data) +
    geom_ribbon(mapping = ribbonAes, data = ribbon1, fill = "red", alpha = 0.3) +
    geom_ribbon(mapping = ribbonAes, data = ribbon2, fill = "darkred", alpha = 0.3) +
    geom_hline(yintercept = val_mean, color = "black", size = 2) +
    geom_ribbon(aes(x=time, ymin = val - sd, ymax = val + sd ), alpha = 0.3, fill = "blue") +
    geom_ribbon(aes(x=time, ymin = val - sd, ymax = val + sd ), alpha = 0.3, fill = "blue") +
    geom_line(aes(x=time, y=val))
}

plotSco <- function(data, targetSco) {
  relevant_data = subset(data, sco==targetSco);

  ggplot(relevant_data, aes(x=time, y=val, ymin = val - sd, ymax = val + sd )) + geom_ribbon(alpha = 0.3, fill = "blue") +  geom_line() 
}

fitBySco <- function(long_data, target_sco, model='interpolate.stan', ...) {
  relevant_data = subset(long_data, sco == target_sco);
  relevant_data = relevant_data[sort.list(relevant_data$time),];
  
  #normalize time
  relevant_data$time = (relevant_data$time - min(relevant_data$time)) / (max(relevant_data$time) - min(relevant_data$time));
  
  data = list(numData = length(relevant_data$sco), y = relevant_data[,"val"], sigma = relevant_data[,"sd"], time = relevant_data[,"time"]);
  return(stan(file =model, data = data, ...));
}

compareFits <- function(data, target_sco) {
  ow <- options("warn");
  options(warn = 1);
  
  cat("No change fit\n");
  noChangeFit = fitBySco(data, target_sco, 'no-change.stan');
  print(plotNoChangeFit(noChangeFit, data, target_sco));
  
  cat("C. synth fit\n");
  csynthFit = fitBySco(data, target_sco, 'constant-synthesis-euler.stan', control=list(adapt_delta = 0.98), iter = 5000);
  print(plotProfileFit(csynthFit, data, target_sco));
  
  cat("Compare > 0 means csynth better");
  
  options(ow);
  compare(loo(extract_log_lik(noChangeFit)), loo(extract_log_lik(csynthFit)));
}