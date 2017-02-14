testConstantSynthesis <- function(init, synthesis, decay, minNoise = 0.5, relativeNoise = 0.2, time = seq(0,1,by = 0.1))
{
  cat(paste0("Target init: ", init, " synth: ", synthesis, " decay: ", decay , "\n"));
  testConstant = longToVal(generateConstantSynthesis(1, time, 3, init, synthesis,decay, minNoise = minNoise, relativeNoise = relativeNoise));
  testConstantSynthesisModel(testConstant, 1, 'constant-synthesis.stan');
  testConstantSynthesisModel(testConstant, 1, 'constant-synthesis-euler.stan');
  testConstantSynthesisModel(testConstant, 1, 'constant-synthesis-derivative.stan', FALSE);
}

testConstantSynthesisModel <- function(testData, sco, model, plot = TRUE)
{
  cat(paste0("Model: ", model), "\n");
  fit = fitBySco(testData, sco, model, control=list(adapt_delta = 0.98), iter = 3000);
  if(plot)
  {
    plotProfileFit(fit, testData, sco);
  }
  
  smmry = summary(fit)$summary;
  init_mean = smmry["init","mean"];
  init_25 = smmry["init","25%"];
  init_75 = smmry["init","75%"];
  synthesis_mean = smmry["synthesis","mean"];
  synthesis_25 = smmry["synthesis","25%"];
  synthesis_75 = smmry["synthesis","75%"];
  decay_mean = smmry["decay","mean"];
  decay_25 = smmry["decay","25%"];
  decay_75 = smmry["decay","75%"];
  
  cat(paste0("Init: ", init_25, " -> ", init_mean, " <- ", init_75, "\n"));
  cat(paste0("Synthesis: ", synthesis_25, " -> ", synthesis_mean, " <- ", synthesis_75, "\n"));
  cat(paste0("Decay: ", decay_25, " -> ", decay_mean, " <- ", decay_75, "\n"));
  cat("\n");
}

testConstantSynthesisModelComparison <- function(init, synthesis, decay, minNoise = 0.5, relativeNoise = 0.2, time = seq(0,1,by = 0.1))
{
  sco = 1;
  testData = longToVal(generateConstantSynthesis(sco, time, 3, init, synthesis,decay, minNoise = minNoise, relativeNoise = relativeNoise));
  compareFits(testData, sco);
}