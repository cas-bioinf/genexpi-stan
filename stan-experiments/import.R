#rawData = read.csv("C:\\Users/MBU/Documents/Genexpi/genenetworks/data/SigQ/Chip ratios normalized repeats.csv");
require(plyr)

chipRatiosNormalizedToLong <- function(rawData)
{
  sco = rawData$sco;
  numData =  length(sco);
  
  timeColumns = c("Timeneg", vapply(seq(from = 0, to = 5.5, by = 0.5), FUN = function(x) { paste0("Time", x) }, FUN.VALUE=""));
  times = c(-0.5, seq(from = 0, to = 5.5, by = 0.5)) * 60;
  
  result = data.frame("sco"= c(),"spot" = c(), "time" = c(), "replicate" = c(), "value" = c());
  for(replicate in 1:3)
  {
    for(timeIdx in 1:length(times)) 
    {
      columnName = paste(timeColumns[timeIdx],replicate,sep="_");
      if(columnName %in% names(rawData)) 
      {
        newRows = data.frame("sco"= sco, "spot" = rawData$spot, "time" = rep(times[timeIdx], numData), "replicate" = rep(replicate,numData), "value" = rawData[, columnName]);
        result = rbind(result, newRows);
      }
      
    }
  }
  
  return(result)
}

longToVal <- function(longData)
{
  ddply(longData, .(sco,time), summarise, val = mean(value), sd = sd(value))
}

#streptomyces_long = chipRatiosNormalizedToLong(rawData);
#streptomyces_val = longToVal(chipRatiosNormalizedToLong(rawData))