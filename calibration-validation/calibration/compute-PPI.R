library(boot)

## Compute bootsrap interval ---
## Copied from validation script with minor modification to work for tidied dataset

# Compute a "predicted percentile interval" with simulated data
# 1. take N samples of size n of my simulated measure
# 2. compute the statistic of the measure in each sample (e.g. mean)
# 3. get the 2.5 and 97.5 quantiles of the empirical distribution of the statistic
getBootCI <- function(simData.dtf, Measure.name, Arm.str, Time.name, Time.val,
  fun, sampleSize.int, n.boot) {
  simData.df <- simData.dtf[simData.dtf$Arm == Arm.str,]
  if (!is.na(Arm.str) & !is.na(Time.val) & nrow(simData.df) > 0)
    simData.df <- simData.df[simData.df[, Time.name] >= Time.val, ]
  print(paste("number of patients:", dim(simData.df)[1]))

  if (!(Measure.name %in% names(simData.df) & nrow(simData.df) > 0))
        return("NA")
    # Create a custom statistic function that uses the sample function
    custom_statistic <- function(data, indices) {
      sample_data <- data[sample(length(data), size = sampleSize.int, replace = TRUE)]
      result <- get(fun)(sample_data)
      return(result)
    }
    # Use the custom statistic function in the boot function
    res.boot <- boot(
      sim = "ordinary",
      data = simData.df[, Measure.name],
      statistic = custom_statistic,
      R = n.boot
    )

  ci <- boot.ci(res.boot, type = "perc")

  if(!is.null(ci)) {
    ci.df <- data.frame(
      Value = get(fun)(simData.df[, Measure.name]),
      LowRange = ci$percent[4],
      HighRange = ci$percent[5],
      RangeType = "95% PPI"
    )
  } else {
    val.mean = mean(simData.df[, Measure.name])
    ci.df <- data.frame(Value = val.mean,
      LowRange = val.mean, HighRange = val.mean, RangeType = "NA")
  }
  sprintf("%.1f [%.1f; %.1f]", ci.df$Value, ci.df$LowRange, ci.df$HighRange)
}


simData.df <- read.csv('vpop-lipo-data/scalars-trimmed-for-orion10-ldlc-tidied.csv')
nb.boot <- 100
set.seed(2024)

getBootCI(simData.df, "relChangeLDLcUncorrectedD90", "ORION10-inclisiran-300mg",
  "timeOfCvDeath", 0.25, "mean", sampleSize.int=781, n.boot=nb.boot)
getBootCI(simData.df, "relChangeLDLcUncorrectedD150", "ORION10-inclisiran-300mg",
  "timeOfCvDeath", 0.41, "mean", sampleSize.int=781, n.boot=nb.boot)
getBootCI(simData.df, "relChangeLDLcUncorrectedD270", "ORION10-inclisiran-300mg",
  "timeOfCvDeath", 0.74, "mean", sampleSize.int=781, n.boot=nb.boot)
getBootCI(simData.df, "relChangeLDLcUncorrectedD330", "ORION10-inclisiran-300mg",
  "timeOfCvDeath", 0.90, "mean", sampleSize.int=781, n.boot=nb.boot)
getBootCI(simData.df, "relChangeLDLcUncorrectedD450", "ORION10-inclisiran-300mg",
  "timeOfCvDeath", 1.23, "mean", sampleSize.int=781, n.boot=nb.boot)
getBootCI(simData.df, "relChangeLDLcUncorrectedD510", "ORION10-inclisiran-300mg",
  "timeOfCvDeath", 1.39, "mean", sampleSize.int=781, n.boot=nb.boot)
getBootCI(simData.df, "relChangeLDLcUncorrectedD540", "ORION10-inclisiran-300mg",
  "timeOfCvDeath", 1.47, "mean", sampleSize.int=781, n.boot=nb.boot)


simData.df <- read.csv('vpop-lipo-data/scalars-trimmed-for-fourier-ldlc-tidied.csv')
nb.boot <- 100
set.seed(2024)

getBootCI(simData.df, "baselineLDLcUncorrectedMass", "FOURIER-placebo",
  "timeOfCvDeath5Y", 0, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW4", "FOURIER-placebo",
  "timeOfCvDeath5Y", 0.07, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW12", "FOURIER-placebo",
  "timeOfCvDeath5Y", 0.23, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW24", "FOURIER-placebo",
  "timeOfCvDeath5Y", 0.46, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW48", "FOURIER-placebo",
  "timeOfCvDeath5Y", 0.92, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW72", "FOURIER-placebo",
  "timeOfCvDeath5Y", 1.38, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW96", "FOURIER-placebo",
  "timeOfCvDeath5Y", 1.84, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW120", "FOURIER-placebo",
  "timeOfCvDeath5Y", 2.30, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW144", "FOURIER-placebo",
  "timeOfCvDeath5Y", 2.76, "median", sampleSize.int=13780, n.boot=nb.boot)


getBootCI(simData.df, "baselineLDLcUncorrectedMass", "FOURIER-evolocumab",
  "timeOfCvDeath5Y", 0, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW4", "FOURIER-evolocumab",
  "timeOfCvDeath5Y", 0.07, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW12", "FOURIER-evolocumab",
  "timeOfCvDeath5Y", 0.23, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW24", "FOURIER-evolocumab",
  "timeOfCvDeath5Y", 0.46, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW48", "FOURIER-evolocumab",
  "timeOfCvDeath5Y", 0.92, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW72", "FOURIER-evolocumab",
  "timeOfCvDeath5Y", 1.38, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW96", "FOURIER-evolocumab",
  "timeOfCvDeath5Y", 1.84, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW120", "FOURIER-evolocumab",
  "timeOfCvDeath5Y", 2.30, "median", sampleSize.int=13780, n.boot=nb.boot)
getBootCI(simData.df, "valueLDLcUncorrectedMassW144", "FOURIER-evolocumab",
  "timeOfCvDeath5Y", 2.76, "median", sampleSize.int=13780, n.boot=nb.boot)

