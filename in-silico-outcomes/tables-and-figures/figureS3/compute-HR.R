library(survival)
library(boot)
library(JinkoStats)

# (This function is copied from SIRIUS data analysis script)
# Estimation of Cox models based on survival.bootstrap outputs
# .. This function takes as input two outputs (control and treated) of survival.bootstrap
# .. ("SimulatedDataAllSamples.dtf" element of the function's output object) with an equal
# .. number of bootstrap iterations. It then estimates a Cox model for each pair of samples, and extracts
# .. the HR. It also checks if the proportional hazards assumption is met for each cox model.
BootstrapCox.fun <- function(BootDataCtrl.dtf,
                             BootDataTrtd.dtf,
                             Year.int = NULL # time of follow-up
){
  # First checking if both input data-frames have the same boot iterations
  if(!all(length(unique(BootDataCtrl.dtf[, "study"])) == length(unique(BootDataTrtd.dtf[, "study"]))) |
     !all(unique(BootDataCtrl.dtf[, "study"]) == unique(BootDataTrtd.dtf[, "study"]))){
    stop("The bootstrap iteration available in both datasets seem to differ, check your inputs.")
  }

  # If "Year.int" is defined, then events happening after Year.int are censored
  if(!is.null(Year.int)){
    BootDataCtrl.dtf$status <- ifelse(test = BootDataCtrl.dtf$time > Year.int, yes = 0, no = BootDataCtrl.dtf$status)
    BootDataTrtd.dtf$status <- ifelse(test = BootDataTrtd.dtf$time > Year.int, yes = 0, no = BootDataTrtd.dtf$status)
  }

  # If "Year.int" is defined, then all times are caped at Year.int
  if(!is.null(Year.int)){
    BootDataCtrl.dtf$time <- ifelse(test = BootDataCtrl.dtf$time > Year.int, yes = Year.int, no = BootDataCtrl.dtf$time)
    BootDataTrtd.dtf$time <- ifelse(test = BootDataTrtd.dtf$time > Year.int, yes = Year.int, no = BootDataTrtd.dtf$time)
  }

  # Creating a data-frame that will store the function's outputs
  Outputs.dtf <- data.frame()

  # Iterating on boot samples
  for (bootIter in unique(unique(BootDataCtrl.dtf[, "study"]))){
    # Selecting the subsets corresponding to the boot iteration of interest
    TempCtrl.dtf <- BootDataCtrl.dtf[BootDataCtrl.dtf[, "study"] == bootIter, c("time", "status")]
    TempTrtd.dtf <- BootDataTrtd.dtf[BootDataTrtd.dtf[, "study"] == bootIter, c("time", "status")]

    # Adding a column to differentiate the control from the treated arm (control = 0, treated = 1)
    TempCtrl.dtf$arm <- 0
    TempTrtd.dtf$arm <- 1
    # Binding both dataframes
    TempDataForCox.dtf <- rbind(TempCtrl.dtf, TempTrtd.dtf)

    # Estimating the Cox model
    TempCox.mdl <- coxph(Surv(time, status) ~ arm, data = TempDataForCox.dtf)
    # Checking if prop. hazards assumption is met
    pValue.dbl <- tryCatch(
      {
        round(survival::cox.zph(fit = TempCox.mdl)$table["GLOBAL", "p"], digits = 3)
      },
      error=function(cond){
        return(NA)
      }
    )



    # Adding the values to the output data-frame
    Outputs.dtf <- rbind(Outputs.dtf, c(bootIter,                                         # Boot iteration
                                        round(x = exp(coef(TempCox.mdl)[1]), digits = 3), # HR
                                        round(x = exp(confint(TempCox.mdl)), digits = 3), # 95% CI
                                        pValue.dbl                                        # p-value of the prop. haz. test
                                        ))

  }
  names(Outputs.dtf) <- c("BootIter", "HR", "HR95LoCI", "HR95HiCI", "PropHaz")

  Outputs.dtf <- Outputs.dtf %>% dplyr::mutate_at(c('HR', 'HR95LoCI', 'HR95HiCI', 'PropHaz'), as.numeric)

  return(Outputs.dtf)
}

# A function to run bootstrap and extract the HR and 95% PPI of HR (arm2 vs arm1)
BootstrapHRAndPPI.fun <- function(SimuData.dtf, ArmOne.var, ArmTwo.var,
  TTE.var, Outcome.var,
  nb.boot, nb.clusters, pct.sample) {

  ArmOne.dtf <- data.frame(cbind(
    time = SimData.dtf[SimData.dtf$Arm == ArmOne.var, TTE.var],
    status = SimData.dtf[SimData.dtf$Arm == ArmOne.var, Outcome.var],
    arm = 1))

  ArmTwo.dtf <- data.frame(cbind(
    time = SimData.dtf[SimData.dtf$Arm == ArmTwo.var, TTE.var],
    status = SimData.dtf[SimData.dtf$Arm == ArmTwo.var, Outcome.var],
    arm = 2))

  # One realization based on Cox model estimation
  Cox.mdl        <- coxph(
    Surv(time, event = status) ~ arm,
    data = rbind(ArmOne.dtf, ArmTwo.dtf))
  HRtrial0.vec   <- exp(coef(Cox.mdl)[1])

  # Hazard ratio prediction intervals computation
  SurvArmOne.list <- Survival.Bootstrap(
    SampleSize   = dim(ArmOne.dtf)[1]*pct.sample,
    NbBootstraps = nb.boot,
    SurvivalData = ArmOne.dtf,
    NbClusters   = nb.clusters)

  SurvArmTwo.list <- Survival.Bootstrap(
    SampleSize   = dim(ArmTwo.dtf)[1]*pct.sample,
    NbBootstraps = nb.boot,
    SurvivalData = ArmTwo.dtf,
    NbClusters   = nb.clusters)

  resCox.dtf <- BootstrapCox.fun(BootDataCtrl.dtf = SurvArmOne.list$SimulatedDataAllSamples.dtf,
                                 BootDataTrtd.dtf = SurvArmTwo.list$SimulatedDataAllSamples.dtf)

  resCox.dtf <- resCox.dtf %>% mutate_at(c('HR', 'HR95LoCI', 'HR95HiCI'), as.numeric)
  width.ppi <- mean(resCox.dtf$HR95HiCI, na.rm = T) - mean(resCox.dtf$HR95LoCI, na.rm = T)
  sprintf('%s: %.3f (%.3f-%.3f) Uncertainty: %s',
    Outcome.var,
    HRtrial0.vec,
    mean(resCox.dtf$HR95LoCI, na.rm = T),
    mean(resCox.dtf$HR95HiCI, na.rm = T),
    ifelse(width.ppi < 0.2, "Low", ifelse(width.ppi < 0.5, "Medium", "High"))
    )
}

AllSimData.dtf <- read.csv('scalars-trimmed-secondary-n2.csv') %>%
  select(Arm, isMACE, timeOfMACE, baselineLDLcQuartile)

SimData.dtf <- AllSimData.dtf %>% filter(baselineLDLcQuartile == 3)

summary.events <- SimData.dtf %>% group_by(Arm) %>%
  select(Arm, isMACE) %>%
  summarize(mean(isMACE)*100)
t(summary.events)

set.seed(2023)

ArmOne.var <- "RWLLT"
ArmTwo.var <- "RWLLT-ICL"
nb.boot <- 100
nb.clusters <- 12
pct.sample <- 13780/204691/4 # size of FOURIER divided by 4


BootstrapHRAndPPI.fun(SimuData.dtf, ArmOne.var, ArmTwo.var,
  "timeOfMACE", "isMACE",
  nb.boot, nb.clusters, pct.sample)



