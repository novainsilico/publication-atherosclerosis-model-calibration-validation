green3.color <- "#58AE7A"
blue2.color <- "#7ABBE8"
dark.blue.color <- "#202A44"
light.grey.color <- "#F5F9FC"
# ctrl.color <- "black"
ctrl.color <- blue2.color
trmt.color <- "orange"
ref.color <- "grey"



## Compute bootsrap interval ---

# Compute a "predicted percentile interval" with simulated data
# 1. take N samples of size n of my simulated measure
# 2. compute the statistic of the measure in each sample (e.g. mean)
# 3. get the 2.5 and 97.5 quantiles of the empirical distribution of the statistic
getBootCI <- function(simData.df, Measure.name, Condition.str, Arm.str, Time.val, fun, sampleSize.int) {
  # extract Condition
  if (!is.na(Condition.str)) simData.df <- applyCondition.fct(simData.df, Condition.str)

  if (!is.na(Arm.str) & !is.na(Time.val) & nrow(simData.df) > 0) {
    simData.df <- simData.df[simData.df[, paste0(Arm.str, "_timeOfCvDeath")] >= Time.val, ]
  }

  if (Measure.name %in% names(simData.df) & nrow(simData.df) > 0) {
      # Create a custom statistic function that uses the sample function
      custom_statistic <- function(data, indices) {
        # print("Indices list:");
        # print(indices);
        sample_data <- data[sample(length(data), size = sampleSize.int, replace = TRUE)]
        result <- get(fun)(sample_data)
        # print(paste("sample size = ", length(sample_data))); 
        # print("Used sample data:");
        # print(sample_data)
        return(result)
      }
      
      # Use the custom statistic function in the boot function
      res.boot <- boot(
        data = simData.df[, Measure.name],
        statistic = custom_statistic,
        R = nBoot
      )
      
    
    ci <- boot.ci(res.boot, type = "perc")

    if(!is.null(ci)){
      ci.df <- data.frame(
        Value = get(fun)(simData.df[, Measure.name]),
        LowRange = ci$percent[4],
        HighRange = ci$percent[5],
        RangeType = "95% PPI"
      )

    }else{
      val.mean = mean(simData.df[, Measure.name])
      ci.df <- data.frame(
        Value = val.mean,
        LowRange = val.mean,
        HighRange = val.mean,
        RangeType = "NaN"
      )
    }

    return(ci.df)
  } else {
    return(data.frame(Value = NA, LowRange = NA, HighRange = NA, RangeType = NA))
  }
}





# Compute bootstrapped inter-quartile range
getBootIQR <- function(simData.df, Measure.name, Condition.str, Arm.str, Time.val) {
  # extract Condition
  if (!is.na(Condition.str)) simData.df <- applyCondition.fct(simData.df, Condition.str)

  if (nrow(simData.df) == 0) {
    warning(paste(Condition.str, "not respected in the simulated dataset"))
    return(data.frame(Value = NA, LowRange = NA, HighRange = NA, RangeType = NA))
  }

  if (!is.na(Arm.str) & !is.na(Time.val) & nrow(simData.df) > 0) {
    simData.df <- simData.df[simData.df[, paste0(Arm.str, "_timeOfCvDeath")] >= Time.val, ]
  }

  if (Measure.name %in% names(simData.df)) {
    # iqr.ci <- expect_warning(ci_IQR(simData.df[, Measure.name], R = nBoot, boot_type = "perc"),"extreme order statistics used as endpoints") # if too much warnings
    iqr.ci <- ci_IQR(simData.df[, Measure.name], R = nBoot, boot_type = "perc")

    iqr.df <- data.frame(
      Value = median(simData.df[, Measure.name]),
      LowRange = median(simData.df[, Measure.name]) - iqr.ci$estimate / 2,
      HighRange = median(simData.df[, Measure.name]) + iqr.ci$estimate / 2,
      RangeType = "Boot-IQR"
    )

    return(iqr.df)
  } else {
    print(paste0(Measure.name, "not in the simulated dataset"))
    return(data.frame(Value = NA, LowRange = NA, HighRange = NA, RangeType = NA))
  }
}


# Parse a subgroup condition ----
# Each condition is defined as "Baseline == X"
# In it's current version, conditions can only be connected with  `And`
# Apply the subgroup condition to a simulated dataset
applyCondition.fct <- function(simData.df, Condition.str) {
  if (stringr::str_length(Condition.str) > 0 | Condition.str == "NA") {
    condition.filter <- parseCondition(Condition.str, simData.df)
    simData.df <- simData.df[condition.filter, ]
  }

  return(simData.df)
}

# Parse a condition side (either a baseline variable or a value)
parseOneConditionSide <- function(side.str, simData.df) {
  # if it's a number
  if (!is.na(suppressWarnings(as.numeric(side.str)))) {
    return(as.numeric(side.str))

    # if it's a baseline variable - return the vector of baseline values
  } else {
    if (side.str %in% names(simData.df)) {
      return(simData.df[, side.str])
    }
  }
}

# Parse one condition defined as "Baseline == X"
# return a Boolean : true if the patient respect the condition
parseOneCondition <- function(cond.str, simData.df) {
  if (cond.str == "TRUE") {
    return(TRUE)
  }
  # Extract the 2 sides of the condition (should be "==")
  if (grepl("==", cond.str)) {
    condSides <- stringr::str_split(cond.str, " == ")[[1]]
    operator <- "=="
  } else {
    if (grepl("<=", cond.str)) {
      condSides <- stringr::str_split(cond.str, " <= ")[[1]]
      operator <- "<="
    } else {
      if (grepl(">=", cond.str)) {
        condSides <- stringr::str_split(cond.str, " >= ")[[1]]
        operator <- ">="
      } else {
        if (grepl(">", cond.str)) {
          condSides <- stringr::str_split(cond.str, " > ")[[1]]
          operator <- ">"
        } else {
          if (grepl("<", cond.str)) {
            condSides <- stringr::str_split(cond.str, " < ")[[1]]
            operator <- "<"
          } else {
            stop(paste("operand not implemented:", cond.str))
          }
        }
      }
    }
  }

  # parse each side of the condition and check the equality
  if (length(condSides == 2)) {
    varLeft.str <- parseOneConditionSide(condSides[1], simData.df)
    varRight.str <- parseOneConditionSide(condSides[2], simData.df)

    if (operator == "==") filter <- varLeft.str == varRight.str
    if (operator == "<=") filter <- varLeft.str < varRight.str | varLeft.str == varRight.str
    if (operator == ">=") filter <- varLeft.str > varRight.str | varLeft.str == varRight.str
    if (operator == "<") filter <- varLeft.str < varRight.str
    if (operator == ">") filter <- varLeft.str > varRight.str

    if (sum(filter) == 0) print(warning(paste("No patient with", cond.str)))
    return(filter)
  } else {
    stop(paste("condition parsing", cond.str))
  }
}

# Parse all the conditions defining the subgroups
parseCondition <- function(conditions.str, simData.df) {
  # Extract conditions separated by And
  if (grepl("And", conditions.str)) {
    cond.split <- stringr::str_split(conditions.str, " And |\\)|\\(")[[1]]
    cond.vec <- Filter(f = function(x) {
      stringr::str_length(x) > 1
    }, cond.split)

    # Parse each condition
    (filterCond.list <- lapply(cond.vec, FUN = function(x) {
      parseOneCondition(x, simData.df)
    }))
    # Apply the `And`
    return(Reduce("&", filterCond.list))
  } else {
    (grepl("Or", conditions.str))
    cond.split <- stringr::str_split(conditions.str, " Or |\\)|\\(")[[1]]
    cond.vec <- Filter(f = function(x) {
      stringr::str_length(x) > 1
    }, cond.split)

    # Parse each condition
    (filterCond.list <- lapply(cond.vec, FUN = function(x) {
      parseOneCondition(x, simData.df)
    }))

    # Apply the `And`
    return(Reduce("|", filterCond.list))
  }
}



# Lipoproteins ----
# Compute the predicted boostrapped inter-quartile range interval for each lipoprotein measures
addPredictedLipValue <- function(refData.df, simData.df) {
  ## Extract info from the Variable name
  refData.df$Arm <- sapply(stringr::str_split(refData.df$Variable, "_"), FUN = function(i) i[1])
  refData.df$Measure <- sapply(stringr::str_split(refData.df$Variable, "_"), FUN = function(i) i[2])
  refData.df$Measure <- gsub(pattern = "[W|D|M|Y][0-9]+", replacement = "", x = refData.df$Measure)


  ## Prepare the prediction dataset
  predData.df <- refData.df[, !names(refData.df) %in% c("Value", "LowRange", "HighRange")]
  predData.df$Ref <- NA

  # Add Origin column
  refData.df$Origin <- "Reference"
  predData.df$Origin <- "Simulation"

  # Compute the predicte interval for each line of observation
  predictedLipdDyn.list <- lapply(c(1:nrow(predData.df)),
    FUN = function(l) {
      line <- predData.df[l, ]
      print(paste(l,"/",nrow(predData.df)))
      if (is.na(line$RangeType)) {
        line$RangeType <- ""
      }
      if (line$ValueType == "Mean") {
        if (line$RangeType == "IQR") {
          stop(print(paste("Mean is not compatible with IQR")))
        }
        value.fun <- "mean"
      } else if (line$ValueType == "Median") {
        value.fun <- "median"
      } else {
        stop(print(paste("Value Type not recognized")))
      }
      if (line$RangeType == "IQR") {
        line <- line[, !names(line) %in% c("RangeType")]
        pred.line <- getBootIQR(
          simData.df = simData.df,
          Measure.name = as.character(line$Variable),
          Condition.str = as.character(line$Condition),
          Arm.str = line$Arm,
          Time.val = line$Time
        )
      } else if (line$RangeType %in% c("CI","95CI")) {
        line <- line[, !names(line) %in% c("RangeType")]
        pred.line <- getBootCI(
          simData.df = simData.df,
          Measure.name = line$Variable,
          Condition.str = line$Condition,
          Arm.str = line$Arm,
          Time.val = line$Time,
          fun = value.fun,
          sampleSize.int = onearm.sampleSize.int
        )
      } else if (line$RangeType == "") {
        line <- line[, !names(line) %in% c("RangeType")] # if no range type, we use 95%CI as default
        pred.line <- getBootCI(
          simData.df = simData.df,
          Measure.name = line$Variable,
          Condition.str = line$Condition,
          Arm.str = line$Arm,
          Time.val = line$Time,
          fun = value.fun,
          sampleSize.int = onearm.sampleSize.int
        )
      } else {
        stop(print(paste("addPredictedLipValue: RangeType not implemented", line$RangeType)))
      }
      res <- cbind(line, pred.line)
      return(res)
    }
  )
  predictedLipdDyn.df <- do.call("rbind", predictedLipdDyn.list)
  return(rbind(refData.df, predictedLipdDyn.df))
}


# Plot observed lipoprotein data with predicted interval
# And compute coverage and precision on each arm
compareLipDyn <- function(allLipDyn.df, Measure.name, Condition.str, ValueType.str) {
  # Extract measure and condition of interest
  toPlot.df <- allLipDyn.df[allLipDyn.df$Measure == Measure.name & allLipDyn.df$ValueType == ValueType.str, ]
  if (!is.na(Condition.str)) toPlot.df <- toPlot.df[toPlot.df$Condition == Condition.str, ]

  # Add pretty arm name for plot
  toPlot.df$Arm.plot <- gsub(paste0(Study.name, "-"), "", toPlot.df$Arm)
  toPlot.df$Arm.plot <- factor(toPlot.df$Arm.plot, levels = gsub(paste0(Study.name, "-"), "", arms))

  # Replace "Reference" origin by the study name
  toPlot.df$Data <- factor(toPlot.df$Origin, levels = c("Reference", "Simulation"))
  toPlot.df$Data <- gsub("Reference", paste("Observed in",Study.name), toPlot.df$Data)
  toPlot.df$Data <- gsub("Simulation", paste("Simulated in",Study.name,"Vpop"), toPlot.df$Data)

  #  Exract the observed and predicted range types
  obsRangeType.str <- unique(toPlot.df[toPlot.df$Origin == "Reference", "RangeType"])
  simRangeType.str <- unique(toPlot.df[toPlot.df$Origin == "Simulation", "RangeType"])
  simValType.str <- unique(toPlot.df[toPlot.df$Origin == "Simulation", "ValueType"])
  
  name.str = variables.metadata.df[variables.metadata.df$Measure == Measure.name,"Variable"]
  unit.str =  variables.metadata.df[variables.metadata.df$Measure == Measure.name,"Unit"]

  if (Condition.str %in% c("","TRUE",NA)) {
    cond.legend <- ""
  } else {
    cond.legend <- paste("Condition:", Condition.str)
  }

  plot <- Plot.NovaTheme(
    ggplot(toPlot.df, aes(x = Time, y = Value, by = Data, col = Arm.plot, linetype = Data)) +
      geom_ribbon(
        data = toPlot.df,
        aes(
          ymin = LowRange, ymax = HighRange,
          col = Arm.plot, fill = Arm.plot, alpha = Data
        ), linetype = "blank"
      ) +
      geom_point(data = toPlot.df, aes(shape = Data), lwd = 1, size = 3) +
      geom_line(data = toPlot.df, lwd = 1) +
      labs(
        title = name.str,
        subtitle = paste(cond.legend, simValType.str, "and", simRangeType.str),
        x = "Time (year)",
        y = paste(name.str,"(",unit.str,")"),
        fill = "Arm",
        col = "Arm"
      ) +
      scale_color_manual(values = c(ctrl.color, trmt.color)) +
      scale_fill_manual(values = c(ctrl.color, trmt.color)) +
      scale_alpha_manual(values = c(0.1, 0.3)) +
      scale_linetype_manual(values = c("dashed", "solid")) +
      scale_shape_manual(values = c(16, 3)) +
      ylim(0,NA)
  )

  print(plot)

  res.list <- list()
  for(arm.name in unique(toPlot.df$Arm)){
    error.df = data.frame(
      error = mean(toPlot.df[toPlot.df$Arm == arm.name & toPlot.df$Origin == "Simulation","Value"] - toPlot.df[toPlot.df$Arm == arm.name & toPlot.df$Origin == "Reference","Value"], na.rm = T)
    )
    names(error.df) = paste("error",arm.name,sep="-")
    res.list <- append(res.list, error.df)

  }


  # Compute coverage and precision
    for (arm.name in unique(toPlot.df$Arm)) {

      # Remove NA value to compute coverage and precision
      timesWithNA = toPlot.df[is.na(toPlot.df$Value),"Time"]
      toPlot.df = toPlot.df[!(toPlot.df$Time %in% timesWithNA),]

      if (nrow(toPlot.df) > 0 & (all(!is.na(toPlot.df[toPlot.df$Arm == arm.name, "LowRange"]))) & (all(!is.na(toPlot.df[toPlot.df$Arm == arm.name, "HighRange"])))) {
        cAndp.df <- Coverage.Precision.Validation(
          SimUpperBound = toPlot.df[toPlot.df$Arm == arm.name & toPlot.df$Origin == "Simulation", "HighRange"],
          SimLowerBound = toPlot.df[toPlot.df$Arm == arm.name & toPlot.df$Origin == "Simulation", "LowRange"],
          TimeVecSim = toPlot.df[toPlot.df$Arm == arm.name & toPlot.df$Origin == "Simulation", "Time"],
          ObsUpperBound = toPlot.df[toPlot.df$Arm == arm.name & toPlot.df$Origin == "Reference", "HighRange"],
          ObsLowerBound = toPlot.df[toPlot.df$Arm == arm.name & toPlot.df$Origin == "Reference", "LowRange"],
          TimeVecObs = toPlot.df[toPlot.df$Arm == arm.name & toPlot.df$Origin == "Simulation", "Time"]
        )
      } else {
        warning(paste(Measure.name, Condition.str, ":no simulation data"))
        cAndp.df <- (data.frame(coverage = NA, precision = NA))
      }
      names(cAndp.df) <- paste(arm.name, names(cAndp.df))
      res.list <- append(res.list, cAndp.df)
    }

  compareRes.df <- do.call("cbind", res.list)


  return(list(plot =plot, data = compareRes.df))
}


## Compare one relative change of lipoproteins
compareOneLipChangeLine <- function(lc.line, simData.df) {
  # apply condition
  if (!is.na(lc.line$Condition)) simData.df <- applyCondition.fct(simData.df, lc.line$Condition)

  # extract simulated values of interest
  ctrl.val <- simData.df[, lc.line$Control.Variable]
  trmt.val <- simData.df[, lc.line$Treatment.Variable]

  # I don't know why I need this conversion but I need it..
  lc.line$Value = as.numeric(lc.line$Value)
  lc.line$HighRange = as.numeric(lc.line$HighRange)
  lc.line$LowRange = as.numeric(lc.line$LowRange)

  # Compute change confidence interval

  if(lc.line$ChangeType == "Relative") simData.df$change <- 100 * (ctrl.val - trmt.val) / ctrl.val
  if(lc.line$ChangeType == "Absolute") simData.df$change <- ctrl.val - trmt.val
    if (is.na(lc.line$RangeType)) {
    lc.line$RangeType <- ""
  }
  if (lc.line$RangeType %in% c("95CI", "")) {
    if (lc.line$ValueType == "Mean") {
      value.fun <- "mean"
    } else if (lc.line$ValueType == "Median") {
      value.fun <- "median"
    } else {
      stop(print(paste("Value Type not recognized")))
    }
    change.df <- getBootCI(simData.df = simData.df, Measure.name = "change", Condition.str = NA, Arm.str = NA, Time.val = NA, fun = value.fun, sampleSize.int = onearm.sampleSize.int)
  }else{
    if (lc.line$RangeType == "IQR"){
      change.df <- getBootIQR(simData.df = simData.df, Measure.name ="change", Condition.str = NA, Arm.str = NA, Time.val = NA)
  }else {
    stop(print(paste("compareOneLipChangeLine: RangeType not implemented", lc.line$RangeType)))
  }}

  if (!is.na(lc.line$LowRange)) {
    # Compute coverage and precision
    coveragePrecison.vec <- Coverage.Precision.Validation(
      SimUpperBound = change.df$HighRange,
      SimLowerBound = change.df$LowRange,
      TimeVecSim = 1,
      ObsUpperBound = lc.line$HighRange,
      ObsLowerBound = lc.line$LowRange,
      TimeVecObs = 1
    )
  } else {
    coveragePrecison.vec <- list(coverage = NA, precision = NA)
  }

  # Gather the results in a data frame
  res.df <- data.frame(
    Control   = lc.line$Control.Variable,
    Treatment = lc.line$Treatment.Variable,
    ChangeType = lc.line$ChangeType,
    Condition = lc.line$Condition,
    Size      = nrow(simData.df),
    ValueType = lc.line$ValueType,
    RangeType = lc.line$RangeType,
    Ref       = lc.line$Ref,
    obs       = paste0(round(lc.line$Value,1)," ",lc.line$Unit," [", round(lc.line$LowRange, 1), "-", round(lc.line$HighRange, 1), "]"),
    sim       = paste0(round(change.df$Value, 1)," ",lc.line$Unit," [", round(change.df$LowRange, 1), "-", round(change.df$HighRange, 1), "]"),
    coverage  = paste(coveragePrecison.vec$coverage, "%"),
    precision = paste(coveragePrecison.vec$precision, "%"),
    simInObs  = between(change.df$Value, lower = lc.line$LowRange, upper = lc.line$HighRange, NAbounds = NA),
    obsInSim  = between(lc.line$Value, lower = change.df$LowRange, upper = change.df$HighRange, NAbounds = NA),
    error     = round(change.df$Value - lc.line$Value, 2)
  )

  return(res.df)
}




## Compare one relative change of lipoproteins
compareOneLipDiffLine <- function(ld.line, simData.df) {
  # apply condition for ref pop and subgroup
  if (!is.na(ld.line$RefCondition)) ref.simData.df <- applyCondition.fct(simData.df, ld.line$RefCondition) else ref.simData.df <- simData.df
  if (!is.na(ld.line$SubGroupCondition)) sg.simData.df <- applyCondition.fct(simData.df, ld.line$SubGroupCondition) else stop("compareOneLipChangeLine: no subgroup condition")
  
  if (nrow(sg.simData.df) > 0 &
      ld.line$Variable %in% names(simData.df) &
      all(!is.na(ref.simData.df[,ld.line$Variable])) &
      all(!is.na(sg.simData.df[,ld.line$Variable]))) {
    if (ld.line$ValueType == "Mean") {
      value.fun <- "mean"
    } else if (ld.line$ValueType == "Median") {
      value.fun <- "median"
    } else {
      stop(print(paste("Value Type not recognized")))
    }
    # extract simulated values of interest
    ref.val <- boot(ref.simData.df[, ld.line$Variable], function(u, i) get(value.fun)(u[i]), R = nBoot)$t
    subgroup.val <- boot(sg.simData.df[, ld.line$Variable], function(u, i) get(value.fun)(u[i]), R = nBoot)$t
    
    # Compute change confidence interval
    diff.vec <- subgroup.val / ref.val
    
    if (ld.line$RangeType %in% c("CI","95CI")) {
      diff.df <- data.frame(
        mean = mean(diff.vec),
        LowRange = quantile(diff.vec, probs = c(0.025), na.rm = T),
        HighRange = quantile(diff.vec, probs = c(0.975), na.rm = T)
      )
    } else {
      stop(print(paste("compareOneLipChangeLine: RangeType not implemented", ld.line$RangeType)))
    }
    
    # Compute coverage and precision
    coveragePrecison.vec <- Coverage.Precision.Validation(
      SimUpperBound = diff.df$HighRange,
      SimLowerBound = diff.df$LowRange,
      TimeVecSim = 1,
      ObsUpperBound = ld.line$HighRange,
      ObsLowerBound = ld.line$LowRange,
      TimeVecObs = 1
    )
    
    # Gather the results in a data frame
    res.df <- data.frame(
      Variable = ld.line$Variable,
      RefCondition = ld.line$RefCondition,
      SubGrouprCondition = ld.line$SubGroupCondition,
      SubGroupSize = nrow(sg.simData.df),
      RangeType = ld.line$RangeType,
      Ref = ld.line$Ref,
      obs = paste0(ld.line$Value, "% [", ld.line$LowRange, "-", ld.line$HighRange, "]"),
      sim = paste0(round(diff.df$mean, 4), "% [", round(diff.df$LowRange, 4), "-", round(diff.df$HighRange, 4), "]"),
      coverage = paste(coveragePrecison.vec$coverage, "%"),
      precision = paste(coveragePrecison.vec$precision, "%"),
      simInObs = between(diff.df$mean, lower = ld.line$LowRange, upper = ld.line$HighRange, NAbounds = NA),
      obsInSim = between(ld.line$Value, lower = diff.df$LowRange, upper = diff.df$HighRange, NAbounds = NA),
      error = round(diff.df$mean - ld.line$Value, 2)
    )
  } else { # No patient in the simulated subgroup
    res.df <- data.frame(
      Variable = ld.line$Variable,
      RefCondition = ld.line$RefCondition,
      SubGrouprCondition = ld.line$SubGroupCondition,
      SubGroupSize = nrow(sg.simData.df),
      RangeType = ld.line$RangeType,
      Ref = ld.line$Ref,
      obs = paste0(ld.line$Value, "% [", round(ld.line$LowRange, 1), "-", round(ld.line$HighRange, 1), "]"),
      sim = NA,
      coverage = NA,
      precision = NA,
      simInObs = NA,
      obsInSim = NA,
      error = NA
    )
  }
  
  return(res.df)
}

# Clinical events ----

## Event rate
compareOneOccurenceLine <- function(occ.line,simData.df){
  # Apply the subgroup condition
  condition.str <- occ.line$Condition
  if (stringr::str_length(condition.str) > 0) {
    condition.filter <- parseCondition(condition.str, simData.df)
    simData.df <- simData.df[condition.filter, ]
  }

  Arm.name = occ.line$Arm
  Event.name = gsub("is","",occ.line$Event)
  TimeVar.name = sort(grep(paste0(Arm.name,"_(t|(obsT))","imeOf", Event.name), names(simData.df),value = T))[1] # priority to obsTime if existing
  occPred.dbl = sum(simData.df[simData.df[,TimeVar.name] <= occ.line$`Time (years)`,paste0(occ.line$Arm,"_",occ.line$Event)]) / nrow(simData.df)

  error.dbl = occPred.dbl - occ.line$Value

  # Gather the results in a data frame
  res.df <- data.frame(
    Outcome = occ.line$Event,
    Arm = occ.line$Arm,
    Time = occ.line$`Time (years)`,
    Condition = occ.line$Condition,
    SampleSize = nrow(simData.df),
    Ref = occ.line$Ref,
    obs = occ.line$Value * 100,
    sim = round(occPred.dbl,3) * 100,
    error = round(error.dbl*100,1),
    ErrorBelowOne = abs(error.dbl) <= 0.01
  )

  return(res.df)
}

twoArmsOccurenceTable = function(oneArm.df, simData.df){
  df = as.data.frame(t(as.data.frame(oneArm.df)))
  compare.df = df[df$Arm == Ctrl.Arm.name, ]
  names(compare.df) <- gsub("obs","obs-control",names(compare.df))
  names(compare.df) <- gsub("sim","sim-control",names(compare.df))
  names(compare.df) <- gsub("error","error-control",names(compare.df))
  names(compare.df) <- gsub("ErrorBelowOne","errorBelowOne-control",names(compare.df))

  head(compare.df)

  compare.df[,"obs-trmt"] = NA
  compare.df[,"sim-trmt"] = NA
  compare.df[,"error-trmt"] = NA
  compare.df[,"errorBelowOne-trmt"] = NA

  for(r in c(1:nrow(compare.df))){
    trmt.line = df[df$Arm == Trmt.Arm.name &
                     df$Outcome == unlist(compare.df[r,"Outcome"]) &
                     df$Time == unlist(compare.df[r,"Time"]) &
                     df$Condition == unlist(compare.df[r,"Condition",])
                   ,]

    if(nrow(trmt.line) == 1){
      compare.df[r,"obs-trmt"] = trmt.line$obs
      compare.df[r,"sim-trmt"] = trmt.line$sim
      compare.df[r,"error-trmt"] = trmt.line$error
      compare.df[r,"errorBelowOne-trmt"] = trmt.line$ErrorBelowOne
    }else{
      # Even if no ref, Compute simulated value
      Event.name = gsub("is","",unlist(compare.df[r,"Outcome"]))
      TimeVar.name = sort(grep(paste0(Trmt.Arm.name,"_(t|(obsT))","imeOf", Event.name), names(simData.df),value = T))[1] # priority to obsTime if existing
      occPred.dbl = sum(simData.df[simData.df[,TimeVar.name] <= unlist(compare.df[r,"Time"]),paste0(Trmt.Arm.name,"_",unlist(compare.df[r,"Outcome"]))]) / nrow(simData.df)

      compare.df[r,"sim-trmt"] = round(occPred.dbl*100, 1)
    }

  }

  return(compare.df)
}

## Hazard ratio ----

compareOneHrLine <- function(hr.line, simData.df, Ctrl.Arm.name, Trmt.Arm.name) {
  # follow-up time
  min.time = hr.line$Time.Min.Bound
  if(is.na(min.time)) min.time = 0
  max.time <- hr.line$Time.Max.Bound
  
  # Event observed
  Event.name <- hr.line$Event

  # Apply the subgroup condition
  condition.str <- hr.line$Condition
  if (stringr::str_length(condition.str) > 0) {
    condition.filter <- parseCondition(condition.str, simData.df)
    simData.df <- simData.df[condition.filter, ]
  }
  # remove patients dead at tmin
  simData.df <- simData.df[simData.df[,paste0(Ctrl.Arm.name,"_timeOfCvDeath")] > min.time & simData.df[,paste0(Trmt.Arm.name,"_timeOfCvDeath")] > min.time,]

  if (nrow(simData.df) > 0) {
    # Define the occurrence and time to variables for the hazard ratio evaluation
    eventAndTime.df <- cbind(
      addVarOutcomeBefore(simData.df = simData.df, Arm.name = Ctrl.Arm.name, Event.name = Event.name, min.time, max.time),
      addVarOutcomeBefore(simData.df = simData.df, Arm.name = Trmt.Arm.name, Event.name = Event.name, min.time, max.time)
    )

    # Compute the simulated hazard ratio
    EventVar.Sim <- paste0("is", Event.name, "Y", min.time, "-", max.time)
    TimeVar.Sim <- paste0("timeOf", Event.name, "Y", min.time, "-", max.time)
  
    simHR.df <- getSimHR(simData.df = eventAndTime.df, Ctrl.Arm.name, Trmt.Arm.name,
      EventVar.Sim, TimeVar.Sim,
      Year.int = max.time,
      bootSampleSize.int = mean(hr.line$SampleSizeControl, hr.line$SampleSizeTreatment)
    )
  } else {
    simHR.df <- data.frame(
      HR = NA,
      HR.min = NA,
      HR.max = NA
    )
  }

  # Compute coverage and precision
  if (!is.na(simHR.df$HR)) {
    coveragePrecison.vec <-
      Coverage.Precision.Validation(
        SimUpperBound = simHR.df[, "HR.max"],
        SimLowerBound = simHR.df[, "HR.min"],
        TimeVecSim = c(1), # fake time data because there is only one point
        ObsUpperBound = hr.line$High95CI,
        ObsLowerBound = hr.line$Low95CI,
        TimeVecObs = c(1) # fake time data because there is only one point
      )
  } else {
    coveragePrecison.vec <- list(coverage = NA, precision = NA)
  }


  # Gather the results in a data frame
  res.df <- data.frame(
    Id = hr.line$Id,
    Outcome = hr.line$Event,
    Condition = hr.line$Condition,
    tMin = hr.line$Time.Min.Bound,
    tMax = hr.line$Time.Max.Bound,
    Size.Sim = nrow(simData.df),
    Size.CtrlObs = hr.line$SampleSizeControl,
    Size.TrmtObs = hr.line$SampleSizeTreatment,
    Ref = hr.line$Ref,
    RangeType = "95CI",
    obs = paste0(round(hr.line$Value, 2), " [", round(hr.line$Low95CI, 2), "-", round(hr.line$High95CI, 2), "]"),
    sim = paste0(round(simHR.df$HR, 2), " [", round(simHR.df$HR.min, 2), "-", round(simHR.df$HR.max, 2), "]"),
    coverage = paste(coveragePrecison.vec$coverage, "%"),
    precision = paste(coveragePrecison.vec$precision, "%"),
    simInObs = if (!is.na(simHR.df$HR)) {
      between(simHR.df$HR, lower = hr.line$Low95CI, upper = hr.line$High95CI, NAbounds = NA)
    } else {
      NA
    },
    obsInSim = if (!is.na(simHR.df$HR)) {
      between(hr.line$Value, lower = simHR.df$HR.min, upper = simHR.df$HR.max, NAbounds = NA)
    } else {
      NA
    },
    error = round(simHR.df$HR - hr.line$Value, 2)
  )
  return(res.df)
}


compareOneSubgroupHrLine <- function(sghr.line, simData.df, Arm.name) {
  # follow-up time
  min.time = sghr.line$Time.Min.Bound
  max.time <- sghr.line$Time.Max.Bound

  # Event observed
  Event.name <- sghr.line$Event

  # apply condition for ref pop and subgroup
  if (!is.na(sghr.line$ConditionNoRisk)) noRF.simData.df <- applyCondition.fct(simData.df, sghr.line$ConditionNoRisk) else noRF.simData.df <- simData.df
  if (!is.na(sghr.line$ConditionRisk)) rf.simData.df <- applyCondition.fct(simData.df, sghr.line$ConditionRisk) else stop("compareOneLipChangeLine: no subgroup condition")
  varName = paste(Arm.name, paste0("is", Event.name), sep = "_")
  timeVarName = sort(grep(paste0("(t|(obsT))","imeOf", Event.name), names(simData.df),value = T))[1] # prioiryt to obsTime if existing

  if(!varName %in% names(simData.df) | !timeVarName %in% names(simData.df)){
    print(paste(varName,"or",timeVarName,"not in simulated dataset"))
  }else{
    if(sghr.line$Low95CI == sghr.line$High95CI){
      print(paste(sghr.line$Id,":",sghr.line$Value,"[",sghr.line$Low95CI,"-",sghr.line$High95CI,"]"))
    }else{
    if (nrow(rf.simData.df) > 0 & sum(noRF.simData.df[, paste(Arm.name, paste0("is", Event.name), sep = "_")]) > 0) {
      # Define the occurence and time to variables for the hazard ratio evaluation
      ref.eventAndTime.df <- addVarOutcomeBefore(simData.df = noRF.simData.df, Arm.name = Arm.name, Event.name = Event.name, min.time, max.time)
      sg.eventAndTime.df <- addVarOutcomeBefore(simData.df = rf.simData.df, Arm.name = Arm.name, Event.name = Event.name, min.time, max.time)

      # Compute the simulated hazard ratio
      EventVar.Sim  <- paste0("is", Event.name, "Y", min.time, "-", max.time)
      TimeVar.Sim   <- paste0("timeOf", Event.name, "Y", min.time, "-", max.time)
      ref.boot      <- getBootHR(ref.eventAndTime.df,
                                 event.name         = paste(Arm.name, EventVar.Sim, sep = "_"),
                                 time.name          = paste(Arm.name, TimeVar.Sim, sep = "_"),
                                 Year.int           = max.time,
                                 bootSampleSize.int = sghr.line$SampleSizeWithout
      )
      sg.boot       <- getBootHR(sg.eventAndTime.df,
                                 event.name         = paste(Arm.name, EventVar.Sim, sep = "_"),
                                 time.name          = paste(Arm.name, TimeVar.Sim, sep = "_"),
                                 Year.int           = max.time,
                                 bootSampleSize.int = sghr.line$SampleSizeWith
      )

      # Computing HRs based on Cox
      resCox.dtf <- BootstrapCox.fun(BootDataCtrl.dtf = ref.boot$SimulatedDataAllSamples.dtf,
                                     BootDataTrtd.dtf = sg.boot$SimulatedDataAllSamples.dtf)

      resCox.dtf <- resCox.dtf %>% mutate_at(c('HR', 'HR95LoCI', 'HR95HiCI'), as.numeric)

      sim.bootHR      <- resCox.dtf$HR
      sim.mean.bootHR <- mean(resCox.dtf$HR, na.rm = T)
      sim.low.bootHR  <- mean(resCox.dtf$HR95LoCI, na.rm = T)
      sim.high.bootHR <- mean(resCox.dtf$HR95HiCI, na.rm = T)
      
      # Compute coverage and precision
        coveragePrecison.vec <-
          Coverage.Precision.Validation(
            SimUpperBound = sim.high.bootHR,
            SimLowerBound = sim.low.bootHR,
            TimeVecSim = c(1), # fake time data because there is only one point
            ObsUpperBound = sghr.line$High95CI,
            ObsLowerBound = sghr.line$Low95CI,
            TimeVecObs = c(1) # fake time data because there is only one point
          )


        # Gather the results in a data frame
        res.df <- data.frame(
          Id = NA, # sghr.line$Id,
          Outcome = sghr.line$Event,
          RiskCondition = sghr.line$ConditionRisk,
          SampleSize = paste("RF:",sghr.line$SampleSizeWith,"no RF:",sghr.line$SampleSizeWithout),
          Ref = sghr.line$Ref,
          RangeType = "CI",
          obs = paste0(round(sghr.line$Value, 2), " [", round(sghr.line$Low95CI, 2), "-", round(sghr.line$High95CI, 2), "]"),
          sim = paste0(round(sim.mean.bootHR, 2), " [", round(sim.low.bootHR, 2), "-", round(sim.high.bootHR, 2), "]"),
          coverage = paste(coveragePrecison.vec$coverage, "%"),
          precision = paste(coveragePrecison.vec$precision, "%"),
          simInObs = between(sim.mean.bootHR, lower = sghr.line$Low95CI, upper = sghr.line$High95CI, NAbounds = NA),
          obsInSim = between(sghr.line$Value, lower = sim.low.bootHR, upper = sim.high.bootHR, NAbounds = NA),
          error = round(sim.mean.bootHR - sghr.line$Value, 2)
        )

        return(res.df)
      }}}

    res.df <- data.frame(
      Id = NA, #sghr.line$Id,
      Outcome = sghr.line$Event,
      RiskCondition = sghr.line$ConditionRisk,
      SampleSize = paste("RF:",sghr.line$SampleSizeWith,"no RF:",sghr.line$SampleSizeWithout),
      Ref = sghr.line$Ref,
      RangeType = "CI",
      obs = paste0(round(sghr.line$Value, 2), " [", round(sghr.line$Low95CI, 2), "-", round(sghr.line$High95CI, 2), "]"),
      sim = NA,
      coverage = NA,
      precision = NA,
      simInObs = NA,
      obsInSim = NA,
      error = NA
    )
}


# Return a data frame with 2 variable:
# - event occurrence before a given time
# - time of event occurrence before a given time
addVarOutcomeBefore <- function(simData.df, Arm.name, Event.name, min.time, max.time) {
  # Define the new variables name
  newEvent.name <- paste0(Arm.name, "_is", Event.name, "Y", min.time, "-", max.time)
  newTime.name <- paste0(Arm.name, "_timeOf", Event.name, "Y", min.time, "-", max.time)

  # Define the event occurrence variable
  var.name <- paste0(Arm.name, "_is", Event.name)
  if (var.name %in% names(simData.df)) {
    newEvent.vec <- simData.df[, paste0(Arm.name, "_is", Event.name)]
  } else {
    stop(paste("addVarOutcomeBefore:", var.name, "not in simulated dataset"))
  }

  timeVarName = sort(grep(paste0(Arm.name,"_(t|(obsT))","imeOf", Event.name), names(simData.df),value = T))[1] #priority to obsTime if existing

  # if the time of event is not between min.time and max.time. event is set at 0.
  newEvent.vec[!between(simData.df[, timeVarName], min.time, max.time)] <- 0

  # Define the time to event variable
  newTime.vec <- simData.df[, timeVarName]
  newTime.vec[!between(simData.df[, timeVarName], min.time, max.time)] <- max.time

  res <- data.frame(event = newEvent.vec, time = newTime.vec)
  names(res) <- c(newEvent.name, newTime.name)

  return(res)
}


getBootHR <- function(data.df,
                      event.name, # outcome variable
                      time.name, # time variable
                      Year.int, # year of interest
                      bootSampleSize.int) {
  
  survSimData.df <- data.frame(
    status = data.df[, event.name],
    time = data.df[, time.name]
  )

  ## Computation of prediction interval for the lifetime distribution function

  SurvArm.list <- Survival.Bootstrap(
    SampleSize = bootSampleSize.int,
    NbBootstraps = nBoot,
    SurvivalData = survSimData.df,
    NbClusters = 4,
    Seed      = seed.int
  )

  # Extraction of survival at time of interest

  SurvAtEOT.df <- data.frame()

  for (boot in unique(SurvArm.list$FullResults.dtf$BootNb)) {
    boot.dtf <- subset(SurvArm.list$FullResults.dtf, BootNb == boot)
    SurvAtEOTrow <- boot.dtf[which.min(abs(boot.dtf$time - Year.int)), ]
    SurvAtEOT.df <- rbind(SurvAtEOT.df, SurvAtEOTrow)
  }
  hrate.vec <- -log(SurvAtEOT.df$surv) / SurvAtEOT.df$time

  # UPDATE: Also returning "SimulatedDataAllSamples.dtf" from Survival.Bootstrap()
  # so that Cox based HR can be computed

  Output.lst <- list(hrate.vec                   = hrate.vec,
                     SimulatedDataAllSamples.dtf = SurvArm.list$SimulatedDataAllSamples.dtf)

  return(Output.lst)
}

# Compute a bootstrapped distribution of hazard rate
getSimHrateBootDistrib <- function(simData.df, # Simulated data
                                   Arm.name, # arm name
                                   EventVar.Sim, # Outcome of interest
                                   TimeVar.Sim, # Time-to-event of interest
                                   Year.int, # Year of interest
                                   bootSampleSize.int) {
  ## Survival data set

  if (any(simData.df[, paste(Arm.name, TimeVar.Sim, sep = "_")] <= 0)) {
    stop("At least one time-to-event value is negative in the reference arm (first). This function can't deal with such values")
  }

  bootRes.lst <- getBootHR(simData.df,
                           event.name = paste(Arm.name, EventVar.Sim, sep = "_"),
                           time.name = paste(Arm.name, TimeVar.Sim, sep = "_"), # time variable
                           Year.int, # year of interest
                           bootSampleSize.int
  )
  return(bootRes.lst)
}

# Compute a hazard ratio from simulated data from bootstrapp
getSimHR <- function(simData.df, # Simulated data
                     Ctrl.Arm.name, # Reference arm
                     Trmt.Arm.name, # Comparator arm
                     EventVar.Sim, # Outcome of interest
                     TimeVar.Sim, # Time-to-event of interest
                     Year.int, # Year of interest
                     bootSampleSize.int) {
  ## Control arm hazard rate bootstrapped distribution
  BootCtrl.res   <- getSimHrateBootDistrib(simData.df, Ctrl.Arm.name, EventVar.Sim, TimeVar.Sim, Year.int, bootSampleSize.int)

    ## Treatment arm hazard rate bootstrapped distribution
  BootTrtd.res   <- getSimHrateBootDistrib(simData.df, Trmt.Arm.name, EventVar.Sim, TimeVar.Sim, Year.int, bootSampleSize.int)

  ## Bootstraps for control arm
  ctrlArm.dtf <- BootCtrl.res$SimulatedDataAllSamples.dtf

  ## Bootstraps for treatment arm
  trtdArm.dtf <- BootTrtd.res$SimulatedDataAllSamples.dtf

  ## Compute empirical distribution of Hazard ratios
  if (sum(BootCtrl.res$hrate.vec > 0)) {
    # Computing HRs based on Cox
    resCox.dtf <- BootstrapCox.fun(BootDataCtrl.dtf = ctrlArm.dtf,
                                   BootDataTrtd.dtf = trtdArm.dtf)

    resCox.dtf <- resCox.dtf %>% mutate_at(c('HR', 'HR95LoCI', 'HR95HiCI'), as.numeric)

    sim.bootHR      <- resCox.dtf$HR
    sim.mean.bootHR <- mean(resCox.dtf$HR, na.rm = T)
    sim.low.bootHR  <- mean(resCox.dtf$HR95LoCI, na.rm = T)
    sim.high.bootHR <- mean(resCox.dtf$HR95HiCI, na.rm = T)

    ## Hazard ratio prediction intervals computation

    # HRci.df <- data.frame(
    #   HR     = mean(HRtrial.vec),
    #   HR.min = quantile(HRtrial.vec, 0.025),
    #   HR.max = quantile(HRtrial.vec, 0.975)
    # )
    HRci.df <- data.frame(
      HR     = mean(resCox.dtf$HR, na.rm = T),
      HR.min = mean(resCox.dtf$HR95LoCI, na.rm = T),
      HR.max = mean(resCox.dtf$HR95HiCI, na.rm = T)
    )
  } else { # No event in control arm
    HRci.df <- data.frame(
      HR     = NA,
      HR.min = NA,
      HR.max = NA
    )
  }
  return(HRci.df)
}


## Curves coverage and precision (+ plot) ----


isOutputInSimData <- function(arm.name, eventName, names.vec){
  if(substr(eventName,0,2) == "is") eventName = substring(eventName,first = 3)
  return(
    (paste(arm.name, paste0("is",eventName), sep = "_") %in% names.vec)
    &
      ((paste(arm.name, paste0("timeOf",eventName), sep = "_") %in% names.vec)
       |
         (paste(arm.name, paste0("obsTimeOf",eventName), sep = "_") %in% names.vec))
  )
}

# Execute the comparison of one life time distribution curve with a observation
compareOneLdcLine <- function(ldc.line, simData.df, Study.name, ldc.inputFolder.str, startTime = 0, subtitle.str = NA) {
  # Get reference Data
  refData.df <- getRefLdc.fun(Study.name, ldc.line$File.name, ldc.inputFolder.str)
  simData.df <- applyCondition.fct(simData.df, ldc.line$Condition)

  # initialize df
  res.df <- data.frame(
    Outcome   = ldc.line$Event.name,
    Condition = ldc.line$Condition,
    Size      = nrow(simData.df),
    Nevent_ref= sum(refData.df$eventId),
    Ref       = ldc.line$Ref
  )
  # Plot LDC and coverage and precision
  # Get a dataframe with lifetime distribution of both simulated and ref data for each arm
  if(nrow(simData.df) > 0){
    survAllArms.list <- lapply(arms, FUN = function(a) {
      surv.df <- getOneArmSurvFitDataBoot95CI(simData.df, refData.df,
                                              Arm.name = a,
                                              EventVar.Sim = ldc.line$Event.name,
                                              TimeVar.Sim = ldc.line$Event.time
      )
    })
    surv.df <- do.call("rbind", survAllArms.list)



    plot <- plotTwoArmsLifetimeDistribution(
      survData.df = surv.df,
      Study.name = Study.name,
      Event.name = ldc.line$Event.name,
      Subtitle = if_else(is.na(subtitle.str),paste(ldc.line$Condition),subtitle.str),
      startTime = startTime
    )
    print(plot)
  }


  for (arm.name in arms) {
    if (sum(simData.df[, paste(arm.name, ldc.line$Event.name, sep = "_")]) > 0) {
      # print(paste("Compare LDC:",ldc.line$File.name,arm.name)) # log
      coveragePrecison.vec <- getOneArmCoverageAndPrecision(arm.name, df = surv.df)


      # Compute Combo long rank test
      if(sum(refData.df$eventId) > 1){
        test.boot <-
          bootTestingOneArm(
            simData.df = simData.df,
            refData.df = refData.df,
            Arm.name = arm.name,
            EventVar.Sim = ldc.line$Event.name,
            TimeVar.Sim = ldc.line$Event.time,
            bootSampleSize = mean(table(refData.df$Arm)),
            nBoot = nBoot,
            testType = "LRCombo"
          )
      }else{
        test.boot = list(ratio = NA)
      }
      # Gather the results in a data frame
      oneArm.res.df <- data.frame(
        coverage  = paste(coveragePrecison.vec$coverage, "%"),
        precision = paste(coveragePrecison.vec$precision, "%"),
        cLogRankTest =  paste(round(test.boot$ratio * 100), "%")
      )
    } else { # ne event in the subgroup
      oneArm.res.df <- data.frame(
        coverage  = NA,
        precision = NA,
        cLogRankTest =  NA
      )
    }


    ## add the arm name in the raw name.
    names(oneArm.res.df) <- paste0(names(oneArm.res.df), gsub(Study.name, "", arm.name))

    ## concatenate data frames
    res.df <- cbind(res.df, oneArm.res.df)
  }

  return(list(data = res.df, plot = plot))
}





# Execute the comparison of one life time distribution curve with a observation
compareOneTimePointFromLdcLine <- function(ldc.line, simData.df, Study.name, ldc.inputFolder.str, timePoint, startTime = 0, subtitle.str = NA) {
  # Get reference Data
  refData.df <- getRefLdc.fun(Study.name, ldc.line$File.name, ldc.inputFolder.str)
  simData.df <- applyCondition.fct(simData.df, ldc.line$Condition)
  
  # initialize df
  res.df <- data.frame(
    Outcome   = ldc.line$Event.name,
    Condition = ldc.line$Condition,
    Size      = nrow(simData.df),
    Nevent_ref= sum(refData.df$eventId),
    Ref       = ldc.line$Ref
  )
  # Plot LDC and coverage and precision
  # Get a dataframe with lifetime distribution of both simulated and ref data for each arm
  if(nrow(simData.df) > 0){
    survAllArms.list <- lapply(arms, FUN = function(a) {
      surv.df <- getOneArmSurvFitDataBoot95CI(simData.df, refData.df,
                                              Arm.name = a,
                                              EventVar.Sim = ldc.line$Event.name,
                                              TimeVar.Sim = ldc.line$Event.time
      )
    })
    surv.df <- do.call("rbind", survAllArms.list)
  }
  
  res.list <- list()
  for(o in unique(surv.df$Origin)){
    for(a in unique(surv.df$Arm)){
      data.df = surv.df[surv.df$Origin == o & surv.df$Arm == a,]
      oneTimePoint.df = data.frame(
        Cumul = approx(x = data.df$Time,y = data.df$Cumul,xout = timePoint)$y,
        Cumul_highStd = approx(x = data.df$Time,y = data.df$Cumul_highStd,xout = timePoint)$y,
        Cumul_LowStd = approx(x = data.df$Time,y = data.df$Cumul_lowStd,xout = timePoint)$y,
        Origin = o,
        Arm = a)
      
      res.list = append(res.list, list(oneTimePoint.df))
  }}
  
  res.df = do.call("rbind",res.list)
  res.df$Time = timePoint
  res.df$Outcome   = ldc.line$Event.name
  res.df$Condition = ldc.line$Condition
  
  return(res.df)
}



# Read the reference raw lifetime distribution data
getRefLdc.fun <- function(Study.name, file.name, ldc.inputFolder.str) {
  filePath.str = paste0(ldc.inputFolder.str, file.name, ".csv")
  print(paste("reading: ",filePath.str))
  return(read.csv(filePath.str))
}



# Compute Lifetime distribution function for simulated and ref data and merge them in a dataset
# Uses the 95% boot strapped CI provided by Survival.Bootstrap
getOneArmSurvFitDataBoot95CI <- function(simData.df, # Simulated Dataset
                                         refData.df, # Reference Dataset
                                         Arm.name, # Scenario name from `Arm` column
                                         EventVar.Sim, # Event measure name from `Arm_EventVar` in simulated data
                                         TimeVar.Sim # Time of observation measure name from `Arm_TimeVar` in simulated data
) {
  if(! paste0(Arm.name, "_", TimeVar.Sim) %in% names(simData.df)) print(paste0(Arm.name, "_", TimeVar.Sim," not in simData"))
  if(! paste0(Arm.name, "_", EventVar.Sim) %in% names(simData.df)) print(paste0(Arm.name, "_", EventVar.Sim," not in simData"))
  simLDC.boot <- Survival.Bootstrap(
    SampleSize = mean(table(refData.df$Arm)), # mean size between arms
    SurvivalData = data.frame(
      time = simData.df[, paste0(Arm.name, "_", TimeVar.Sim)],
      status = simData.df[, paste0(Arm.name, "_", EventVar.Sim)]
    ),
    NbBootstraps = nBoot,
    Replace = T,
    NbClusters = 4,
    Seed       = seed.int
  )
  sim.surv <- simLDC.boot$PredictionIntevalLT.dtf

  # Fit Lifetime distribution function from reference data
  ref.surv <- survfit(Surv(
    time = refData.df[refData.df$Arm == Arm.name, ]$time,
    event = refData.df[refData.df$Arm == Arm.name, ]$eventId
  ) ~ 1)

  # Create a data frame with Lifetime distribution function data
  res.df <- rbind(
    data.frame(
      Time          = sim.surv$time,
      Cumul         = sim.surv$surv,
      Cumul_highStd = sim.surv$upper,
      Cumul_lowStd  = sim.surv$lower,
      Arm           = Arm.name,
      Origin        = "Simulation"
    ),
    data.frame(
      Time = sim.surv$time,
      Cumul = stats::approx(
        x = ref.surv$time,
        y = 1 - ref.surv$surv,
        xout = sim.surv$time
      )$y,
      Cumul_highStd = stats::approx(
        x = ref.surv$time,
        y = 1 - ref.surv$lower,
        xout = sim.surv$time
      )$y,
      Cumul_lowStd = stats::approx(
        x = ref.surv$time,
        y = 1 - ref.surv$upper,
        xout = sim.surv$time
      )$y,
      Arm = Arm.name,
      Origin = "Reference"
    )
  )

  return(res.df)
}


# Plot Lifetime distribution function Data by arms and origin
plotTwoArmsLifetimeDistribution <- function(
    survData.df, # Simulation dataset with only the arm of interest data
    Study.name, # Pretty study name for labels
    Event.name, # Pretty event name for labels
    Subtitle = NA,
    startTime = 0 # time to start the follow-up if we want to "cut" the beginning of the curve
    ){
  survData.df$Arm <- gsub(paste0(Study.name, "-"), "", survData.df$Arm)
  survData.df$Arm <- factor(survData.df$Arm, levels = gsub(paste0(Study.name, "-"), "", arms))

  survData.df$Origin <- gsub("Reference", Study.name, survData.df$Origin)
  survData.df$Origin <- factor(survData.df$Origin, levels = c(Study.name, "Simulation"))

  survData.df$ColorType <- paste(survData.df$Arm, survData.df$Origin, sep = "-")
  survData.df$ColorType[grep(Study.name, survData.df$ColorType)] <- Study.name

  if(startTime > 0){
    Subtitle = paste(Subtitle, "start at", startTime, "year")

    for(a in unique(survData.df$Arm)){
      for(o in unique(survData.df$Origin)){

        startTimeAO = max(survData.df[survData.df$Time <= startTime & survData.df$Arm == a & survData.df$Origin == o,]$Time)
        cumulAtStartTimeAO = survData.df[survData.df$Time == startTimeAO & survData.df$Arm == a & survData.df$Origin == o,"Cumul"]
        # print(cumulAtStartTimeAO)
        survData.df[survData.df$Arm == a & survData.df$Origin == o,]$Cumul = survData.df[survData.df$Arm == a & survData.df$Origin == o,]$Cumul - cumulAtStartTimeAO
        survData.df[survData.df$Arm == a & survData.df$Origin == o,]$Cumul_highStd = survData.df[survData.df$Arm == a & survData.df$Origin == o,]$Cumul_highStd - cumulAtStartTimeAO
        survData.df[survData.df$Arm == a & survData.df$Origin == o,]$Cumul_lowStd = survData.df[survData.df$Arm == a & survData.df$Origin == o,]$Cumul_lowStd - cumulAtStartTimeAO
        survData.df[survData.df$Arm == a & survData.df$Origin == o,]$Time = survData.df[survData.df$Arm == a & survData.df$Origin == o,]$Time - startTimeAO
      }
    }
  }
  survData.df = survData.df[survData.df$Time >= 0,]

  summary(survData.df)

  pretty.name = variables.metadata.df[variables.metadata.df$Measure == Event.name,"Variable"]
  pretty.unit = variables.metadata.df[variables.metadata.df$Measure == Event.name,"Unit"]
  
  # Convert in percentage
  survData.df$Cumul = 100 * survData.df$Cumul
  survData.df$Cumul_highStd = 100 * survData.df$Cumul_highStd
  survData.df$Cumul_lowStd = 100 * survData.df$Cumul_lowStd
  
  Plot.NovaTheme(
    ggplot(survData.df, aes(x = Time, y = Cumul)) +
      geom_line(data = survData.df, aes(x = Time, y = Cumul, col = Arm, by = Origin, linetype = Origin)) +
      geom_ribbon(
        data = survData.df[survData.df$Origin == Study.name, ],
        aes(
          ymin = Cumul_lowStd, ymax = Cumul_highStd,
          col = Arm, fill = ColorType
        ), alpha = 0.5, linetype = "blank"
      ) +
      geom_ribbon(
        data = survData.df[survData.df$Origin == "Simulation", ],
        aes(
          ymin = Cumul_lowStd, ymax = Cumul_highStd,
          col = Arm, fill = ColorType
        ), alpha = 0.5, linetype = "blank"
      ) +
      labs(
        title = pretty.name, # title = paste("Lifetime distribution of", pretty.name),
        subtitle = Subtitle,
        x = "Time (year)",
        y = paste(pretty.name, "(%)"),
        fill = NULL,
        linetype = NULL
      ) +
      scale_color_manual(values = c(ctrl.color, trmt.color)) +
      scale_fill_manual(values = c(ref.color, trmt.color, ctrl.color)) +
      scale_linetype_manual(values = c("dashed", "solid")) +
      facet_wrap(~Arm) +
      guides(col = FALSE)
  )
}


# Compute coverage and precision for a given arm
getOneArmCoverageAndPrecision <- function(a, df) {
  Coverage.Precision.Validation(
    SimUpperBound = df[df$Arm == a & df$Origin == "Simulation", "Cumul_highStd"],
    SimLowerBound = df[df$Arm == a & df$Origin == "Simulation", "Cumul_lowStd"],
    TimeVecSim = df[df$Arm == a & df$Origin == "Simulation", "Time"],
    ObsUpperBound = df[df$Arm == a & df$Origin == "Reference", "Cumul_highStd"],
    ObsLowerBound = df[df$Arm == a & df$Origin == "Reference", "Cumul_lowStd"],
    TimeVecObs = df[df$Arm == a & df$Origin == "Reference", "Time"]
  )
}



# Compute validation bootstrap functions on one arm
bootTestingOneArm <- function(simData.df, refData.df, Arm.name, EventVar.Sim, TimeVar.Sim, bootSampleSize, nBoot, testType) {
  # Simulated data
  survSimData.df <- data.frame(
    status = simData.df[, paste(Arm.name, EventVar.Sim, sep = "_")],
    time   = simData.df[, paste(Arm.name, TimeVar.Sim, sep = "_")],
    study  = "Simulation"
  )

  # ref data
  survRefData.df <- data.frame(
    status = refData.df[refData.df$Arm == Arm.name, "eventId"],
    time   = refData.df[refData.df$Arm == Arm.name, "time"],
    study  = "Reference"
  )

  test <- tryCatch(
    {
      Validation.BootTesting(
        group1             = survSimData.df,
        group2             = survRefData.df,
        sampleSize         = bootSampleSize,
        replace.grp1       = T,
        replace.grp2       = T,
        alpha              = 0.05,
        reps               = nBoot,
        test               = testType,
        non.significant    = T,
        VarOfInterest.surv = "time",
        Status.surv        = "status",
        Strata.surv        = "study",
      )
    },
    error=function(cond) {
      return(list(ratio = NA))
    }
  )

  return(test)
}

# Outcomes and time-to-event data

outcomes.vec   <- c("isMACE", "isCvDeath", "isMI", "isIscStroke", "isMajorAdverseLimbEvent")
outcomes_n.vec <- c("3P-MACE", "CV death", "MI", "Ischemic Stroke", "MALE")
times.vec      <- c("timeOfMACE", "timeOfCvDeath", "obsTimeOfMI", "obsTimeOfIscStroke", "obsTimeOfMajorAdverseLimbEvent")


# Compute the HR of a risk factors defined as a condition.
# Plot cumulative incidence per subgroup
rf.hr = function(simData.df, Condition.str, arm.name, outcome="isMACE", study="FOURIER"){
  event.time <- times.vec[which(outcome == outcomes.vec)]

  TempDataForCox.dtf = data.frame(
    time = simData.df[,paste0(arm.name, "_", event.time)],
    status = simData.df[,paste0(arm.name, "_", outcome)],
    rf = parseCondition(Condition.str, simData.df)
  )
  if (nrow(TempDataForCox.dtf[!TempDataForCox.dtf$rf,]) < 1) {return(data.frame(HR="NA", Range="NA"))}
  if (nrow(TempDataForCox.dtf[TempDataForCox.dtf$rf,]) < 1) {return(data.frame(HR="NA", Range="NA"))}

  cox.mdl <- coxph(Surv(time, status) ~ rf, data = TempDataForCox.dtf)

  rangeHR <-round(exp(confint(cox.mdl)), digits = 3)
  rfhr.df = data.frame(
    HR = round(exp(coef(cox.mdl)[1]), digits = 3),
    Range = paste("[",rangeHR[1],"-",rangeHR[2],"]")
    )

  surv_fit_1 <- survfit(Surv(time, status) ~ 1, data = TempDataForCox.dtf[TempDataForCox.dtf$rf,])
  surv_fit_0 <- survfit(Surv(time, status) ~ 1, data = TempDataForCox.dtf[!TempDataForCox.dtf$rf,])


  surv_fit.df = rbind(
    data.frame(
      time = surv_fit_1$time,
      cumul = 1 - surv_fit_1$surv,
      low = 1 - surv_fit_1$upper,
      up = 1 - surv_fit_1$lower,
      rf = 1
    ),
    data.frame(
      time = surv_fit_0$time,
      cumul = 1 - surv_fit_0$surv,
      low = 1 - surv_fit_0$upper,
      up = 1 - surv_fit_0$lower,
      rf = 0
    )
  )
  surv_fit.df$rf = factor(surv_fit.df$rf,levels = c("0","1"))
  # Create a ggplot object
  plot = Plot.NovaTheme(
    ggplot(surv_fit.df, aes(x = time, y = cumul, color = rf, group = rf, fill = rf)) +
      geom_line() +
      geom_ribbon(aes(ymin = low, ymax = up), linetype = "blank", alpha =0.1) +
      labs(title = paste(outcome, "cumulative incidence"), subtitle = paste("in placebo arm in", study),
           x = "Time", y = "Cumulative Incidence", color = Condition.str, fill = Condition.str) +
      theme(legend.position = "bottom")  # Set legend position to the bottom
  )

  # print(plot)

  return(rfhr.df)
}

# Compute Vpop table ----
getVpopTable <- function(subsample.margins.df, simData.df, condition.line = NULL) {
  dedicatedTarget.bool = T # Flag to check if a dedicated target exist. True by default for entire pop

  # Apply the subgroup condition
  if (!is.null(condition.line)) {
    condition.str = condition.line$Condition
    condition.name = condition.line$Category

    condition.filter <- parseCondition(condition.str, simData.df)
    simData.df <- simData.df[condition.filter, ]

    # check if there is dedicated targets
    if(length(grep(paste0("(",condition.name,")"), names(subsample.margins.df))) > 1){
      subsample.margins.df$Variable = subsample.margins.df[,paste0("Variable (",condition.name,")")]
      subsample.margins.df$Target = subsample.margins.df[,paste0("Target (",condition.name,")")]
      subsample.margins.df$TargetValue = subsample.margins.df[,paste0("TargetValue (",condition.name,")")]
      subsample.margins.df$SummaryType = subsample.margins.df[,paste0("SummaryType (",condition.name,")")]
      subsample.margins.df$Ref = subsample.margins.df[,paste0("Ref (",condition.name,")")]

      dedicatedTarget.bool = T
    }else{
      subsample.margins.df$Target = ""
      subsample.margins.df$TargetValue = ""
      subsample.margins.df$Ref = ""

      dedicatedTarget.bool = F
    }


  }else{
    condition.str = ""
  }
  if (nrow(simData.df) == 0) {
    return(print(paste("No patient in subgroup:", condition.str)))
  }

  vpop.table.df <- subsample.margins.df
  # Compute the simulated value
  vpop.table.df[vpop.table.df$SummaryType == "Proportion", "SimulatedValue"] <-
    sapply(
      vpop.table.df[vpop.table.df$SummaryType == "Proportion", "Variable"],
      FUN = function(var.name) {
        round(mean(simData.df[, var.name]) * 100, 2)
      }
    )

  if ("FOURIER-OLE-switch_isMACE" %in% vpop.table.df$Variable) {
    vpop.table.df[vpop.table.df$Variable %in% c("FOURIER-OLE-switch_isMACE", "FOURIER-OLE-evolocumab_isMACE"), "SimulatedValue"] <-
      sapply(
        vpop.table.df[vpop.table.df$Variable %in% c("FOURIER-OLE-switch_isMACE", "FOURIER-OLE-evolocumab_isMACE"), "Variable"],
        FUN = function(var.name) {
          # print(var.name)
          round(mean(simData.df[, stringr::str_replace(var.name, "is", "timeOf")] <= 2.2) * 100, 2)
        }
      )
  }

  vpop.table.df[vpop.table.df$SummaryType == "MeanSD", "SimulatedValue"] <-
    sapply(
      vpop.table.df[vpop.table.df$SummaryType == "MeanSD", "Variable"],
      FUN = function(var.name) {
        # print(var.name)
        round(mean(simData.df[, var.name]), 2)
      }
    )

  vpop.table.df[vpop.table.df$SummaryType == "MedianIQR", "SimulatedValue"] <-
    sapply(
      vpop.table.df[vpop.table.df$SummaryType == "MedianIQR", "Variable"],
      FUN = function(var.name) {
        # print(var.name)
        round(median(simData.df[, var.name]), 2)
      }
    )
  # Compute pretty output
  vpop.colName <- paste("Virtual Population (N=", nrow(simData.df), ")")
  vpop.table.df[, vpop.colName] <- ""
  vpop.table.df[vpop.table.df$SummaryType == "MedianIQR", vpop.colName] <-
    sapply(
      vpop.table.df[vpop.table.df$SummaryType == "MedianIQR", "Variable"],
      FUN = function(var.name) {
        median <- round(median(simData.df[, var.name]), 2)
        quantiles <- round(quantile(simData.df[, var.name], probs = c(0.25, 0.75)), 1)
        return(paste0(median, " (", quantiles[1], "-", quantiles[2], ")"))
      }
    )

  vpop.table.df[vpop.table.df$SummaryType == "MeanSD", vpop.colName] <-
    sapply(
      vpop.table.df[vpop.table.df$SummaryType == "MeanSD", "Variable"],
      FUN = function(var.name) {
        mean <- round(mean(simData.df[, var.name]), 2)
        sd <- round(sd(simData.df[, var.name]), 2)
        return(paste0(mean, " (", sd, ")"))
      }
    )

  vpop.table.df[vpop.table.df$SummaryType == "Proportion", vpop.colName] <-
    vpop.table.df[vpop.table.df$SummaryType == "Proportion", "SimulatedValue"]



  # Compute the error
  if(dedicatedTarget.bool){
    vpop.table.df$Error <- vpop.table.df$SimulatedValue - vpop.table.df$TargetValue
    vpop.table.df$Acceptable <- abs(vpop.table.df$Error) <= vpop.table.df$AbsoluteTol

    vpop.table.df$Acceptable = as.character(vpop.table.df$Acceptable)
    vpop.table.df$Acceptable = replace_na(vpop.table.df$Acceptable,replace="NOCHECK")
    vpop.table.df <- replace(vpop.table.df, is.na(vpop.table.df), "")

  }

  # Remove Target if no check so that it is not shown
  # vpop.table.df[vpop.table.df$Acceptable == "NOCHECK","Target"] = ""

  # Target is defined only on the entire population
  if (condition.str == "" | dedicatedTarget.bool) {
    col.names <- c("Description", "Target", "Ref", vpop.colName, "Acceptable")
  } else {
    col.names <- c("Description", vpop.colName)
  }

  return(vpop.table.df[, col.names])
}
