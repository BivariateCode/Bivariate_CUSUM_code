rm(list = ls(all.names = TRUE))
library(parallel)
library(data.table)

#### Helper functions #########
make_prior_df <- function(priorsMat, rhoVals, label) {
  data.frame(
    prior_set = label,
    prior_row = seq_len(nrow(priorsMat)),
    prior1    = priorsMat[,1],
    prior2    = priorsMat[,2],
    prior3    = priorsMat[,3],
    prior4    = priorsMat[,4],
    rho1      = rhoVals[1],
    rho2      = rhoVals[2]
  )
}
make_method_Biv8_stable <- function(p1, priorParams, rhoVec, condCdfs_stable, tab) {
  list(
    name = "Biv8_Stable_CUSUM",
    rl = function(maxIter, h) {
      Biv8_Stable_CUSUM(maxIter = maxIter, h = h, p1 = p1, priorParams = priorParams, rhoVec = rhoVec, condCdfs = condCdfs_stable, startingValues = tab[sample.int(nrow(tab), 1), ])
    }
  )
}
make_ARL_policy <- function(targetARL0, fewerIter = FALSE) {
  stopifnot(is.numeric(targetARL0), length(targetARL0) == 1)
  
  if (targetARL0 >= 1450 && targetARL0 <= 1550) {
    if (fewerIter) {
      Ns <- c(350, 400)
      eps <- c(300, 250)
    } else {
      Ns <- c(500, 4000, 45000)
      eps <- c(67, 24, 7)
    }
    
  } else if (targetARL0 >= 995 && targetARL0 <= 1050) {
    if (fewerIter) {
      Ns <- c(300, 350, 600)
      eps <- c(300, 250, 100)
    } else {
      Ns <- c(500, 4000, 45000)
      eps <- c(62, 21, 5)
    }
    
  } else if (targetARL0 >= 450 && targetARL0 <= 550) {
    if (fewerIter) {
      Ns <- c(200, 250)
      eps <- c(300, 250)
    } else {
      Ns <- c(500, 4000, 38500)
      eps <- c(53, 15, 3)
    }
    
  } else if (targetARL0 >= 24 && targetARL0 <= 65) {
    if (fewerIter) {
      Ns <- c(100, 200)
      eps <- c(10, 9)
    } else {
      Ns <- c(225, 2500, 10000)
      eps <- c(5, 1, 0.25)
    }
    
  } else if (targetARL0 >= 99 && targetARL0 <= 250) {
    if (fewerIter) {
      Ns <- c(100, 200)
      eps <- c(18, 12)
    } else {
      Ns <- c(500, 4000, 38500)
      eps <- c(50, 11, 1.5)
    }
    
  } else {
    stop("No ARL policy defined for targetARL0 = ", targetARL0)
  }
  
  list(
    Ns = Ns,
    epsilonVec = eps
  )
}
make_clims_job <- function(method, targetARL0, h_init, fewerIter = FALSE, maxIterMult = 50, maxOuterIter = 300, keep_history = TRUE, keep_RLs = TRUE, logger = NULL) {
  stopifnot(is.list(method), is.function(method$rl))
  
  policy <- make_ARL_policy(targetARL0, fewerIter)
  
  list(
    # solver inputs
    h_init       = h_init,
    targetARL0   = targetARL0,
    RL_funct     = method$rl,
    
    # policy
    Ns           = policy$Ns,
    epsilonVec   = policy$epsilonVec,
    
    # controls
    maxIterMult  = maxIterMult,
    maxOuterIter = maxOuterIter,
    
    # diagnostics
    keep_history = keep_history,
    keep_RLs     = keep_RLs,
    
    # side effects
    logger       = logger
  )
}
make_clims_job_Biv8 <- function(row, condCdfs_stable, tab, h_init, fewerIter = TRUE, logfile = NULL, logger = NULL) {
  
  method <- make_method_Biv8_stable(
    p1          = row$p1,
    priorParams = c(row$prior1, row$prior2, row$prior3, row$prior4),
    rhoVec      = c(row$rho1, row$rho2),
    condCdfs_stable = condCdfs_stable,
    tab         = tab
  )
  
  job = make_clims_job(
    method       = method,
    targetARL0   = row$ARL0,
    h_init       = h_init,
    fewerIter    = fewerIter,
    logger = logger
  )
  
  job$logfile <- logfile
  return(job)
}
stable_fname <- function(row) {
  sprintf(
    "Biv8_stable_pr_%.3g_%.3g_%.3g_%.3g_rho_%.3g_%.3g_p1_%.2f.rds",
    row$prior1, row$prior2, row$prior3, row$prior4,
    row$rho1, row$rho2,
    row$p1
  )
}
IC_Rl_fname <- function(row)  {
  sprintf(
    "Biv8_stable_pr_%.3g_%.3g_%.3g_%.3g_rho_%.3g_%.3g_p1_%.2f_arl_%s.txt",
    row$prior1, row$prior2, row$prior3, row$prior4,
    row$rho1, row$rho2,
    row$p1, row$ARL0
  )
}
make_file_logger <- function(filepath) {
  first_write <- TRUE
  force(filepath)
  
  function(...) {
    cat(
      ...,
      "\n",
      file   = filepath,
      append = !first_write
    )
    first_write <<- FALSE
    invisible(NULL)
  }
}
########### Main functions ###########
Biv8_Stable_CUSUM <- function(maxIter, h, p1, priorParams, rhoVec, condCdfs, startingValues)   {
  # p1 is P(V = 1 | V != 2). We assume L = 1 case (this first obs) occurs every other obs deterministically.
  # Doing the above is a very reasonable approximation for all of the IC cases considered in the Zwetsloot paper
  
  Cvec = startingValues[1:8]
  Svec1 = startingValues[9:16]
  Svec2 = startingValues[17:24]
  Svec3 = startingValues[25:32]
  Tvec1 = startingValues[33:40]
  Tvec2 = startingValues[41:48]
  Tvec3 = startingValues[49:56]
  
  alpha0_inc = priorParams[1]
  beta0_inc = priorParams[2]
  alpha0_dec = priorParams[3]
  beta0_dec = priorParams[4]
  
  incVal = rhoVec[1]
  decVal = rhoVec[2]
  
  l1_defaultVals = c(rep(incVal,4), rep(decVal, 4))
  l2_defaultVals = c(rep(incVal,2), rep(decVal,2), rep(incVal,2), rep(decVal,2))
  temp = c(incVal, decVal, incVal, decVal)
  l3_defaultVals = c(temp, temp)
  
  lambda_vec1 = l1_defaultVals
  
  lambda_vec2 = l2_defaultVals
  
  lambda_vec3 = l3_defaultVals
  
  lamPosInd1 = 1:4
  lamNegInd1 = 5:8
  lamPosInd2 = c(1,2,5,6)
  lamNegInd2 = setdiff(1:8, lamPosInd2)
  lamPosInd3 = seq(1,7,2)
  lamNegInd3 = setdiff(1:8, lamPosInd3)
  
  b = 0 # Doesn't matter whether the starting value is 0 or 1 because either way the first update will return Svec = = Tvec = 0vec for first iteration since Cvec = 0vec
  
  GetEstimates <- function(Cvec, Svec, Tvec, lamPosInd, lamNegInd, l_defaultVals, Z, initVals = F)  {
    posInd = (Cvec > 0)
    lambda_vec = rep(NA, 8)
    
    if (initVals == F)  { # When itit values is false, we start from a point where we have already updated Svec, Tvec (at this iter we did not update Cvec though) and hence all we need to do is use Svec, Tvec to compute lambdaVec. We DO NOT update anything in the case where initVals == T
      Svec[posInd] = Svec[posInd] + Z
      Svec[!posInd] = 0
      
      Tvec[posInd] = Tvec[posInd] + 1
      Tvec[!posInd] = 0
    }
    
    # Update Estimates 
    # lambda_vec[lamPosInd] <- pmax(incVal, (alpha0_inc + Tvec[lamPosInd])/(Svec[lamPosInd] + beta0_inc) )
    # lambda_vec[lamNegInd] <- pmin(decVal, (alpha0_dec + Tvec[lamNegInd])/(Svec[lamNegInd] + beta0_dec) )
    
    inc_ratio <- (alpha0_inc + Tvec[lamPosInd]) / (Svec[lamPosInd] + beta0_inc)
    inc_ratio <- ifelse(is.finite(inc_ratio), inc_ratio, incVal)
    lambda_vec[lamPosInd] <- pmax(incVal, inc_ratio)
    
    dec_ratio <- (alpha0_dec + Tvec[lamNegInd]) / (Svec[lamNegInd] + beta0_dec)
    dec_ratio <- ifelse(is.finite(dec_ratio), dec_ratio, decVal)
    lambda_vec[lamNegInd] <- pmin(decVal, dec_ratio)
    
    lambda_vec[!posInd] =  l_defaultVals[!posInd] #Only needed if rho != alpha0/beta0 - aka the prior expectation. IOW this is because the prior expectation is NOT equal to the minimum size shift we care about
    
    if (initVals) {
      return(lambda_vec)
    } else  {
      return(list(Svec, Tvec, lambda_vec))
    }
  }
  ResetEstimates <- function(Cvec, Svec, Tvec, lamPosInd, lamNegInd, l_defaultVals)  {
    posInd = (Cvec > 0)
    lambda_vec = rep(NA, 8)
    # cat("\nInside ResetEstimates: \n Before: \n")
    # cat("Svec: ", Svec, "\n")
    # cat("Tvec: ", Tvec, "\n")
    
    Svec[!posInd] = 0
    Tvec[!posInd] = 0
    
    # Update Estimates 
    ## The two lines below don't work when the alpha and beta priors are both 0
    # lambda_vec[lamPosInd] <- pmax(incVal, (alpha0_inc + Tvec[lamPosInd])/(Svec[lamPosInd] + beta0_inc) )
    # lambda_vec[lamNegInd] <- pmin(decVal, (alpha0_dec + Tvec[lamNegInd])/(Svec[lamNegInd] + beta0_dec) )
    
    inc_ratio <- (alpha0_inc + Tvec[lamPosInd]) / (Svec[lamPosInd] + beta0_inc)
    inc_ratio <- ifelse(is.finite(inc_ratio), inc_ratio, incVal)
    lambda_vec[lamPosInd] <- pmax(incVal, inc_ratio)
    
    dec_ratio <- (alpha0_dec + Tvec[lamNegInd]) / (Svec[lamNegInd] + beta0_dec)
    dec_ratio <- ifelse(is.finite(dec_ratio), dec_ratio, decVal)
    lambda_vec[lamNegInd] <- pmin(decVal, dec_ratio)
    
    lambda_vec[!posInd] =  l_defaultVals[!posInd] #Only needed if rho != alpha0/beta0 - aka the prior expectation. IOW this is because the prior expectation is NOT equal to the minimum size shift we care about
    
    # cat("After: \n")
    # cat("Svec: ", Svec, "\n")
    # cat("Tvec: ", Tvec, "\n\n")
    
    return(list(Svec, Tvec, lambda_vec))
  }
  
  
  lambda_vec1 = GetEstimates(Cvec, Svec1, Tvec1, lamPosInd1, lamNegInd1, l1_defaultVals, Z = NA, initVals = T)
  lambda_vec2 = GetEstimates(Cvec, Svec2, Tvec2, lamPosInd2, lamNegInd2, l2_defaultVals, Z = NA, initVals = T)
  lambda_vec3 = GetEstimates(Cvec, Svec3, Tvec3, lamPosInd3, lamNegInd3, l3_defaultVals, Z = NA, initVals = T)
  
  for (i in 1:maxIter)  {
    ### Start of operations done with previous values ###
    if (i%%2 == 0)  { # These are offset because we update using previous value to preserve martingale property
      # L = GetEstimates(Cvec, Svec1, Tvec1, lambda_vec1, lamPosInd1, lamNegInd1, l1_defaultVals, Z1)
      L_res1 = GetEstimates(Cvec, Svec1, Tvec1, lamPosInd1, lamNegInd1, l1_defaultVals, Z1)
      L_res2 = ResetEstimates(Cvec, Svec2, Tvec2, lamPosInd2, lamNegInd2, l2_defaultVals)
      L_res3 = ResetEstimates(Cvec, Svec3, Tvec3, lamPosInd3, lamNegInd3, l3_defaultVals)
      Svec1 = L_res1[[1]]; Tvec1 = L_res1[[2]]; lambda_vec1 = L_res1[[3]]
      Svec2 = L_res2[[1]]; Tvec2 = L_res2[[2]]; lambda_vec2 = L_res2[[3]]
      Svec3 = L_res3[[1]]; Tvec3 = L_res3[[2]]; lambda_vec3 = L_res3[[3]]
    } else if (i%%2 == 1 & i != 1) { # Since this was already done when obtaining stable dist. Therefore we generate and update Svec1
      if (b == 0) {
        L_res1 = ResetEstimates(Cvec, Svec1, Tvec1, lamPosInd1, lamNegInd1, l1_defaultVals)
        L_res2 = GetEstimates(Cvec, Svec2, Tvec2, lamPosInd2, lamNegInd2, l2_defaultVals, Z2)
        L_res3 = ResetEstimates(Cvec, Svec3, Tvec3, lamPosInd3, lamNegInd3, l3_defaultVals)
        Svec1 = L_res1[[1]]; Tvec1 = L_res1[[2]]; lambda_vec1 = L_res1[[3]]
        Svec2 = L_res2[[1]]; Tvec2 = L_res2[[2]]; lambda_vec2 = L_res2[[3]]
        Svec3 = L_res3[[1]]; Tvec3 = L_res3[[2]]; lambda_vec3 = L_res3[[3]]
      } else  {
        L_res1 = ResetEstimates(Cvec, Svec1, Tvec1, lamPosInd1, lamNegInd1, l1_defaultVals)
        L_res2 = ResetEstimates(Cvec, Svec2, Tvec2, lamPosInd2, lamNegInd2, l2_defaultVals)
        L_res3 = GetEstimates(Cvec, Svec3, Tvec3, lamPosInd3, lamNegInd3, l3_defaultVals, Z3)
        Svec1 = L_res1[[1]]; Tvec1 = L_res1[[2]]; lambda_vec1 = L_res1[[3]]
        Svec2 = L_res2[[1]]; Tvec2 = L_res2[[2]]; lambda_vec2 = L_res2[[3]]
        Svec3 = L_res3[[1]]; Tvec3 = L_res3[[2]]; lambda_vec3 = L_res3[[3]]
      }
    }
    ### End of operations done with previous values ###
    
    if (i%%2 == 1)  {
      Z1 = rexp(1, rate = 1)
      Cvec = pmax(0, Cvec + log(lambda_vec1) + (1 - lambda_vec1)*Z1 )
    } else if (i%%2 == 0) {
      b = rbinom(1, 1, p1)
      if (b == 0) {
        Z2 = rexp(1, rate = 1)
        Cvec = pmax(0, Cvec + log(lambda_vec2) + (1 - lambda_vec2)*Z2 )
      } else  {
        Z3 = rexp(1, rate = 1)
        Cvec = pmax(0, Cvec + log(lambda_vec3) + (1 - lambda_vec3)*Z3 )
      }
    }
    
    # Compute p-values
    pvals = 1 - vapply(seq_along(condCdfs), function(j) condCdfs[[j]](Cvec[j]), numeric(1))
    
    # Compute ChiSq values and max selection
    ChiSqVals = -log(pvals)
    ChiSqVal = max(ChiSqVals)
    
    if (ChiSqVal > h) { # I don't think we can use the 1/ARL0 ChiSq quantile to find clim b/c the chi-sq stats are not temporally independent
      break
    }
  }
  return(i)
}
LogInterp <- function(hvals, arls, targetARL)  {
  # Using the form h = k*ln(a) + b
  k = (hvals[2] - hvals[1])/(log(arls[2]) - log(arls[1]))
  b = hvals[1] - k*log(arls[1])
  return(k*log(targetARL) + b)
}
ObtainClims <- function(job) {
  # ---------- 0) Validate + pull fields ----------
  stopifnot(is.list(job))
  
  # Required
  h          <- as.numeric(job$h_init)
  targetARL0 <- as.numeric(job$targetARL0)
  RL_funct   <- job$RL_funct
  Ns         <- job$Ns
  epsilonVec <- job$epsilonVec
  
  if (!is.function(RL_funct)) stop("job$RL_funct must be a function.")
  if (!is.numeric(Ns) || length(Ns) < 1) stop("job$Ns must be a numeric vector of length >= 1.")
  if (!is.numeric(epsilonVec) || length(epsilonVec) < 1) stop("job$epsilonVec must be a numeric vector of length >= 1.")
  if (length(epsilonVec) != length(Ns)) stop("job$Ns and job$epsilonVec must have the same length.")
  if (!is.finite(h) || h <= 0) stop("job$h_init must be positive and finite.")
  if (!is.finite(targetARL0) || targetARL0 <= 0) stop("job$targetARL0 must be positive and finite.")
  
  # Optional knobs (with defaults)
  maxIterMult   <- if (!is.null(job$maxIterMult)) as.numeric(job$maxIterMult) else 50
  maxOuterIter  <- if (!is.null(job$maxOuterIter)) as.integer(job$maxOuterIter) else 300
  maxInnerIter  <- if (!is.null(job$maxInnerIter)) as.integer(job$maxInnerIter) else NA_integer_
  keep_RLs      <- isTRUE(job$keep_RLs)   # store RLs for final iteration
  keep_history  <- isTRUE(job$keep_history)
  
  # Optional logger (injected side-effect)
  logger <- job$logger
  if (!is.function(logger)) {
    logger <- function(...) invisible(NULL)  # no-op
  }
  
  # Derived
  maxIter <- targetARL0 * maxIterMult
  
  # ---------- 1) State variables (same as your code) ----------
  rlBlowUp <- FALSE
  
  ARL_h  <- c(0, 0)   # last two ARL estimates
  hvals  <- c(0, 0)   # last two h values
  nvals  <- c(0, 0)   # sample sizes used for those ARLs
  sizeIndex <- 1
  
  Estim_ARL0 <- -1
  sd_RLs     <- NA_real_
  j_last     <- NA_integer_
  RLs_last   <- NULL
  
  # Optional tracking of the entire outer loop progress
  if (keep_history) {
    hist <- data.frame(
      iter = integer(0),
      N = integer(0),
      eps = numeric(0),
      h = numeric(0),
      arl = numeric(0),
      j = integer(0),
      blowup = logical(0),
      sizeIndex = integer(0),
      stringsAsFactors = FALSE
    )
  }
  
  logger("ObtainClims start | targetARL0=", targetARL0, " | h_init=", h)
  
  # ---------- 2) Outer loop: propose h, estimate ARL(h), update h ----------
  for (i in seq_len(maxOuterIter)) {
    rlBlowUp <- FALSE
    ARL_h_recurs <- 0
    
    N_now <- Ns[sizeIndex]
    eps_now <- epsilonVec[sizeIndex]
    
    # Preallocate as you did (matrix not needed; numeric is fine)
    RLs <- rep.int(-1, N_now)
    
    # ---------- 2a) Inner loop: simulate RLs ----------
    # (Your reason for a for-loop still applies: we may break on blow-up)
    j_end <- 0L
    for (j in seq_len(N_now)) {
      j_end <- j
      
      # KEY ARCH CHANGE: wrapper must make this signature valid
      RLs[j] <- RL_funct(maxIter = maxIter, h = h)
      
      if (j == 1) {
        ARL_h_recurs <- RLs[j]
      } else {
        ARL_h_recurs <- ((j - 1) * ARL_h_recurs + RLs[j]) / j
      }
      
      # Blow-up: same idea, but we record and break
      if (RLs[j] == maxIter) {
        rlBlowUp <- TRUE
        logger("Run length blow-up at iter=", i, " j=", j, " h=", h, " maxIter=", maxIter)
        break
      }
      
      # Optional safety cap independent of N
      if (!is.na(maxInnerIter) && j >= maxInnerIter) {
        logger("Inner cap hit at iter=", i, " j=", j, " h=", h)
        break
      }
    }
    
    # Trim RLs to realized length (j_end)
    RLs <- RLs[seq_len(j_end)]
    j_last <- j_end
    
    logger("iter=", i,
           " | sizeIndex=", sizeIndex,
           " | N=", N_now,
           " | eps=", eps_now,
           " | h=", h,
           " | ARL=", ARL_h_recurs,
           " | j=", j_end,
           " | blowup=", rlBlowUp)
    
    if (keep_history) {
      hist <- rbind(hist, data.frame(
        iter = i, N = N_now, eps = eps_now, h = h, arl = ARL_h_recurs,
        j = j_end, blowup = rlBlowUp, sizeIndex = sizeIndex,
        stringsAsFactors = FALSE
      ))
    }
    
    # ---------- 2b) Update the rolling (h, ARL) pairs ----------
    tol <- 1e-12  # numerical safety
    
    if (ARL_h[1] == 0) {
      
      # first ever entry
      hvals[1]  <- h
      ARL_h[1]  <- ARL_h_recurs
      nvals[1]  <- N_now
      
    } else if (ARL_h[2] == 0) {
      
      if (abs(h - hvals[1]) > tol) {
        # second distinct h
        hvals[2] <- h
        ARL_h[2] <- ARL_h_recurs
        nvals[2] <- N_now
      } else {
        # same h as before → refine MC estimate only
        ARL_h[1] <- ARL_h_recurs
        nvals[1] <- N_now
      }
      
    } else {
      
      if (abs(h - hvals[2]) > tol) {
        # genuine new h → shift window
        hvals[1] <- hvals[2]
        ARL_h[1] <- ARL_h[2]
        nvals[1] <- nvals[2]
        
        hvals[2] <- h
        ARL_h[2] <- ARL_h_recurs
        nvals[2] <- N_now
      } else {
        # same h as most recent → refine only
        ARL_h[2] <- ARL_h_recurs
        nvals[2] <- N_now
      }
      
    }
    
    # if (ARL_h[1] == 0) {
    #   ARL_h[1] <- ARL_h_recurs
    #   hvals[1] <- h
    #   nvals[1] <- N_now
    # } else if (ARL_h[2] == 0) {
    #   ARL_h[2] <- ARL_h_recurs
    #   hvals[2] <- h
    #   nvals[2] <- N_now
    # } else {
    #   hvals[1] <- hvals[2]
    #   ARL_h[1] <- ARL_h[2]
    #   nvals[1] <- nvals[2]
    #   
    #   hvals[2] <- h
    #   ARL_h[2] <- ARL_h_recurs
    #   nvals[2] <- N_now
    # }
    
    # ---------- 2c) Check convergence at current epsilon ----------
    # You had an additional constraint "rlBlowUp == F" — keep it.
    if (!rlBlowUp && abs(ARL_h_recurs - targetARL0) < eps_now) {
      if (sizeIndex == length(epsilonVec)) {
        # Final stage: accept
        Estim_ARL0 <- ARL_h_recurs
        sd_RLs <- sd(RLs)
        if (keep_RLs) RLs_last <- RLs
        logger("Finished | h=", h, " | ARL=", Estim_ARL0, " | sd=", sd_RLs, " | j=", j_end)
        break
      } else {
        # Increase sample size for a more accurate check
        sizeIndex <- sizeIndex + 1L
        # logger("Promoting to next sample size | new sizeIndex=", sizeIndex)
        # Keep h the same for the next stage (matches your intent)
        next
      }
    }
    
    # ---------- 2d) Choose next h (your three cases) ----------
    if (hvals[2] == 0) {
      # Only one point so far: your log rule-of-thumb
      h <- (hvals[1] / log(ARL_h[1])) * log(targetARL0)
      
    } else if ((hvals[1] < hvals[2] && ARL_h[1] > ARL_h[2]) ||
               (hvals[1] > hvals[2] && ARL_h[1] < ARL_h[2])) {
      # Negative slope issue / unreliable estimates
      logger("Case1 (negative slope / noisy estimates)")
      
      if (nvals[1] < nvals[2]) {
        h <- (hvals[2] / log(ARL_h[2])) * log(targetARL0)
      } else {
        h <- rnorm(1, mean = h, sd = (h / 10))
      }
      
    } else {
      # Normal path: log interpolation
      logger("Case8 (LogInterp)")
      h1 <- LogInterp(hvals, ARL_h, targetARL0)
      if (!is.finite(h1) || h1 <= 0) {
        h1 <- rnorm(1, mean = h, sd = (h / 10))
      }
      h <- h1
    }
  }
  
  # ---------- 3) If not converged ----------
  # Keep your convention: h = -1 indicates failure
  if (Estim_ARL0 == -1) {
    Estim_ARL0 <- mean(RLs)
    sd_RLs <- sd(RLs)
    h <- -1
    if (keep_RLs) RLs_last <- RLs
    logger("Failed to converge within maxOuterIter. Returning h=-1.")
  }
  
  # ---------- 4) Return structured result (no side effects) ----------
  out <- list(
    h          = signif(h, 4),
    p1         = job$p1,
    prior1 = job$prior1,
    prior2 = job$prior2,
    prior3 = job$prior3,
    prior4 = job$prior4,
    rho1 = job$rho1,
    rho2 = job$rho2,
    targetARL0 = targetARL0,
    Estim_ARL0 = signif(Estim_ARL0, 4),
    sd_ARL0    = signif(sd_RLs, 4),
    nRuns      = as.integer(j_last),
    sizeIndex  = sizeIndex,
    methodName = job$methodName
  )
  
  if (keep_RLs) out$RLs_last <- RLs_last
  if (keep_history) out$history <- hist
  
  return(out)
}


####### Main #######
setwd("path/to/folder/AdaptiveCUSUM_SimulationCode/") # Make this the path to your folder where this folder exists

#### Setup All Jobs ####

priorsMat1 <- rbind(
  c(22.05, 21, 9.5, 10)  # This is the prior set that was used in the paper
)
rhoVals1 = c(1.05, .95)

MOBE_tarARL0 = c(53.25, 53.25, 24.55, 24.55)
MOBW_tarARL0 = c(62, 60, 26, 26)
Gumbel_tarARL0 = c(53.25, 62, 24.55, 26)


unique_p1s = c(.1, .2, .5)

## These here are the true p1s (aka P(V = 1| V != 2)
# MOBE_p1s = c(.5, .5, .25, .2)
# MOBW_p1s = c(.5, .5, .1, .01)
# Gumbel_p1s = c(.5, .5, .25, .1)

## I will use these p1s instead though
MOBE_p1s = c(.5, .5, .2, .2)
MOBW_p1s = c(.5, .5, .1, .1)
Gumbel_p1s = c(.5, .5, .2, .1)

pairs <- rbind(
  cbind(p1 = MOBE_p1s,   ARL0 = MOBE_tarARL0),
  cbind(p1 = MOBW_p1s,   ARL0 = MOBW_tarARL0),
  cbind(p1 = Gumbel_p1s, ARL0 = Gumbel_tarARL0)
)
unique_pairs <- unique(pairs)

prior_df <- rbind(
  make_prior_df(priorsMat1, rhoVals1, "P1")
)

clim_grid <- merge(prior_df, unique_pairs)

##### Load the clims file #####
clims_dir <- file.path("dataIC")
dir.create(clims_dir, showWarnings = FALSE, recursive = TRUE)

clims_file <- file.path(clims_dir, "Clims_CUSUM.txt")

if (file.exists(clims_file)) {
  clims_tab <- fread(clims_file)
} else {
  clims_tab <- data.table(
    h          = numeric(),
    p1  = numeric(),
    prior1 = numeric(),
    prior2 = numeric(),
    prior3 = numeric(),
    prior4 = numeric(),
    rho1 = numeric(),
    rho2 = numeric(),
    Tar_ARL0   = numeric(),
    Est_ARL0   = numeric(),
    sd_RLs    = numeric(),
    tot_RLs   = integer(),
    functName = character()
  )
}


IC_rls_dir <- file.path("dataIC", "IC_Rls")
if (!dir.exists(IC_rls_dir)) dir.create(IC_rls_dir, recursive = TRUE)
AllIC_DataFiles <- list.files(IC_rls_dir)

######### Run Cluster ########
jobs <- split(clim_grid, seq_len(nrow(clim_grid))) # Need to first turn clim_grid into a list

cl <- makeCluster(length(jobs))
clusterExport(
  cl,
  c(
    "Biv8_Stable_CUSUM",
    "make_method_Biv8_stable",
    "make_file_logger",
    "LogInterp",
    "make_ARL_policy",
    "make_clims_job",
    "stable_fname",
    "IC_Rl_fname",
    "make_clims_job_Biv8",
    "ObtainClims", 
    "IC_rls_dir"
  ),
  envir = environment()
)

results = parLapply(cl, jobs, function(row) {

  # load stable distribution for THIS row
  stable_obj <- readRDS(file.path("dataIC/CUSUM_dist/", stable_fname(row)))
  tab = stable_obj$stable_dist
  cdfs_stable <- lapply(1:8, function(i) {
    Ftemp <- ecdf(tab[, i])
    function(x) (Ftemp(x) - Ftemp(0)) / (1 - Ftemp(0))
  })
  logfile = file.path(
    IC_rls_dir,
    IC_Rl_fname(row)
  )
  job <- make_clims_job_Biv8(
    row              = row,
    condCdfs_stable  = cdfs_stable,
    tab              = tab,   # or wherever it lives
    h_init           = 2,
    fewerIter = F,
    logfile = logfile,
    logger       = make_file_logger(logfile)
  )
  job$p1 = row$p1
  job$methodName = "Biv8_Stable_CUSUM"
  job$prior1 = row$prior1
  job$prior2 = row$prior2
  job$prior3 = row$prior3
  job$prior4 = row$prior4
  job$rho1 = row$rho1
  job$rho2 = row$rho2
  ObtainClims(job)
})

####### Record Output #########


new_rows <- rbindlist(lapply(seq_along(results), function(k) {
  x <- results[[k]]
  data.table(
    h          = x$h,
    p1 = x$p1,
    prior1 = x$prior1,
    prior2 = x$prior2,
    prior3 = x$prior3,
    prior4 = x$prior4,
    rho1 = x$rho1,
    rho2 = x$rho2,
    Tar_ARL0   = x$targetARL0,
    Est_ARL0   = x$Estim_ARL0,
    sd_RLs    = x$sd_ARL0,
    tot_RLs   = x$nRuns,
    functName = x$methodName
  )
}), fill = TRUE)

if (nrow(new_rows) > 0) {
  clims_tab <- rbind(clims_tab, new_rows, fill = TRUE)
  fwrite(clims_tab, clims_file, sep = "\t")
}
