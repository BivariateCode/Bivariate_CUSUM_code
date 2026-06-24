rm(list = ls(all.names = TRUE))
library(parallel)
library(data.table)

#### Helper functions #########
make_method_EWMA <- function(lambda, ICparams) {
  list(
    name = "EWMA_XieEtAl",
    rl = function(maxIter, h) {
      DoEWMA_XieEtAl(maxIter = maxIter, h = h, lambda = lambda, ICparams = ICparams)
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
get_ICparams_from_lmat <- function(lmat, scenario, ic_cols = .ic_cols_Gumb) {
  stopifnot(is.matrix(lmat))
  stopifnot(scenario %in% 1:4)
  col <- ic_cols[scenario]
  # Return exactly what you described: lmat[1:3, ic_col]
  # (vector of length 3: theta1, theta2, delta)
  lmat[1:3, col]
}
make_EWMA_grid_tab <- function(lmat_Gumb, Gumbel_tarARL0, lambdaVals, ic_cols = .ic_cols_Gumb) {
  stopifnot(length(Gumbel_tarARL0) == 4)
  stopifnot(is.numeric(lambdaVals), length(lambdaVals) >= 1)
  
  grid <- expand.grid(
    scenario = 1:4,
    lambda   = lambdaVals,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # add IC column index + IC params + target ARL0
  grid$ic_col <- ic_cols[grid$scenario]
  ICparams_mat <- t(vapply(
    grid$scenario,
    function(s) get_ICparams_from_lmat(lmat_Gumb, s, ic_cols = ic_cols),
    numeric(3)
  ))
  colnames(ICparams_mat) <- c("theta1", "theta2", "delta")
  grid <- cbind(grid, ICparams_mat)
  
  grid$targetARL0 <- Gumbel_tarARL0[grid$scenario]
  grid
}
# Analogue of make_clims_job_Biv8, but for EWMA
make_clims_job_EWMA <- function(row, lmat_Gumb, Gumbel_tarARL0, h_init, fewerIter = TRUE, maxIterMult = 50, maxOuterIter = 300, keep_history = TRUE, keep_RLs = TRUE, logfile = NULL, logger = NULL, ic_cols = .ic_cols_Gumb) {
  
  stopifnot(is.list(row) || is.data.frame(row))
  scenario <- as.integer(row[["scenario"]])
  lambda   <- as.numeric(row[["lambda"]])
  
  ICparams <- get_ICparams_from_lmat(lmat_Gumb, scenario, ic_cols = ic_cols)
  targetARL0 <- Gumbel_tarARL0[scenario]
  
  method <- make_method_EWMA(lambda = lambda, ICparams = ICparams)
  
  if (is.null(logger) && !is.null(logfile)) {
    logger <- make_file_logger(logfile)
  }
  
  job <- make_clims_job(
    method       = method,
    targetARL0   = targetARL0,
    h_init       = h_init,
    fewerIter    = fewerIter,
    maxIterMult  = maxIterMult,
    maxOuterIter = maxOuterIter,
    keep_history = keep_history,
    keep_RLs     = keep_RLs,
    logger       = logger
  )
  
  # Optional: attach metadata (helps debugging / saving results)
  job$logfile   <- logfile
  job$scenario  <- scenario
  job$lambda    <- lambda
  job$ICparams  <- ICparams
  job$methodName = "EWMA_XieEtAl"
  
  job
}
# Build the full list of 20 jobs
make_EWMA_jobs <- function(lmat_Gumb, Gumbel_tarARL0, lambdaVals, h_init, fewerIter = TRUE, logfile_dir = NULL, logger = NULL, ic_cols = .ic_cols_Gumb, ...) {
  
  tab <- make_EWMA_grid_tab(
    lmat_Gumb       = lmat_Gumb,
    Gumbel_tarARL0  = Gumbel_tarARL0,
    lambdaVals      = lambdaVals,
    ic_cols         = ic_cols
  )
  
  # allow scalar h_init or length(tab) vector
  if (length(h_init) == 1) h_init <- rep(h_init, nrow(tab))
  stopifnot(length(h_init) == nrow(tab))
  
  jobs <- vector("list", nrow(tab))
  for (i in seq_len(nrow(tab))) {
    lf <- NULL
    if (!is.null(logfile_dir)) {
      # simple, deterministic logfile naming
      lf <- file.path(
        logfile_dir,
        sprintf("clims_EWMA_scen%d_lambda%s.log",
                tab$scenario[i],
                format(tab$lambda[i], scientific = FALSE, trim = TRUE))
      )
    }
    
    jobs[[i]] <- make_clims_job_EWMA(
      row          = tab[i, ],
      lmat_Gumb    = lmat_Gumb,
      Gumbel_tarARL0 = Gumbel_tarARL0,
      h_init       = h_init[i],
      fewerIter    = fewerIter,
      logfile      = lf,
      logger       = logger,
      ic_cols      = ic_cols,
      ...
    )
  }
  
  list(tab = tab, jobs = jobs)
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
######## Sampling Functions ########
ParseXvec <- function(xvals) {
  if (round(xvals[1], 7) == round(xvals[2], 7)) {
    v = 2
    x = y = xvals[1]
  } else  {
    minInd = which.min(xvals)
    x = xvals[minInd]
    y = xvals[-minInd]
    if (minInd == 1)  {
      v = 0
    } else  {
      v = 1
    }
  }
  return(c(x,y,v))
}
lambertWp_loc <- function(x) {
  if (!is.numeric(x))
    stop("Argument 'x' must be a numeric (real) vector.")
  
  if (length(x) == 1) {
    if (x  < -1/exp(1)) return(NaN)
    if (x == -1/exp(1)) return(-1)
    
    # compute first iteration of $W_0$
    if (x <= 1) {
      eta <- 2 + 2*exp(1)*x;
      f2 <- 3*sqrt(2) + 6 - (((2237+1457*sqrt(2))*exp(1) - 4108*sqrt(2) - 5764)*sqrt(eta)) /
        ((215+199*sqrt(2))*exp(1) - 430*sqrt(2)-796)
      f1 <- (1-1/sqrt(2))*(f2+sqrt(2));
      w0 <- -1 + sqrt(eta)/(1 + f1*sqrt(eta)/(f2 + sqrt(eta)));
    } else {
      w0 = log( 6*x/(5*log( 12/5*(x/log(1+12*x/5)) )) )
    }
    
    # w0 <- 1
    w1 <- w0 - (w0*exp(w0)-x)/((w0+1)*exp(w0)-(w0+2)*(w0*exp(w0)-x)/(2*w0+2))
    
    # iter = 1
    # while(abs(w1-w0) > 1e-11) { #Originally was 1e-15 which caused problems. 1e-13 should be sufficiently large but using 1e-11 to be extra safe
    #   w0 <- w1
    #   cat("w0: ", w0, "\n")
    #   w1 <- w0 - (w0*exp(w0)-x)/((w0+1)*exp(w0)-(w0+2)*(w0*exp(w0)-x)/(2*w0+2))
    #   cat("w1: ", w1, "\n")
    #   cat("absDiff: ", abs(w1-w0), "\n")
    #   iter = iter + 1
    # }
    # cat("-------Total Iterations required: ", iter, "-------\n")
    
    for (i in 1:10) {
      w0 <- w1
      w1 <- w0 - (w0*exp(w0)-x)/((w0+1)*exp(w0)-(w0+2)*(w0*exp(w0)-x)/(2*w0+2))
      if (abs(w1-w0) < 1e-11) {
        break
      }
    }
    # print(i)
    return(w1)
  } else {
    sapply(x, lambertWp_loc)
  }
}
SampleGumbel <- function(BivDistParams)  {
  theta1 = BivDistParams[1]
  theta2 = BivDistParams[2]
  delta = BivDistParams[3]
  
  delRecip = 1/delta
  
  if (delta == 1) {
    X1 = rexp(1, 1/theta1)
    X2 = rexp(1, 1/theta2)
    return(c(X1, X2))
  }
  
  # For X1 -------
  U1 = runif(1)
  X1 = -theta1*log(1 - U1)
  ### For X2 -----------
  h <- function(x) {
    a = (1/(delta - 1))
    b = (1/delta) - 1
    num = (1-x)
    den = exp(X1/theta1)*((X1/theta1)^b)
    return( (num/den)^a )
  }
  ginv <- function(x)  {
    a = 1 - delta
    b = lambertWp_loc((delta*x^delta)/a)
    c = ((a*b)/delta)^(1/delta)
    e = (X1/theta1)^(1/delta)
    return( theta2*((c - e)^delta) )
  }
  
  U2 = runif(1)
  # t1 = h(U2)
  # X2 = ginv(t1)
  X2 = ginv(h(U2))
  
  return(c(X1, X2))
}
GetGumb_Exp <- function(BivDistParams) {
  theta1 = BivDistParams[1]
  theta2 = BivDistParams[2]
  delta = BivDistParams[3]
  
  mu1 = theta1
  mu2 = theta2
  # return(sprintf("( %s , %s )", round(mu1,1), round(mu2,1) ))
  return(c(mu1, mu2))
}
########### Main functions ###########
DoEWMA_XieEtAl <- function(maxIter, h, lambda, ICparams)   {
  # There is no decorrelating, etc here and so this depends directly on the Gumbel parameters
  theta1 = ICparams[1]
  theta2 = ICparams[2]
  delta = ICparams[3]
  
  RunsCounter = TimeCounter = 0
  
  num = 2*(gamma(delta + 1))^2
  den = gamma(2*delta + 1)
  offDiag = ((num/den)-1)*theta1*theta2
  CovMat = matrix(c(theta1^2, offDiag, offDiag, theta2^2), nrow = 2, ncol = 2)
  SigmaInv = solve(CovMat)
  c = (2 - lambda)/lambda
  
  muVec = matrix(GetGumb_Exp(ICparams), nrow = 2, ncol = 1)
  zVec = matrix(rep(0,2), nrow = 2, ncol = 1)
  
  for (i in 1:maxIter)  {
    xVec = SampleGumbel(ICparams)
    xVec = matrix(xVec, nrow = 2, ncol = 1) # To turn xvals into a column vector as required for proper matrix mult
    parsedVals = ParseXvec(xVec)
    
    zVec = lambda*(xVec - muVec) + (1-lambda)*zVec
    E = c*t(zVec)%*%SigmaInv%*%zVec  # The EWMA statistic
    
    ## Recall that V != 2 ever in the Gumbel case. Hence since we don't have to worry about the X1 = X2 case we can just operate as below
    RunsCounter = RunsCounter + 2
    TimeCounter = TimeCounter + parsedVals[2] 
    
    
    if (E > h) {
      break
    }
  }
  return(RunsCounter)
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
  
  ICparams = job$ICparams
  theta1 = ICparams[1]
  theta2 = ICparams[2]
  delta = ICparams[3]
  out <- list(
    h          = signif(h, 4),
    lambda         = job$lambda,
    theta1     = theta1,
    theta2     = theta2,
    delta      = delta,
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

############# Main #############
setwd("path/to/folder/EWMA_Xie/") # Set to your path to this folder
### Tweaking params ####
lambdaVals = c(.1, .2) # Include whatever lambda values you want

#### Distribution Parameters ####

theta1_vec_Gumb = c(5,7.5,10,7.5,10,5,5,2.5,1, 5,7.5,10,7.5,10,5,5,2.5,1,  5,7.5,10,7.5,10,5,5,3.5,2.5, 5,7.5,10,7.5,10,5,5,3.5,2.5)
theta2_vec_Gumb = c(5,5,5,7.5,10,2.5,1,2.5,1, 5,5,5,7.5,10,2.5,1,2.5,1,   15,15,15,22.5,30,10.5,7.5,7.5,7.5, 15,15,15,22.5,30,10.5,7.5,7.5,7.5)
delta_vec = c(rep(1,9), rep(.5, 9), rep(1,9), rep(.5, 9))
scenario = c(rep(1,9), rep(2,9), rep(3,9), rep(4,9))
lmat_Gumb = matrix( c(theta1_vec_Gumb, theta2_vec_Gumb, delta_vec, scenario), nrow = 4, byrow = T)

Gumbel_tarARL0 = c(53.25, 62, 24.55, 26)


.ic_cols_Gumb <- c(1, 10, 19, 28)
h_init = 1

#### Create the jobs ###

# build the 20 jobs
ewma_pack <- make_EWMA_jobs(
  lmat_Gumb        = lmat_Gumb,
  Gumbel_tarARL0   = Gumbel_tarARL0,
  lambdaVals       = lambdaVals,
  h_init           = h_init,          
  fewerIter        = FALSE,
  logfile_dir      = "dataIC/IC_Rls/"
)

tab_ewma  <- ewma_pack$tab
jobs_ewma <- ewma_pack$jobs  

######## Create the Clims file ##########

clims_dir <- file.path("dataIC")
dir.create(clims_dir, showWarnings = FALSE, recursive = TRUE)

clims_file <- file.path(clims_dir, "Clims_EWMA.txt")

if (file.exists(clims_file)) {
  clims_tab <- fread(clims_file)
} else {
  clims_tab <- data.table(
    h          = numeric(),
    lambda  = numeric(),
    theta1     = numeric(),
    theta2     = numeric(),
    delta      = numeric(),
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

cl <- makeCluster(length(jobs_ewma))
clusterExport(
  cl,
  c(
    "DoEWMA_XieEtAl",
    "make_method_EWMA",
    "make_file_logger",
    "LogInterp",
    "make_ARL_policy",
    "make_clims_job",
    "IC_Rl_fname",
    "ObtainClims",
    "ParseXvec",
    "lambertWp_loc",
    "SampleGumbel",
    "GetGumb_Exp",
    "IC_rls_dir"
  ),
  envir = environment()
)

results <- parLapply(cl, jobs_ewma, ObtainClims)
stopCluster(cl)

new_rows <- rbindlist(lapply(seq_along(results), function(k) {
  x <- results[[k]]
  data.table(
    h          = x$h,
    lambda = x$lambda,
    theta1 = x$theta1,
    theta2 = x$theta2,
    delta = x$delta,
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



