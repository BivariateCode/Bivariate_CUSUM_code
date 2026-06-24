# Simulates the OC runlengths for the Adaptive CUSUM Chart proposed by Parakulam and Li
rm(list = ls(all.names = TRUE))
require(parallel)
library(data.table)

#### Helper functions ####
make_prior_df <- function(priorsMat, rhoVals) {
  data.frame(
    prior1    = priorsMat[,1],
    prior2    = priorsMat[,2],
    prior3    = priorsMat[,3],
    prior4    = priorsMat[,4],
    rho1      = rhoVals[1],
    rho2      = rhoVals[2]
  )
}
make_dist_obj <- function(name, GetExpVal, SamplePdf, Fx_v, Fy_xv, lmat, dist_pairs) {
  stopifnot(is.character(name), is.function(SamplePdf), is.function(Fx_v), is.function(Fy_xv))
  stopifnot(is.matrix(lmat), nrow(lmat) >= 4) # you said last row is scenario
  list(
    name = name,
    GetExpVal = GetExpVal,
    SamplePdf = SamplePdf,
    Fx_v = Fx_v,
    Fy_xv = Fy_xv,
    lmat = lmat,
    dist_pairs = as.data.frame(dist_pairs)
  )
}
extract_cols_info <- function(lmat) {
  stopifnot(is.matrix(lmat), ncol(lmat) == 36)
  
  scen_row <- lmat[nrow(lmat), ]
  scenario <- as.integer(scen_row)
  
  # IC columns are fixed positions
  ic_cols <- c(1, 10, 19, 28)
  is_ic <- seq_len(ncol(lmat)) %in% ic_cols
  
  list(
    scenario = scenario,
    is_ic = is_ic
  )
}
scenario_to_pairs <- function(dist_pairs) {
  stopifnot(is.matrix(dist_pairs) || is.data.frame(dist_pairs))
  dp <- as.data.frame(dist_pairs)
  stopifnot(all(c("p1", "ARL0") %in% names(dp)))
  stopifnot(nrow(dp) == 4)
  
  data.frame(
    scenario = 1:4,
    p1 = as.numeric(dp$p1),
    Tar_ARL0 = as.numeric(dp$ARL0),
    stringsAsFactors = FALSE
  )
}
get_h_for_pair <- function(h_table, Tar_ARL0, p1, prior1, prior2, prior3, prior4, rho1, rho2, functName = "Biv8_Stable_CUSUM") {
  ht <- as.data.frame(h_table)
  
  need <- c("h", "Tar_ARL0", "p1", "prior1", "prior2", "prior3", "prior4", "rho1", "rho2", "functName")
  miss <- setdiff(need, names(ht))
  if (length(miss) > 0) stop("h_table missing columns: ", paste(miss, collapse = ", "))
  
  w <- which(ht$Tar_ARL0 == Tar_ARL0 & ht$p1 == p1 & ht$prior1 == prior1 & ht$prior2 == prior2 & ht$prior3 == prior3 & ht$prior4 == prior4 & ht$rho1 == rho1 & ht$rho2 == rho2 & ht$functName == functName)
  if (length(w) == 0) stop("No h found for Tar_ARL0=", Tar_ARL0, ", p1=", p1, ", functName=", functName)
  if (length(w) > 1) stop("Multiple h rows found for Tar_ARL0=", Tar_ARL0, ", p1=", p1, ", functName=", functName)
  
  as.numeric(ht$h[w])
}
stable_fname <- function(prior1, prior2, prior3, prior4, rho1, rho2, p1) {
  sprintf(
    "Biv8_stable_pr_%.3g_%.3g_%.3g_%.3g_rho_%.3g_%.3g_p1_%.2f.rds",
    prior1, prior2, prior3, prior4,
    rho1, rho2,
    p1
  )
}
get_condCdfs_by_scenario <- function(priorParams, rhoVec, dist_obj) {
  
  p1s <- dist_obj$p1
  stopifnot(length(p1s) == 4)
  
  out <- vector("list", length(p1s))
  
  for (k in seq_along(p1s)) {
    
    stable_obj <- readRDS(
      file.path(
        "dataIC/CUSUM_dist",
        stable_fname(
          priorParams[1], priorParams[2],
          priorParams[3], priorParams[4],
          rhoVec[1], rhoVec[2],
          p1s[k]
        )
      )
    )
    
    tab <- stable_obj$stable_dist
    
    cdfs_stable <- lapply(1:8, function(i) {
      Ftemp <- ecdf(tab[, i])
      function(x) (Ftemp(x) - Ftemp(0)) / (1 - Ftemp(0))
    })
    
    out[[k]] <- list(
      tab          = tab,
      cdfs_stable  = cdfs_stable
    )
  }
  
  names(out) <- paste0("scenario_", seq_along(out))
  out
}

######## Simulation Helpers ########
run_one_oc_shift <- function(j,dist_obj,info, h_by_scen,condCdfs_by_scen, nRuns, maxIter, tau,priorParams,rhoVec) {
  
  lmat <- dist_obj$lmat
  last_row <- nrow(lmat)
  ic_cols <- c(1, 10, 19, 28)
  
  ## identify scenario
  scen <- info$scenario[j]
  h <- unname(h_by_scen[as.character(scen)])
  
  ## extract IC / OC params
  IC_col <- ic_cols[scen]
  ICparams <- as.numeric(lmat[1:(last_row - 1), IC_col])
  OCparams <- as.numeric(lmat[1:(last_row - 1), j])
  
  mean_label <- dist_obj$GetExpVal(OCparams)
  
  ## ---- NEW PART ----
  stable_obj_scen <- condCdfs_by_scen[[scen]]
  tab_stable      <- stable_obj_scen$tab
  cdfs_stable     <- stable_obj_scen$cdfs_stable
  ## ------------------
  
  # cat("h is: \n", h, "\n")
  # cat("tab_stable: \n")
  # print(tab_stable[1:5, 1:5])
  # f1 = cdfs_stable[[1]]
  # cat("f1(1) is: ", f1(1), "\n")
  
  sim <- simulate_OC_replicates(
    nRuns        = nRuns,
    maxIter      = maxIter,
    h            = h,
    tau          = tau,
    priorParams  = priorParams,
    rhoVec       = rhoVec,
    ICparams     = ICparams,
    OCparams     = OCparams,
    SamplePdf    = dist_obj$SamplePdf,
    Fx_v         = dist_obj$Fx_v,
    Fy_xv        = dist_obj$Fy_xv,
    condCdfs     = cdfs_stable,
    tab_stable   = tab_stable     # ← passed, not sampled yet
  )
  
  data.frame(
    scenario = scen,
    is_ic    = info$is_ic[j],
    mean     = mean_label,
    ARL      = round(mean(sim$RL), 1),
    sd_rl    = round(sd(sim$RL), 1),
    ATS      = round(mean(sim$TS), 1),
    sd_TS    = round(sd(sim$TS), 1),
    stringsAsFactors = FALSE
  )
}
simulate_OC_replicates <- function(nRuns, maxIter, h, tau, priorParams, rhoVec, ICparams, OCparams, SamplePdf, Fx_v, Fy_xv, condCdfs, tab_stable) {
  
  RLs <- numeric(0)
  TSs <- numeric(0)
  n_attempts <- 0L
  
  n_stable <- nrow(tab_stable)
  
  while (length(RLs) < nRuns) {
    
    ## ---- NEW: draw starting values from stable distribution ----
    idx <- sample.int(n_stable, size = 1L)
    startingValues <- as.numeric(tab_stable[idx, ])
    ## -----------------------------------------------------------
    
    out <- Biv8_CUSUM(
      maxIter        = maxIter,
      h              = h,
      tau            = tau,
      priorParams    = priorParams,
      rhoVec         = rhoVec,
      ICparams       = ICparams,
      OCparams       = OCparams,
      SamplePdf      = SamplePdf,
      Fx_v           = Fx_v,
      Fy_xv          = Fy_xv,
      condCdfs       = condCdfs,
      startingValues = startingValues
    )
    
    n_attempts <- n_attempts + 1L
    
    RL <- out[1]
    TS <- out[2]
    falseAlarm <- out[3]
    
    if (!falseAlarm) {
      RLs <- c(RLs, RL)
      TSs <- c(TSs, TS)
    }
  }
  
  list(
    RL = RLs,
    TS = TSs,
    n_attempts = n_attempts
  )
}
build_oc_table_one_dist <- function(dist_obj, h_table, condCdfs_by_scen, nRuns, maxIter, tau, priorParams, rhoVec, functName="Biv8_Stable_CUSUM", keep_ic_rows=TRUE, ncores=24) {
  
  stopifnot(is.list(dist_obj), is.matrix(dist_obj$lmat))
  
  lmat <- dist_obj$lmat
  
  info <- extract_cols_info(lmat)
  scen_pairs <- scenario_to_pairs(dist_obj$dist_pairs)
  
  ## scenario -> h lookup (now includes priors + rhos)
  scen_pairs$h <- mapply(
    FUN = function(a,p)
      get_h_for_pair(
        h_table  = h_table,
        Tar_ARL0 = a,
        p1       = p,
        prior1   = priorParams[1],
        prior2   = priorParams[2],
        prior3   = priorParams[3],
        prior4   = priorParams[4],
        rho1     = rhoVec[1],
        rho2     = rhoVec[2],
        functName = functName
      ),
    a = scen_pairs$Tar_ARL0,
    p = scen_pairs$p1
  )
  
  h_by_scen <- setNames(scen_pairs$h, scen_pairs$scenario)
  
  rows <- parallel::mclapply(
    X        = seq_len(ncol(lmat)),
    FUN      = run_one_oc_shift,
    dist_obj = dist_obj,
    info     = info,
    h_by_scen = h_by_scen,
    condCdfs_by_scen = condCdfs_by_scen,
    nRuns    = nRuns,
    maxIter  = maxIter,
    tau      = tau,
    priorParams = priorParams,
    rhoVec      = rhoVec,
    mc.cores = ncores
  )
  
  tab <- do.call(rbind, rows)
  
  if (!keep_ic_rows) {
    tab <- tab[!tab$is_ic, , drop = FALSE]
  }
  
  tab <- tab[order(tab$scenario), ]
  rownames(tab) <- NULL
  tab
}
####### Logging: ##########
oc_fname <- function(prior1,prior2,prior3,prior4,rho1,rho2,tau,arl500) {
  if (arl500) {
    sprintf(
      "ARL500_OC_pr_%.3g_%.3g_%.3g_%.3g_rho_%.3g_%.3g_tau%s.txt",
      prior1, prior2, prior3, prior4,
      rho1, rho2,
      tau
    )
  } else  {
    sprintf(
      "OC_pr_%.3g_%.3g_%.3g_%.3g_rho_%.3g_%.3g_tau%s.txt",
      prior1, prior2, prior3, prior4,
      rho1, rho2,
      tau
    )
  }
}
###### Distribution Related Functions #####
Fx_v_Gumb <- function(x, BivDistParams)  {
  theta1 = BivDistParams[1]
  theta2 = BivDistParams[2]
  delta = BivDistParams[3]
  
  delRecip = 1/delta
  a = (1/theta1)^delRecip + (1/theta2)^delRecip
  return(1 - exp(-x*(a^delta)))
}
Fy_xv_Gumb <- function(x, y, v, BivDistParams)  {
  theta1 = BivDistParams[1]
  theta2 = BivDistParams[2]
  delta = BivDistParams[3]
  
  delRecip = 1/delta
  a = 1 - delRecip
  temp = (1/theta1)^delRecip + (1/theta2)^delRecip
  beta = temp^delta
  psi = ((x/theta1)^delRecip + (x/theta2)^delRecip) ^(delta-1)
  
  if (v == 0) {
    Cxy = (x/theta1)^(delRecip) + (y/theta2)^(delRecip)
  } else  {
    Cxy = (x/theta2)^(delRecip) + (y/theta1)^(delRecip)
  }
  t1 = Cxy^(delta-1)
  t2 = Cxy^(delta)
  
  num = exp(-t2)*t1
  den = psi*exp(-beta*x)
  
  return(1 - (num/den))
}
Fx_v_MOBE <- function(x, BivDistParams)  { # invariant of v
  l1 = BivDistParams[1]
  l2 = BivDistParams[2]
  l3 = BivDistParams[3]
  
  Lambda = l1 + l2 + l3
  return(1 - exp(-Lambda*x) )
}
Fy_xv_MOBE <- function(x, y, v, BivDistParams) {
  l1 = BivDistParams[1]
  l2 = BivDistParams[2]
  l3 = BivDistParams[3]
  
  d = x - y
  if (v == 0) {
    return(1 - exp(l2*d + l3*d))
  } else if (v == 1)  {
    return(1 - exp(l1*d + l3*d))
  }
}
Fx_v_MOBW <- function(x, BivDistParams)  { # invariant of v
  l1 = BivDistParams[1]
  l2 = BivDistParams[2]
  l3 = BivDistParams[3]
  eta = BivDistParams[4]
  Lambda = l1 + l2 + l3
  return(1 - exp(-Lambda*(x^eta)) )
}
Fy_xv_MOBW <- function(x, y, v, BivDistParams) {
  l1 = BivDistParams[1]
  l2 = BivDistParams[2]
  l3 = BivDistParams[3]
  eta = BivDistParams[4]
  d = (x^eta) - (y^eta)
  if (v == 0) {
    return(1 - exp(l2*d + l3*d))
  } else if (v == 1)  {
    return(1 - exp(l1*d + l3*d))
  }
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
SampleMOBE <- function(BivDistParams) { # Marshal Olkin Bivariate Exponential sampling
  l1 = BivDistParams[1]
  l2 = BivDistParams[2]
  l3 = BivDistParams[3]
  
  Z1 = rexp(1, l1)
  Z2 = rexp(1, l2)
  if (l3 == 0)  {
    Z3 = Inf
  } else  {
    Z3 = rexp(1, l3)
  }
  X1 = min(Z1, Z3)
  X2 = min(Z2, Z3)
  return(c(X1, X2))
}
SampleMOBW <- function(BivDistParams)  {
  l1 = BivDistParams[1]
  l2 = BivDistParams[2]
  l3 = BivDistParams[3]
  eta = BivDistParams[4]
  
  b1 = (1/l1)^(1/eta)
  b2 = (1/l2)^(1/eta)
  b3 = (1/l3)^(1/eta)
  Z1 = rweibull(1, eta, b1)
  Z2 = rweibull(1, eta, b2)
  if (l3 == 0)  {
    Z3 = Inf
  } else  {
    Z3 = rweibull(1, eta, b3)
  }
  X1 = min(Z1, Z3)
  X2 = min(Z2, Z3)
  return(c(X1, X2))
}

# LCL and UCL as given in Bivariate Paper
MOBE_Clims_1 <- function(l1, l2, l3, alpha)  {
  Lambda = l1 + l2 + l3
  LCL = -log(1 - (alpha/2))/Lambda
  UCL = -log(alpha/2)/Lambda
  return(c(LCL, UCL))
}
MOBE_Clims_2 <- function(x, v, l1, l2, l3, alpha)  {
  if (v == 0) {
    c = l2 + l3
  } else if (v == 1)  {
    c = l1 + l3
  }
  
  LCL = x - log(1 - (alpha/2))/c
  UCL = x - log(alpha/2)/c
  return(c(LCL, UCL))
}
MOBW_Clims_1 <- function(l1, l2, l3, eta, alpha)  {
  Lambda = l1 + l2 + l3
  LCL = (-log(1 - (alpha/2))/Lambda)^(1/eta)
  UCL = (-log(alpha/2)/Lambda)^(1/eta)
  return(c(LCL, UCL))
}
MOBW_Clims_2 <- function(x, v, l1, l2, l3, eta, alpha)  {
  if (v == 0) {
    c = l2 + l3
  } else if (v == 1)  {
    c = l1 + l3
  }
  
  LCL = (x^eta - log(1 - (alpha/2))/c)^(1/eta)
  UCL = (x^eta - log(alpha/2)/c)^(1/eta)
  return(c(LCL, UCL))
}
Gumb_Clims_1 <- function(theta1, theta2, delta, alpha)  {
  delRecip = 1/delta
  a = 1 - delRecip
  temp = (1/theta1)^delRecip + (1/theta2)^delRecip
  beta = temp^delta
  
  
  LCL = -log(1- (alpha/2))/beta
  UCL = -log(alpha/2)/beta
  return(c(LCL, UCL))
}

MOBE_ConvertARL <- function(BivDistParams, ARL) {
  l1 = BivDistParams[1]
  l2 = BivDistParams[2]
  l3 = BivDistParams[3]
  
  # cat("The ARL is :", ARL, "\n")
  
  L = l1 + l2 + l3
  t1 = (l2/(L^2)) + (l2/(L*(l1 + l3))) + (l1/(L^2)) + (l1/(L*(l2 + l3)))
  TBE = .5*t1 + l3/(L^2)
  # ARL = ATS/TBE
  ATS = ARL*TBE
  return(ATS)
}

## The functions below find the probability of Obs1 (namely .5+P(V == 2)), Prob V == 1 conditional that V != 2
MOBE_ProbV1_cond <- function(BivDistParams) {
  l1 = BivDistParams[1]
  l2 = BivDistParams[2]
  l3 = BivDistParams[3]
  
  L = l1 + l2 + l3
  gamma1 = (l1 + l2)/L
  gamma2 = l3/L
  t1 = 1/(2*gamma1 + gamma2)
  t2 = l2/(l1 + l2)
  v = c(t1, t2)
  names(v) = c("Obs1_prob", "V1_Prob_obs2")
  return(v)
}
MOBW_ProbV1_cond <- function(BivDistParams) {
  l1 = BivDistParams[1]
  l2 = BivDistParams[2]
  l3 = BivDistParams[3]
  
  L = l1 + l2 + l3
  gamma1 = (l1 + l2)/L
  gamma2 = l3/L
  t1 = 1/(2*gamma1 + gamma2)
  t2 = l2/(l1 + l2)
  v = c(t1, t2)
  names(v) = c("Obs1_prob", "V1_Prob_obs2")
  return(v)
}
Gumb_ProbV1_cond <- function(BivDistParams) {
  # In the Gumbel case we don't have the V == 2 case
  theta1 = BivDistParams[1]
  theta2 = BivDistParams[2]
  delta = BivDistParams[3]
  t1 = theta1^(-1/delta)
  t2 = ((1/theta1)^(1/delta)) + ((1/theta2)^(1/delta))
  return(1 - (t1/t2))
}

GetMOBE_Exp <- function(BivDistParams) {
  l1 = BivDistParams[1]
  l2 = BivDistParams[2]
  l3 = BivDistParams[3]
  
  mu1 = round( 1/(l1 + l3), 1)
  mu2 = round( 1/(l2 + l3), 1)
  return( paste("(", mu1, ",", mu2, ")") )
}
GetMOBW_Exp <- function(BivDistParams) {
  l1 = BivDistParams[1]
  l2 = BivDistParams[2]
  l3 = BivDistParams[3]
  eta = BivDistParams[4]
  
  r = 1/eta
  mu1 = gamma(1 + r)*((1/(l1 + l3))^r)
  mu2 = gamma(1 + r)*((1/(l2 + l3))^r)
  return(sprintf("( %s , %s )", round(mu1,1), round(mu2,1) ))
}
GetGumb_Exp <- function(BivDistParams) {
  theta1 = BivDistParams[1]
  theta2 = BivDistParams[2]
  delta = BivDistParams[3]
  
  mu1 = theta1
  mu2 = theta2
  return(sprintf("( %s , %s )", round(mu1,1), round(mu2,1) ))
}
#### Core functions ######
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
Biv8_CUSUM <- function(maxIter, h , tau, priorParams, rhoVec, ICparams, OCparams, SamplePdf, Fx_v, Fy_xv, condCdfs, startingValues) {
  # p1 is P(V = 1 | V != 2)
  
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
  
  GetEstimates <- function(Cvec, Svec, Tvec, lamPosInd, lamNegInd, l_defaultVals, Z, initVals = F)  {
    posInd = (Cvec > 0)
    lambda_vec = rep(NA, 8)
    
    if (initVals == F)  {
      Svec[posInd] = Svec[posInd] + Z
      Svec[!posInd] = 0
      
      Tvec[posInd] = Tvec[posInd] + 1
      Tvec[!posInd] = 0
    }
    
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
  
  # ----- START REFACTORED ARL SCRIPT LOOP -----
  
  # --- Initializations for a SINGLE ARL RUN ---
  # V_prev = -1 # We'll use 'L_type_of_prev_obs' and 'Z_val_of_prev_obs' instead
  Z_val_of_prev_obs = NA
  L_type_of_prev_obs = NA # Will be 1, 2, or 3
  
  RunsCounter = 0       # Counts processed R_t values
  TimeCounter = 0
  pair_idx = 0          # To count X_i pairs
  
  # ARL1 specific
  RunsAtTau = 0
  TimeAtTau = 0
  
  # State for managing pairs and their observations
  is_first_obs_of_pair = TRUE
  V_of_current_pair = -1
  # To store parts of X_i needed for the second observation
  # (parsedVals[1] for X_ord1, parsedVals[2] for X_ord2)
  X_ord1_current_pair = NA
  Y_ord2_current_pair = NA # Only strictly needed if V != 2
  
  for (rt_iter_idx in 1:maxIter) { # This loop processes one R_t at a time
    # --- 1. Update Estimates (based on the *previous* R_t) ---
    if (!is.na(Z_val_of_prev_obs)) {
      if (L_type_of_prev_obs == 1) {
        # Update for L=1 type observation
        L_res1 = GetEstimates(Cvec, Svec1, Tvec1, lamPosInd1, lamNegInd1, l1_defaultVals, Z_val_of_prev_obs)
        L_res2 = ResetEstimates(Cvec, Svec2, Tvec2, lamPosInd2, lamNegInd2, l2_defaultVals)
        L_res3 = ResetEstimates(Cvec, Svec3, Tvec3, lamPosInd3, lamNegInd3, l3_defaultVals)
        Svec1 = L_res1[[1]]; Tvec1 = L_res1[[2]]; lambda_vec1 = L_res1[[3]]
        Svec2 = L_res2[[1]]; Tvec2 = L_res2[[2]]; lambda_vec2 = L_res2[[3]]
        Svec3 = L_res3[[1]]; Tvec3 = L_res3[[2]]; lambda_vec3 = L_res3[[3]]
      } else if (L_type_of_prev_obs == 2) {
        # Update for L=2 type observation
        L_res1 = ResetEstimates(Cvec, Svec1, Tvec1, lamPosInd1, lamNegInd1, l1_defaultVals)
        L_res2 = GetEstimates(Cvec, Svec2, Tvec2, lamPosInd2, lamNegInd2, l2_defaultVals, Z_val_of_prev_obs)
        L_res3 = ResetEstimates(Cvec, Svec3, Tvec3, lamPosInd3, lamNegInd3, l3_defaultVals)
        Svec1 = L_res1[[1]]; Tvec1 = L_res1[[2]]; lambda_vec1 = L_res1[[3]]
        Svec2 = L_res2[[1]]; Tvec2 = L_res2[[2]]; lambda_vec2 = L_res2[[3]]
        Svec3 = L_res3[[1]]; Tvec3 = L_res3[[2]]; lambda_vec3 = L_res3[[3]]
      } else if (L_type_of_prev_obs == 3) {
        # Update for L=3 type observation
        L_res1 = ResetEstimates(Cvec, Svec1, Tvec1, lamPosInd1, lamNegInd1, l1_defaultVals)
        L_res2 = ResetEstimates(Cvec, Svec2, Tvec2, lamPosInd2, lamNegInd2, l2_defaultVals)
        L_res3 = GetEstimates(Cvec, Svec3, Tvec3, lamPosInd3, lamNegInd3, l3_defaultVals, Z_val_of_prev_obs)
        Svec1 = L_res1[[1]]; Tvec1 = L_res1[[2]]; lambda_vec1 = L_res1[[3]]
        Svec2 = L_res2[[1]]; Tvec2 = L_res2[[2]]; lambda_vec2 = L_res2[[3]]
        Svec3 = L_res3[[1]]; Tvec3 = L_res3[[2]]; lambda_vec3 = L_res3[[3]]
      }
    }
    # --- 2. Generate Current R_t (Z_current_rt) and its L-type ---
    Z_current_rt = NA
    L_type_current_rt = NA
    time_increment_current_rt = NA
    
    if (is_first_obs_of_pair) {
      pair_idx = pair_idx + 1
      
      # Check for tau *before* sampling (tau is pair-based)
      # This is for ARL1 calculations: record state *before* the first OC pair's *first* R_t
      if (pair_idx == tau & tau > 1) { # if tau is 1, handled by initial value where RunsAtTau = 1 opposed to setting it to 0
        RunsAtTau = RunsCounter # State *before* this OC observation
        TimeAtTau = TimeCounter
      }
      
      current_params = if (pair_idx < tau) ICparams else OCparams
      xvals = SamplePdf(current_params)
      parsedVals = ParseXvec(xvals) # Assuming: [X_ord(1), X_ord(2), V]
      
      X_ord1_current_pair = parsedVals[1]
      Y_ord2_current_pair = parsedVals[2] # Store even if V=2, for consistency
      V_of_current_pair = parsedVals[3]
      
      U1 = Fx_v(X_ord1_current_pair, ICparams) # Transform using IC params
      Z_current_rt = qexp(U1, rate = 1)
      L_type_current_rt = 1
      time_increment_current_rt = X_ord1_current_pair
      
    } else { # Processing the second observation of the pair
      # This means V_of_current_pair must have been 0 or 1
      U2 = Fy_xv(x = X_ord1_current_pair, y = Y_ord2_current_pair, v = V_of_current_pair, ICparams) # IC params
      Z_current_rt = qexp(U2, rate = 1)
      L_type_current_rt = if (V_of_current_pair == 0) 2 else 3
      time_increment_current_rt = Y_ord2_current_pair - X_ord1_current_pair
    }
    
    # Handle Inf values (as in your original ARL script)
    # This needs to be adapted if Z_current_rt is Inf, it means this R_t causes a stop.
    if (is.infinite(Z_current_rt)) {
      RunsCounter = RunsCounter + 1 # Count this problematic R_t
      TimeCounter = TimeCounter + time_increment_current_rt
      # How 'falseAlarm' is determined here depends on 'pair_idx' vs 'tau'
      # If this structure is used, an infinite Z value usually means a problem with parameters or transformations.
      # For now, we'll assume it leads to a break. The return logic will handle ARL.
      break # Exit the rt_iter_idx loop
    }
    
    # --- 3. Update CUSUM ---
    current_lambda_for_cusum = NA
    if (L_type_current_rt == 1) {
      current_lambda_for_cusum = lambda_vec1
    } else if (L_type_current_rt == 2) {
      current_lambda_for_cusum = lambda_vec2
    } else { # L_type_current_rt == 3
      current_lambda_for_cusum = lambda_vec3
    }
    Cvec = pmax(0, Cvec + log(current_lambda_for_cusum) + (1 - current_lambda_for_cusum) * Z_current_rt)
    
    # --- Bookkeeping for this R_t ---
    RunsCounter = RunsCounter + 1
    TimeCounter = TimeCounter + time_increment_current_rt
    
    # Store details for next iteration's estimate update
    Z_val_of_prev_obs = Z_current_rt
    L_type_of_prev_obs = L_type_current_rt
    
    # --- 4. Update State for Next R_t Generation ---
    if (is_first_obs_of_pair) {
      if (V_of_current_pair == 2) {
        is_first_obs_of_pair = TRUE # Next R_t will be from a new pair
      } else {
        is_first_obs_of_pair = FALSE # Next R_t will be second obs of current pair
      }
    } else { # Was processing second obs of a pair
      is_first_obs_of_pair = TRUE # Next R_t will be from a new pair
    }
    
    # --- 5. Signaling Check (occurs after every R_t) ---
    pvals = 1 - vapply(seq_along(condCdfs), function(j) condCdfs[[j]](Cvec[j]), numeric(1))
    cusumStndzd = -log(pvals)

    
    if (max(cusumStndzd) > h) {
      break # Exit the rt_iter_idx loop (Signal)
    }
    
  } # End of rt_iter_idx loop (processing one R_t at a time)
  
  ####### My post loop Calculations #######
  
  falseAlarm = FALSE
  if (pair_idx < tau) { # Loop broke while still processing IC pairs
    falseAlarm = TRUE
  }
  # For IC ARL0 (OCparams == ICparams), we want the total RunsCounter
  # Your original script reset RunsAtTau for IC case:
  is_IC_case = if (sum(abs(ICparams - OCparams)) < 1e-6) TRUE else FALSE
  if (is_IC_case) {
    ARL_val = RunsCounter
    ATS_val = TimeCounter
    falseAlarm = F
  } else {
    # OC case: ARL1 is number of OC observations until signal
    # RunsCounter is total R_t. RunsAtTau is R_t *before* first OC observation.
    ARL_val = RunsCounter - RunsAtTau
    ATS_val = TimeCounter - TimeAtTau
  }
  
  
  return(c(ARL_val, ATS_val, falseAlarm))
}


############# Main #############
setwd("/pathToFolder/AdaptiveCUSUM_SimulationCode") # Path to this folder
### Tweaking params ####
totSim = 1e4
ncpus = 24
maxIter = 8000
tau = 1 # The first OC data point. So tau = 1 means no IC data points
arl500 = F

#### Distribution Parameters ####

l1_vec_MOBE = c(.2, .133, .1, .133, .1, .2, .2, .4, 1, .164, .103, .073, .109, .081, .145, .891, .327, .818, .2, .133, .1, .133, .1, .2, .2, .286, .4, .176, .115, .085, .117, .088, .173, .170, .248, .352)
l2_vec_MOBE = c(.2, .2, .2, .133, .1, .4, 1, .4, 1, .164, .17, .173, .109, .081, .345, .091, .327, .818, .067, .067, .067, .044, .033, .095, .133, .133, .133, .042, .048, .052, .028, .021, .068, .103, .095, .085)
l3_vec_MOBE = c(rep(0,9), .036, .030, .027, .024, .018, .055, .109, .073, .182, rep(0,9), .024, .018, .015, .016, .012, .027, .030, .038, .048)
scenario = c(rep(1,9), rep(2,9), rep(3,9), rep(4,9))
lmat_MOBE = matrix( c(l1_vec_MOBE, l2_vec_MOBE, l3_vec_MOBE, scenario), nrow = 4, byrow = T)

theta1_vec_Gumb = c(5,7.5,10,7.5,10,5,5,2.5,1, 5,7.5,10,7.5,10,5,5,2.5,1,  5,7.5,10,7.5,10,5,5,3.5,2.5, 5,7.5,10,7.5,10,5,5,3.5,2.5)
theta2_vec_Gumb = c(5,5,5,7.5,10,2.5,1,2.5,1, 5,5,5,7.5,10,2.5,1,2.5,1,   15,15,15,22.5,30,10.5,7.5,7.5,7.5, 15,15,15,22.5,30,10.5,7.5,7.5,7.5)
delta_vec = c(rep(1,9), rep(.5, 9), rep(1,9), rep(.5, 9))
scenario = c(rep(1,9), rep(2,9), rep(3,9), rep(4,9))
lmat_Gumb = matrix( c(theta1_vec_Gumb, theta2_vec_Gumb, delta_vec, scenario), nrow = 4, byrow = T)

l1_vec_MOBW = c(0.0314, 0.0140, 0.0079, 0.0140, 0.0079, 0.0314, 0.0314, 0.1257, 0.7854, 0.0257, 0.0098, 0.0043, 0.0114, 0.0064, 0.0171, 0.0107, 0.1028, 0.6426, 0.0314, 0.0140, 0.0079, 0.0140, 0.0079, 0.0314, 0.0314, 0.0641, 0.1257, 0.0282, 0.0124, 0.0068, 0.0126, 0.0070, 0.0279, 0.0273, 0.0570, 0.1130)
l2_vec_MOBW = c(0.0314, 0.0314, 0.0314, 0.0140, 0.0079, 0.1257, 0.7854, 0.1257, 0.7854, 0.0257, 0.0273, 0.0278, 0.0114, 0.0064, 0.1114, 0.1756, 0.1028, 0.6426, 0.0035, 0.0035, 0.0035, 0.0016, 0.0009, 0.0071, 0.0140, 0.0140, 0.0140, 3.17e-04, 1.90e-03, 2.46e-03, 1.41e-04, 7.93e-05, 3.62e-03, 9.83e-03, 6.86e-03, 1.270e-03)
l3_vec_MOBW = c(rep(0,9), 0.0057, 0.0041, 0.0036, 0.0025, 0.0014, 0.0143, 0.0207, 0.0228, 0.1428, rep(0,9), 0.0032, 0.0016, 0.0010, 0.0014, 0.0008, 0.0035, 0.0041, 0.0070, 0.0127)
eta = rep(2,36)
scenario = c(rep(1,9), rep(2,9), rep(3,9), rep(4,9))
lmat_MOBW = matrix( c(l1_vec_MOBW, l2_vec_MOBW, l3_vec_MOBW, eta, scenario), nrow = 5, byrow = T)

##### Priors,  ARL and p1 combos #######

priorsMat1 <- rbind(
  c(22.05, 21, 9.5, 10)  # This is the prior set that was used in the paper
)
rhoVals1 = c(1.05, .95)

if (arl500) {
  MOBE_tarARL0 = rep(500, 4)
  MOBW_tarARL0 = rep(500, 4)
  Gumbel_tarARL0 = rep(500, 4)
} else  {
  MOBE_tarARL0 = c(53.25, 53.25, 24.55, 24.55)
  MOBW_tarARL0 = c(62, 60, 26, 26)
  Gumbel_tarARL0 = c(53.25, 62, 24.55, 26)
}

unique_p1s = c(.1, .2, .5)

## These here are the true p1s (aka P(V = 1| V != 2) but NOT used
# MOBE_p1s = c(.5, .5, .25, .2)
# MOBW_p1s = c(.5, .5, .1, .01)
# Gumbel_p1s = c(.5, .5, .25, .1)

## We use these p1s instead aka we approximate to reduce the number of clims and stable dists required and it gives same results
MOBE_p1s = c(.5, .5, .2, .2)
MOBW_p1s = c(.5, .5, .1, .1)
Gumbel_p1s = c(.5, .5, .2, .1)

MOBE_pairs = as.data.frame(cbind(p1 = MOBE_p1s,   ARL0 = MOBE_tarARL0))
MOBW_pairs = cbind(p1 = MOBW_p1s,   ARL0 = MOBW_tarARL0)
Gumb_pairs = cbind(p1 = Gumbel_p1s, ARL0 = Gumbel_tarARL0)

MOBE_pairs
Gumb_pairs
MOBW_pairs

pairs <- rbind(
  cbind(p1 = MOBE_p1s,   ARL0 = MOBE_tarARL0),
  cbind(p1 = MOBW_p1s,   ARL0 = MOBW_tarARL0),
  cbind(p1 = Gumbel_p1s, ARL0 = Gumbel_tarARL0)
)
unique_pairs <- unique(pairs)
unique_pairs

prior_df <- rbind(
  make_prior_df(priorsMat1, rhoVals1)
  # make_prior_df(priorsMat2, rhoVals2)
)
prior_df

dist_MOBE <- make_dist_obj(name = "MOBE", GetExpVal = GetMOBE_Exp, SamplePdf = SampleMOBE, Fx_v = Fx_v_MOBE, Fy_xv = Fy_xv_MOBE, lmat = lmat_MOBE, dist_pairs = MOBE_pairs)
dist_MOBW <- make_dist_obj(name = "MOBW", GetExpVal = GetMOBW_Exp, SamplePdf = SampleMOBW, Fx_v = Fx_v_MOBW, Fy_xv = Fy_xv_MOBW, lmat = lmat_MOBW, dist_pairs = MOBW_pairs)
dist_Gumbel <- make_dist_obj(name = "Gumbel", GetExpVal = GetGumb_Exp, SamplePdf = SampleGumbel, Fx_v = Fx_v_Gumb, Fy_xv = Fy_xv_Gumb, lmat = lmat_Gumb, dist_pairs = Gumb_pairs)


clims = read.table("dataIC/Clims_CUSUM.txt", header = T)

########## Obtain OC runlengths ##########

dist_oc_dir = "OC_data/MOBE"
for (k in seq_len(nrow(prior_df))) {
  priorParams <- as.numeric(prior_df[k, c("prior1","prior2","prior3","prior4")])
  rhoVec      <- as.numeric(prior_df[k, c("rho1","rho2")])
  
  condCdfs_by_scen <- get_condCdfs_by_scenario(
    priorParams = priorParams,
    rhoVec      = rhoVec,
    dist_obj    = dist_MOBE$dist_pairs   # or just dist name if you prefer
  )
  
  tab_MOBE <- build_oc_table_one_dist(
    dist_obj          = dist_MOBE,
    h_table           = clims,
    priorParams       = priorParams,
    rhoVec            = rhoVec,
    condCdfs_by_scen  = condCdfs_by_scen,
    nRuns             = totSim,
    maxIter           = maxIter,
    tau               = tau,
    functName         = "Biv8_Stable_CUSUM",
    ncores            = ncpus
  )
  
  fname <- oc_fname(
    priorParams[1], priorParams[2],
    priorParams[3], priorParams[4],
    rhoVec[1], rhoVec[2], tau = tau, arl500 = arl500
  )
  
  write.table(
    tab_MOBE,
    file = file.path(dist_oc_dir, fname),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

dist_oc_dir = "OC_data/MOBW"
for (k in seq_len(nrow(prior_df))) {
  priorParams <- as.numeric(prior_df[k, c("prior1","prior2","prior3","prior4")])
  rhoVec      <- as.numeric(prior_df[k, c("rho1","rho2")])

  condCdfs_by_scen <- get_condCdfs_by_scenario(
    priorParams = priorParams,
    rhoVec      = rhoVec,
    dist_obj    = dist_MOBW$dist_pairs   # or just dist name if you prefer
  )

  tab_MOBW <- build_oc_table_one_dist(
    dist_obj          = dist_MOBW,
    h_table           = clims,
    priorParams       = priorParams,
    rhoVec            = rhoVec,
    condCdfs_by_scen  = condCdfs_by_scen,
    nRuns             = totSim,
    maxIter           = maxIter,
    tau               = tau,
    functName         = "Biv8_Stable_CUSUM",
    ncores            = ncpus
  )

  fname <- oc_fname(
    priorParams[1], priorParams[2],
    priorParams[3], priorParams[4],
    rhoVec[1], rhoVec[2], tau = tau, arl500 = arl500
  )

  write.table(
    tab_MOBW,
    file = file.path(dist_oc_dir, fname),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

dist_oc_dir = "OC_data/Gumbel"
for (k in seq_len(nrow(prior_df))) {
  priorParams <- as.numeric(prior_df[k, c("prior1","prior2","prior3","prior4")])
  rhoVec      <- as.numeric(prior_df[k, c("rho1","rho2")])

  condCdfs_by_scen <- get_condCdfs_by_scenario(
    priorParams = priorParams,
    rhoVec      = rhoVec,
    dist_obj    = dist_Gumbel$dist_pairs   # or just dist name if you prefer
  )

  tab_Gumbel <- build_oc_table_one_dist(
    dist_obj          = dist_Gumbel,
    h_table           = clims,
    priorParams       = priorParams,
    rhoVec            = rhoVec,
    condCdfs_by_scen  = condCdfs_by_scen,
    nRuns             = totSim,
    maxIter           = maxIter,
    tau               = tau,
    functName         = "Biv8_Stable_CUSUM",
    ncores            = ncpus
  )

  fname <- oc_fname(
    priorParams[1], priorParams[2],
    priorParams[3], priorParams[4],
    rhoVec[1], rhoVec[2], tau = tau, arl500 = arl500
  )

  write.table(
    tab_Gumbel,
    file = file.path(dist_oc_dir, fname),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

