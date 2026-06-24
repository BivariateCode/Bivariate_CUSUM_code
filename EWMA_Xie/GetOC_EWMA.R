rm(list = ls(all.names = TRUE))
library(parallel)
library(data.table)

get_h_for_pair <- function(h_table, Tar_ARL0, lambda, theta1, theta2, delta, functName = "EWMA_XieEtAl") {
  ht <- as.data.frame(h_table)
  
  need <- c("h", "Tar_ARL0", "lambda", "theta1", "theta2", "delta", "functName")
  miss <- setdiff(need, names(ht))
  if (length(miss) > 0) stop("h_table missing columns: ", paste(miss, collapse = ", "))
  
  w <- which(ht$Tar_ARL0 == Tar_ARL0 & ht$lambda == lambda & ht$theta1 == theta1 & ht$theta2 == theta2 & ht$delta == delta & ht$functName == functName)
  if (length(w) == 0) stop("No h found for Tar_ARL0=", Tar_ARL0, ", lambda=", lambda, ", functName=", functName)
  if (length(w) > 1) stop("Multiple h rows found for Tar_ARL0=", Tar_ARL0, ", lambda=", lambda, ", functName=", functName)
  
  as.numeric(ht$h[w])
}
get_ICparams_for_scenario <- function(lmat_Gumb, scenario, ic_cols = .ic_cols_Gumb) {
  lmat_Gumb[1:3, ic_cols[scenario]]
}
get_scenario_of_col <- function(col_idx) {
  ((col_idx - 1) %/% 9) + 1
}
GetGumb_Exp_str <- function(BivDistParams) {
  theta1 = BivDistParams[1]
  theta2 = BivDistParams[2]
  delta = BivDistParams[3]
  
  mu1 = theta1
  mu2 = theta2
  return(sprintf("( %s , %s )", round(mu1,1), round(mu2,1) ))
  # return(c(mu1, mu2))
}
# Main: run all lambdas for one scenario (IC fixed by scenario), across that scenario's 9 columns
run_EWMA_OC_one_scenario <- function(scenario, lmat_Gumb, clims, Gumbel_tarARL0, lambdaVals, totSim, maxIter, functName = "EWMA_XieEtAl") {
  
  # IC settings for this scenario
  ICparams  <- get_ICparams_for_scenario(lmat_Gumb, scenario, ic_cols = .ic_cols_Gumb)
  targetARL0 <- Gumbel_tarARL0[scenario]
  
  # columns belonging to this scenario
  cols <- ((scenario - 1) * 9 + 1):(scenario * 9)
  
  out_list <- vector("list", length(cols) * length(lambdaVals))
  k <- 1L
  
  for (col_idx in cols) {
    OCparams <- lmat_Gumb[1:3, col_idx]
    mean_str <- GetGumb_Exp_str(OCparams)  # e.g. "( 5 , 5 )"
    is_ic    <- (col_idx == .ic_cols_Gumb[scenario])
    
    # unpack theta values for h lookup (use the IC params for theta1/theta2/delta)
    th1 <- ICparams[1]
    th2 <- ICparams[2]
    del <- ICparams[3]
    
    for (lambda in lambdaVals) {
      
      # get corresponding h for this scenario + lambda + IC params
      h <- get_h_for_pair(
        h_table  = clims,
        Tar_ARL0 = targetARL0,
        lambda   = lambda,
        theta1   = th1,
        theta2   = th2,
        delta    = del,
        functName = functName
      )
      
      # run totSim replicates; each returns c(total_runs_until_alarm, time_to_signal)
      mat <- replicate(
        totSim,
        DoEWMA_XieEtAl_OC(
          maxIter  = maxIter,
          h        = h,
          lambda   = lambda,
          ICparams = ICparams,
          OCparams = OCparams
        )
      )
      
      # replicate() will return a 2 x totSim matrix if your function returns length-2 numeric
      mat <- as.matrix(mat)
      stopifnot(nrow(mat) == 2)
      
      rl <- as.numeric(mat[1, ])
      ts <- as.numeric(mat[2, ])
      
      out_list[[k]] <- data.frame(
        scenario = scenario,
        col_idx  = col_idx,
        lambda   = lambda,
        is_ic    = is_ic,
        mean     = mean_str,
        ARL      = round(mean(rl),1),
        sd_rl    = round(sd(rl),1),
        ATS      = round(mean(ts),1),
        sd_TS    = round(sd(ts),1),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  
  do.call(rbind, out_list)
}
# Driver: parallelize over scenarios only
run_EWMA_OC_all <- function(lmat_Gumb, clims, Gumbel_tarARL0, lambdaVals, totSim, maxIter = 8000, ncpus = 24, functName = "EWMA_XieEtAl") {
  
  scenarios <- 1:4
  mc <- min(length(scenarios), ncpus)
  
  res_list <- mclapply(
    scenarios,
    run_EWMA_OC_one_scenario,
    lmat_Gumb      = lmat_Gumb,
    clims          = clims,
    Gumbel_tarARL0 = Gumbel_tarARL0,
    lambdaVals     = lambdaVals,
    totSim         = totSim,
    maxIter        = maxIter,
    functName      = functName,
    mc.cores       = mc
  )
  
  out <- do.call(rbind, res_list)
  
  # optional: order to match your display (scenario then IC first then by mean then lambda)
  # out <- out[order(out$scenario, !out$is_ic, out$mean, out$lambda), ]
  # rownames(out) <- NULL
  out <- out[order(out$scenario, out$col_idx, out$lambda), ]
  rownames(out) <- NULL
  out
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
#### Core functions ######
DoEWMA_XieEtAl_OC <- function(maxIter, h, lambda, ICparams, OCparams)   {
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
    xVec = SampleGumbel(OCparams) 
    xVec = matrix(xVec, nrow = 2, ncol = 1)
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
  return(c(RunsCounter, TimeCounter))
}


############# Main #############
setwd("path/to/folder/EWMA_Xie/") # Set to your path to this folder
### Tweaking params ####
totSim = 1e4
ncpus = 24
maxIter = 8000
lambdaVals = c(.01, .05, .1, .2)


#### Distribution Parameters ####

theta1_vec_Gumb = c(5,7.5,10,7.5,10,5,5,2.5,1, 5,7.5,10,7.5,10,5,5,2.5,1,  5,7.5,10,7.5,10,5,5,3.5,2.5, 5,7.5,10,7.5,10,5,5,3.5,2.5)
theta2_vec_Gumb = c(5,5,5,7.5,10,2.5,1,2.5,1, 5,5,5,7.5,10,2.5,1,2.5,1,   15,15,15,22.5,30,10.5,7.5,7.5,7.5, 15,15,15,22.5,30,10.5,7.5,7.5,7.5)
delta_vec = c(rep(1,9), rep(.5, 9), rep(1,9), rep(.5, 9))
scenario = c(rep(1,9), rep(2,9), rep(3,9), rep(4,9))
lmat_Gumb = matrix( c(theta1_vec_Gumb, theta2_vec_Gumb, delta_vec, scenario), nrow = 4, byrow = T)

Gumbel_tarARL0 = c(53.25, 62, 24.55, 26)
.ic_cols_Gumb <- c(1, 10, 19, 28)

clims = read.table("dataIC/Clims_EWMA.txt", header = T)
get_h_for_pair(h_table = clims, Tar_ARL0 = 24.55, lambda = .1, theta1 = 5, theta2 = 15, delta = 1, functName = "EWMA_XieEtAl")
get_h_for_pair(h_table = clims, Tar_ARL0 = 62, lambda = .1, theta1 = lmat_Gumb[1, 10], theta2 = lmat_Gumb[2, 10], delta = lmat_Gumb[3, 10], functName = "EWMA_XieEtAl")

oc_tab <- run_EWMA_OC_all(
  lmat_Gumb      = lmat_Gumb,
  clims          = clims,
  Gumbel_tarARL0 = Gumbel_tarARL0,
  lambdaVals     = lambdaVals,
  totSim         = totSim,
  maxIter        = maxIter,
  ncpus          = ncpus
)

oc_tabs <- split(oc_tab, oc_tab$lambda)
for (nm in names(oc_tabs)) {
  cat("\n", "lambda is:", nm, "\n", sep = "")
  tab <- oc_tabs[[nm]]
  tab$lambda <- NULL      # remove lambda column
  rownames(tab) <- NULL   # remove row names (1,6,11,...)
  tab$col_idx <- NULL
  write.table(tab, file = sprintf("OC_Rls/Gumb_lambda_%s.txt", nm), quote = F, row.names = F, col.names = T, sep = "\t") 
}