## This script runs the Shewhart Chart put forth in Zwetsloot et al's 2023 paper titled: A real-time monitoring approach for bivariate event data

rm(list = ls(all.names = TRUE))
library(data.table)
library(parallel)

######## Helper functions ########
make_dist_obj <- function(name, GetExpVal, SamplePdf, Fx_v, Fy_xv, lmat, dist_arls) {
  stopifnot(is.character(name), is.function(SamplePdf), is.function(Fx_v), is.function(Fy_xv))
  stopifnot(is.matrix(lmat), nrow(lmat) >= 4) # you said last row is scenario
  list(
    name = name,
    GetExpVal = GetExpVal,
    SamplePdf = SamplePdf,
    Fx_v = Fx_v,
    Fy_xv = Fy_xv,
    lmat = lmat,
    dist_arls = dist_arls
  )
}

# Run one column of lmat
RunShewhartCol <- function(c, lmat, dist_obj, totSim, maxIter) {
  
  nr = nrow(lmat)
  scen = lmat[nr, c]
  
  # IC columns are fixed
  ic_col = c(1, 10, 19, 28)[scen]
  
  ICparams = lmat[1:(nr-1), ic_col]
  OCparams = lmat[1:(nr-1), c]
  
  alpha = 1 / dist_obj$dist_arls[scen]
  
  out = replicate(
    totSim,
    Shewart_BivPaper(
      maxIter  = maxIter,
      alpha    = alpha,
      ICparams = ICparams,
      OCparams = OCparams,
      SamplePdf = dist_obj$SamplePdf,
      Fx_v      = dist_obj$Fx_v,
      Fy_xv     = dist_obj$Fy_xv
    )
  )
  
  # out is 2 x totSim
  mean_vec = dist_obj$GetExpVal(lmat[1:(nr-1), c])
  
  data.frame(
    scenario = scen,
    is_ic    = (c %in% c(1, 10, 19, 28)),
    mean     = mean_vec,
    ARL      = round(mean(out[1, ]),1),
    sd_rl    = round(sd(out[1, ]),1),
    ATS      = round(mean(out[2, ]),1),
    sd_TS    = round(sd(out[2, ]),1),
    stringsAsFactors = FALSE
  )
}
# Run all columns for one distribution
RunShewhartDist <- function(dist_obj, totSim, ncpus, maxIter) {
  
  lmat = dist_obj$lmat
  nc = ncol(lmat)
  
  res_list = mclapply(
    X = 1:nc,
    FUN = function(c) {
      cat("Running", dist_obj$name, "column", c, "\n")
      RunShewhartCol(
        c = c,
        lmat = lmat,
        dist_obj = dist_obj,
        totSim = totSim,
        maxIter = maxIter
      )
    },
    mc.cores = ncpus
  )
  
  out = do.call(rbind, res_list)
  rownames(out) = NULL
  out
}

################ Distrbution Related Functions ############

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

## Get the expectations of X1 and X2 (NOT X,Y)
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

########## Core Function ##############

Shewart_BivPaper <- function(maxIter, alpha, ICparams, OCparams, SamplePdf, Fx_v, Fy_xv)  {
  RunsCounter = 0
  TimeCounter = 0
  
  # Note** alpha = 1/targetARLO NOT 2/targetARLO
  for (i in 1:maxIter)  {
    xvals = SamplePdf(OCparams)
    parsedVals = ParseXvec(xvals)
    
    TimeCounter = TimeCounter + parsedVals[1]
    RunsCounter = RunsCounter + 1
    if (Fx_v(x = parsedVals[1], ICparams) < (alpha/2) | Fx_v(x = parsedVals[1], ICparams) > (1 - (alpha/2)) )  {
      break
    }
    if (parsedVals[3] == 2) {
      next
    }
    
    TimeCounter = TimeCounter + parsedVals[2] - parsedVals[1]
    RunsCounter = RunsCounter + 1
    
    if (Fy_xv(x = parsedVals[1], y = parsedVals[2], v = parsedVals[3], ICparams) < (alpha/2) | Fy_xv(x = parsedVals[1], y = parsedVals[2], v = parsedVals[3], ICparams) > (1 - (alpha/2)) )  {
      break
    }
  }
  
  return( c(RunsCounter, TimeCounter) )
}


######### Main #############

setwd("path/to/dir/RunShewhartChart/RunShewhartChart") # Set path to the folder RunShewhartChart on your respective directory
### Tweaking params ####
totSim = 1e4
ncpus = 24
maxIter = 8000

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

MOBE_tarARL0 = c(53.25, 53.25, 24.55, 24.55)
MOBW_tarARL0 = c(62, 60, 26, 26)
Gumbel_tarARL0 = c(53.25, 62, 24.55, 26)

dist_MOBE <- make_dist_obj(name = "MOBE", GetExpVal = GetMOBE_Exp, SamplePdf = SampleMOBE, Fx_v = Fx_v_MOBE, Fy_xv = Fy_xv_MOBE, lmat = lmat_MOBE, dist_arls = MOBE_tarARL0)
dist_MOBW <- make_dist_obj(name = "MOBW", GetExpVal = GetMOBW_Exp, SamplePdf = SampleMOBW, Fx_v = Fx_v_MOBW, Fy_xv = Fy_xv_MOBW, lmat = lmat_MOBW, dist_arls = MOBW_tarARL0)
dist_Gumbel <- make_dist_obj(name = "Gumbel", GetExpVal = GetGumb_Exp, SamplePdf = SampleGumbel, Fx_v = Fx_v_Gumb, Fy_xv = Fy_xv_Gumb, lmat = lmat_Gumb, dist_arls = Gumbel_tarARL0)

############## Run Simulations #########
RunShewhartCol(1, dist_Gumbel$lmat, dist_Gumbel, totSim = 10, maxIter = 8000)

MOBE_Shewhart_Tab = RunShewhartDist(
  dist_obj = dist_MOBE,
  totSim = totSim,
  ncpus = ncpus,
  maxIter = maxIter
)

MOBW_Shewhart_Tab = RunShewhartDist(
  dist_obj = dist_MOBW,
  totSim = totSim,
  ncpus = ncpus,
  maxIter = maxIter
)

Gumbel_Shewhart_Tab = RunShewhartDist(
  dist_obj = dist_Gumbel,
  totSim = totSim,
  ncpus = ncpus,
  maxIter = maxIter
)

############ Make Tables #############

write.table(MOBE_Shewhart_Tab,
            file = "OC_Rls/MOBE_Shewhart_Tab.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(MOBW_Shewhart_Tab,
            file = "OC_Rls/MOBW_Shewhart_Tab.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(Gumbel_Shewhart_Tab,
            file = "OC_Rls/Gumbel_Shewhart_Tab.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

