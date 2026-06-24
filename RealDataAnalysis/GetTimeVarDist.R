# Very simple code to obtain the time varying distribution of the CUSUM statistic

rm(list = ls(all.names = TRUE))
require(parallel)
library(data.table)

ExpCUSUM_Biv8_timeVary <- function(maxIter, p1, j, priorParams= c(22.05, 21, 9.5, 10), rhoVec = c(1.05, .95), totSim)   {
  
  startBool = if (j == 1) T else F
  endBool = if (j == totSim) T else F
  
  Cvec = Svec1 = Svec2 = Svec3 = Tvec1 = Tvec2 = Tvec3 = rep(0,8)
  
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
  Z1 = 0
  
  # temp = c(incVal, decVal, incVal, decVal)
  lambda_vec2 = l2_defaultVals
  Z2 = 0
  
  lambda_vec3 = l3_defaultVals
  Z3 = 0
  
  lamPosInd1 = 1:4
  lamNegInd1 = 5:8
  lamPosInd2 = c(1,2,5,6)
  lamNegInd2 = setdiff(1:8, lamPosInd2)
  lamPosInd3 = seq(1,7,2)
  lamNegInd3 = setdiff(1:8, lamPosInd3)
  
  b = 0 # Doesn't matter whether the starting value is 0 or 1 because either way the first update will return Svec = = Tvec = 0vec for first iteration since Cvec = 0vec
  
  GetEstimates <- function(Cvec, Svec, Tvec, lamPosInd, lamNegInd, l_defaultVals, Z)  {
    posInd = (Cvec > 0)
    lambda_vec = rep(NA, 8)
    
    Svec[posInd] = Svec[posInd] + Z
    Svec[!posInd] = 0
    
    Tvec[posInd] = Tvec[posInd] + 1
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
    
    return(list(Svec, Tvec, lambda_vec))
  }
  ResetEstimates <- function(Cvec, Svec, Tvec, lamPosInd, lamNegInd, l_defaultVals)  {
    posInd = (Cvec > 0)
    lambda_vec = rep(NA, 8)
    
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
    
    return(list(Svec, Tvec, lambda_vec))
  }
  
  for (i in 1:maxIter)  {
    ### Start of operations done with previous values ###
    if (i%%2 == 0)  { # These are offset because we update using previous value to preserve martingale property
      L_res1 = GetEstimates(Cvec, Svec1, Tvec1, lamPosInd1, lamNegInd1, l1_defaultVals, Z1)
      L_res2 = ResetEstimates(Cvec, Svec2, Tvec2, lamPosInd2, lamNegInd2, l2_defaultVals)
      L_res3 = ResetEstimates(Cvec, Svec3, Tvec3, lamPosInd3, lamNegInd3, l3_defaultVals)
      Svec1 = L_res1[[1]]; Tvec1 = L_res1[[2]]; lambda_vec1 = L_res1[[3]]
      Svec2 = L_res2[[1]]; Tvec2 = L_res2[[2]]; lambda_vec2 = L_res2[[3]]
      Svec3 = L_res3[[1]]; Tvec3 = L_res3[[2]]; lambda_vec3 = L_res3[[3]]
    } else if (i%%2 == 1) {
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
    
    if (i%%2 == 1 & i != maxIter)  { # By not updating Cvec here we don't need to output the lambda_vecs and Z1, Z2, Z3
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
    if (startBool)  {
      cat(Cvec, sep = "\t", end = "\n", file = sprintf("dataIC/TimVarDist/TimeVar_BivAdpExp_p1_%s_t_%s.txt", p1, i), append = F)
    } else if (endBool) {
      cat(Cvec, sep = "\t", end = "", file = sprintf("dataIC/TimVarDist/TimeVar_BivAdpExp_p1_%s_t_%s.txt", p1, i), append = T)
    } else  {
      cat(Cvec, sep = "\t", end = "\n", file = sprintf("dataIC/TimVarDist/TimeVar_BivAdpExp_p1_%s_t_%s.txt", p1, i), append = T)
    }
  }
  # return(c(Cvec, Svec1, Svec2, Svec3, Tvec1, Tvec2, Tvec3))
}

############## Main ###########
totSim = 1e4
p1 = .45 # p1 is P(V = 1 | V != 2)

setwd("/path/to/folder/RealDataAnalysis/") # Set to your working directory aka only change path/to/folder to wherever the folder RealDataAnalysis is

m = sapply(1:totSim, function(j) ExpCUSUM_Biv8_timeVary(maxIter = 400, p1 = p1, j = j, totSim = totSim) )