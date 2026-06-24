rm(list = ls(all.names = TRUE))
library(data.table)

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

Shewart_BivPaper <- function(maxIter, alpha, data, ICparams, Fx_v, Fy_xv)  {
  # Alpha is just 1/targetARLO NOT 2/targetARLO
  maxIter = min(maxIter, nrow(data))
  
  RunsCounter = 0
  TimeCounter = 0
  
  probs = c()
  for (i in 1:maxIter)  {
    xvals = data[i,]
    parsedVals = ParseXvec(xvals)
    TimeCounter = TimeCounter + parsedVals[1]
    RunsCounter = RunsCounter + 1
    cat("Fx_v(X|V): ", Fx_v(x = parsedVals[1], ICparams), "\n")
    probs[RunsCounter] = Fx_v(x = parsedVals[1], ICparams)
    if (Fx_v(x = parsedVals[1], ICparams) < (alpha/2) | Fx_v(x = parsedVals[1], ICparams) > (1 - (alpha/2)) )  {
      break
    }
    if (parsedVals[3] == 2) {
      next
    }
    
    # clims2 = MOBE_Clims_2(parsedVals[1], parsedVals[3], IC1, IC2, IC3, ICeta, alpha)
    TimeCounter = TimeCounter + parsedVals[2] - parsedVals[1]
    RunsCounter = RunsCounter + 1
    probs[RunsCounter] = Fy_xv(x = parsedVals[1], y = parsedVals[2], v = parsedVals[3], ICparams)
    cat("Fy_xv(Y|X,V): ", Fy_xv(x = parsedVals[1], y = parsedVals[2], v = parsedVals[3], ICparams), "\n")
    if (Fy_xv(x = parsedVals[1], y = parsedVals[2], v = parsedVals[3], ICparams) < (alpha/2) | Fy_xv(x = parsedVals[1], y = parsedVals[2], v = parsedVals[3], ICparams) > (1 - (alpha/2)) )  {
      break
    }
  }
  
  totalTime = 0
  for (i in 1:nrow(data)) {
    xvals = data[i,]
    parsedVals = ParseXvec(xvals)
    totalTime = totalTime + parsedVals[2]
  }
  if (TimeCounter == totalTime) {
    cat("\nNo alarm triggered***\n\n")
    cat("The total observations, time and vectors are: ", RunsCounter, TimeCounter, i, "\n\n")
  } else  {
    cat("The ARL, ATS, total Vectors are: ", RunsCounter, TimeCounter, i, "\n")
  }
  return( c(RunsCounter, TimeCounter, i) )
}
Biv8_CUSUM_TimVar_data <- function(maxIter, h, data, priorParams = c(22.05, 21, 9.5, 10), rhoVec = c(1.05, .95), ICparams, Fx_v, Fy_xv, cdfsList, dist_map)   {
  
  # Modified to handle data from a data set opposed to generating them through simulation
  
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
  lambda_vec2 = l2_defaultVals
  
  lambda_vec3 = l3_defaultVals
  
  lamPosInd1 = 1:4
  lamNegInd1 = 5:8
  lamPosInd2 = c(1,2,5,6)
  lamNegInd2 = setdiff(1:8, lamPosInd2)
  lamPosInd3 = seq(1,7,2)
  lamNegInd3 = setdiff(1:8, lamPosInd3)
  
  b = 0 
  
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
  
  # ----- START REFACTORED ARL SCRIPT LOOP -----
  
  # --- Initializations for a SINGLE ARL RUN ---
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
  # (You used parsedVals[1] for X_ord1, parsedVals[2] for X_ord2)
  X_ord1_current_pair = NA
  Y_ord2_current_pair = NA # Only strictly needed if V != 2
  
  storedVals = matrix(ncol = 4)
  Cmat = matrix(ncol = 8)
  
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
      
      xvals = data[pair_idx, ]
      parsedVals = ParseXvec(xvals)
      cat("The parsedVals are: ", parsedVals, "\n")
      
      X_ord1_current_pair = parsedVals[1]
      Y_ord2_current_pair = parsedVals[2] # Store even if V=2, for consistency
      V_of_current_pair = parsedVals[3]
      
      U1 = Fx_v(X_ord1_current_pair, ICparams) # Transform using IC params
      Z_current_rt = qexp(U1, rate = 1)
      L_type_current_rt = 1
      time_increment_current_rt = X_ord1_current_pair
      if (rt_iter_idx == 1) {
        storedVals[rt_iter_idx, ] = c(X_ord1_current_pair, U1, Z_current_rt, L_type_current_rt)
      } else  {
        storedVals = rbind(storedVals, c(X_ord1_current_pair, U1, Z_current_rt, L_type_current_rt))
      }
      
      
    } else { # Processing the second observation of the pair
      # This means V_of_current_pair must have been 0 or 1
      U2 = Fy_xv(x = X_ord1_current_pair, y = Y_ord2_current_pair, v = V_of_current_pair, ICparams) # IC params
      Z_current_rt = qexp(U2, rate = 1)
      L_type_current_rt = if (V_of_current_pair == 0) 2 else 3
      time_increment_current_rt = Y_ord2_current_pair - X_ord1_current_pair
      
      storedVals = rbind(storedVals, c(Y_ord2_current_pair, U2, Z_current_rt, L_type_current_rt))
    }
    
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
    if (RunsCounter <= 400) {
      condCdfs <- cdfsList[[ dist_map[RunsCounter] ]]
    } else  {
      condCdfs = cdfsList[["stable"]]
    }
    # Compute p-values
    pvals = 1 - vapply(seq_along(condCdfs), function(j) condCdfs[[j]](Cvec[j]), numeric(1))
    cusumStndzd = -log(pvals)
    cat("RunsCounter and TimeCounter are: ", RunsCounter, TimeCounter, "\n")
    cat("Cvec: ", Cvec, "\n")
    cat("cusumStdzd: ", cusumStndzd, "\n\n")
    
    if (rt_iter_idx == 1) {
      Cmat[rt_iter_idx, ] = cusumStndzd
    } else  {
      Cmat = rbind(Cmat, cusumStndzd)
    }
    
    if (max(cusumStndzd) > h | (pair_idx >= nrow(data) & is_first_obs_of_pair)) {
      break
    }
    
  } # End of rt_iter_idx loop (processing one R_t at a time)
  
  
  ARL_val = RunsCounter
  ATS_val = TimeCounter
  
  totalTime = 0
  for (i in 1:nrow(data)) {
    xvals = data[i,]
    parsedVals = ParseXvec(xvals)
    totalTime = totalTime + parsedVals[2]
  }
  if (TimeCounter == totalTime) {
    cat("\nNo alarm triggered***\n\n")
    cat("The total observations, time and vectors are: ", RunsCounter, TimeCounter, pair_idx, "\n\n")
  } else  {
    cat("The ARL, ATS, total Vectors are: ", RunsCounter, TimeCounter, pair_idx, "\n")
  }
  
  # cat("ARL, ATS, VectorInd: ", c(ARL_val, ATS_val, pair_idx), "\n")
  colnames(storedVals) = c("ObsVal", "unifVal", "TransfVal", "Ltype")
  # return(c(ARL_val, ATS_val, pair_idx))
  return(list(storedVals = storedVals, Cmat = Cmat))
  
}

####### Main ######
setwd("/path/to/folder/RealDataAnalysis/") # Set to your working directory where the folder RealDataAnalysis exists
L = readRDS(file = "DataAndParamsList.RDS")

ICdata = L$ICdata
OCdata = L$OCdata
ICparams = L$ICparams


Shewart_BivPaper(200, alpha = 1/50, data = OCdata, ICparams = ICparams, Fx_v = Fx_v_MOBW, Fy_xv = Fy_xv_MOBW)

###### For My CUSUM Statistic #########

############ Load Stable Distribution ############

stable_obj <- readRDS(file.path("dataIC/CUSUM_dist/Biv8_stable_pr_22.1_21_9.5_10_rho_1.05_0.95_p1_0.50.rds"))
tab = stable_obj$stable_dist
cdfs_stable <- lapply(1:8, function(i) {
  Ftemp <- ecdf(tab[, i])
  function(x) (Ftemp(x) - Ftemp(0)) / (1 - Ftemp(0))
})
####### Prepare the list of Time Varying CDFs ########
# 1. Define the time points I will load
snapshot_times <- c(1:20, 25, seq(30, 100, by = 10), seq(150, 200, by = 50))
snapshot_times
# 2. stored_cdfs will be indexed by character name
time_var_dir = file.path("dataIC", "TimeVarDist")
TimeVar_cdfs <- list()
for (t in snapshot_times) {
  # Construct path
  fpath = file.path(time_var_dir, sprintf("TimeVar_BivAdpExp_p1_0.45_t_%s.txt", t))
  
  tab <- as.matrix(fread(fpath))
  cdfs = sapply(1:8, function(c) ecdf(tab[,c]), simplify = F)  # Technically for the mean only CUSUM you only need the first two components of this. The second two elements will be the CDFs for delta
  condCdfs = lapply(1:8, function(i) {
    Ftemp = cdfs[[i]]
    function(x) (Ftemp(x) - Ftemp(0))/(1 - Ftemp(0))
  })
  TimeVar_cdfs[[as.character(t)]] <- condCdfs
}
warnings()
TimeVar_cdfs[["stable"]] = cdfs_stable

dist_map <- character(200) # For t=1 to 100

for (t in 1:200) {
  # findInterval returns the index of the closest lower bound in snapshot_times
  # e.g., if t=22, it finds index of '20'
  idx <- findInterval(t, snapshot_times)
  
  # If t is less than the first snapshot (shouldn't happen with 1:10), handle safely
  if (idx == 0) idx <- 1 
  
  # Map t to the specific snapshot name
  dist_map[t] <- as.character(snapshot_times[idx])
}
# dist_map
############# Run the CUSUM charting Statistic ###########
h = 2.933 # Control limit that gives you and ARL0 of 50
output = Biv8_CUSUM_TimVar_data(maxIter = 200, h = 2.933, data = OCdata, ICparams = ICparams, Fx_v = Fx_v_MOBW, Fy_xv = Fy_xv_MOBW, cdfsList = TimeVar_cdfs, dist_map = dist_map)


