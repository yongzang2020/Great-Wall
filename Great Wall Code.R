#######################################################################################################################
## UFO: Great Wall: A Generalized Phase I-II Dose-Optimization Design for Drug Combination Trials Maximizing Survival Benefit
##          by Yan Han, Yingjie Qiu and Yong Zang
##
## This file contains two parts for simulation studies and trial implementation of Great Wall designs to find the ODC
## Part 1: simulation studies
## get.Great_Wall.oc() is a function used to generate operating characteristics for the proposed designs
## Part 2: trial implementation
## Decision_stageI() is used to get dose escalation decision for incoming cohort of the stage I of Great Wall design 
## generate_path_S1_values() is used to get dose finding path at stage I  
## decision_end_stageI() is used to get decision by the end of stage I
## decision_end_stageII() is used to get decision by the end of stage II 
## get_surv_est() is used to get survival estimates by the end of stage III  
## decision_end_stageIII() is used to get decision by the end of stage III  


########################################################################################################################

##########################################################################################################
## Function to generate operating characteristics of the UFO design 
## To use this function, library "Iso" and "survival" should be installed in R. 
##
##  Arguments:
## maxn1: maximum sample size for stage I
## maxn2: maximum sample size for stage II
## maxn3: maximum sample size for stage III
## phiT: upper limit of toxicity rate 
## rhoT: rhoT times phiT be the toxicity rate that is deemed overly toxic 
## true_matrix_tox : True toxicity rate of drug combination 
## true_matrix_eff : True efficacy rate of drug combination 
## true_matrix_surv : True PFS probability of drug combination 
## alpha: percentage of disease progression outcome occured in latter half of the assessment window
## maxt: survival outcome follow up time
## phiE: lower limit of efficacy rate 
## tau_E: cutoff probability for futile combination
## alpha_E: hyper parameter for beta prior of efficacy 
## beta_E: hyper parameter for beta prior of efficacy 
## cohortsize1: sample size for each cohort at stage I
## sim_trial: the number of simulated trial
## rho: assocaition parameter between efficacy and toxicity at each combination
## gamma:  flexibility of candidate set selection at end of stage II
## phi_s: lower limite for survival measurement

## Utility_score: utility score for each possible outcome [(T,E): (0,0), (0,1), (1,0), (1,1)]


#########################################################################################################


get.Great_Wall.oc<-function(maxn1,maxn2,maxn3, phiT,rhoT,true_matrix_tox, 
                     true_matrix_eff, true_matrix_surv,
                     surv_dist,alpha,maxt,
                     phiE, tau_E,alpha_E,beta_E,
                     cohortsize1,sim_trial,rho,gamma, phi_s, Utility_score)
{
  simulation_N <- array(0, dim = c(nrow(true_matrix_tox), 
                                   ncol(true_matrix_tox), 
                                   sim_trial))
  simulation_Tox <- array(0, dim = c(nrow(true_matrix_tox), 
                                     ncol(true_matrix_tox), 
                                     sim_trial))
  simulation_Eff <- array(0, dim = c(nrow(true_matrix_tox), 
                                     ncol(true_matrix_tox), 
                                     sim_trial))
  simulation_Uti <- array(NA, dim = c(nrow(true_matrix_tox), 
                                      ncol(true_matrix_tox), 
                                      sim_trial))
  #  simulation_Surv <- array(NA, dim = c(nrow(true_matrix_surv), 
  #                                   ncol(true_matrix_surv), 
  #                                   sim_trial))
  d_opt = data.frame(row= rep(0,sim_trial) ,col=rep(0,sim_trial)) 
  
  max_cor_size_s1 = maxn1/cohortsize1
  
  ################ load function ############
  library(survival)
  ####### Gumbel Copula ######
  Gumbel=function(ttox,teff,c){
    ## ttox: marginal toxicity probability
    ## teff: marginal efficacy probability
    ## c: association parameter between tox and eff, c>0 indicates a positive correlation
    ndose=length(ttox) ## dose level
    out=matrix(rep(0,4*ndose),nrow=4) ## joint tox-eff probability matrix; row is dose level, column=(tox=0,eff=0; tox=0,eff=1; tox=1,eff=0; tox=1,eff=1)
    for (j in 1: ndose){
      out[1,j]=(1-ttox[j])*(1-teff[j])+ttox[j]*teff[j]*(1-ttox[j])*(1-teff[j])*( exp(c)-1 )/(exp(c)+1)
      out[2,j]=(1-ttox[j])*(teff[j])-ttox[j]*teff[j]*(1-ttox[j])*(1-teff[j])*( exp(c)-1 )/(exp(c)+1)
      out[3,j]=(ttox[j])*(1-teff[j])-ttox[j]*teff[j]*(1-ttox[j])*(1-teff[j])*( exp(c)-1 )/(exp(c)+1)
      out[4,j]=(ttox[j])*(teff[j])+ttox[j]*teff[j]*(1-ttox[j])*(1-teff[j])*( exp(c)-1 )/(exp(c)+1)
    }
    return(out)
  }
  
  
  Gumbel2=function(tox_matrix, eff_matrix, c){
    ## tox_matrix: matrix of toxicity probabilities
    ## eff_matrix: matrix of efficacy probabilities
    ## c: association parameter between tox and eff, c>0 indicates a positive correlation
    
    nrow_mat = nrow(tox_matrix)
    ncol_mat = ncol(tox_matrix)
    
    ## Initialize 3D array to store joint probabilities
    joint_probs_array = array(0, dim = c( nrow_mat, ncol_mat,4))
    
    for (i in 1:nrow_mat){
      for (j in 1:ncol_mat){
        ttox = tox_matrix[i,j]
        teff = eff_matrix[i,j]
        
        joint_probs_array[i, j,1] = (1-ttox)*(1-teff) + ttox*teff*(1-ttox)*(1-teff)*(exp(c)-1)/(exp(c)+1)
        joint_probs_array[i, j,2] = (1-ttox)*teff - ttox*teff*(1-ttox)*(1-teff)*(exp(c)-1)/(exp(c)+1)
        joint_probs_array[i, j,3] = ttox*(1-teff) - ttox*teff*(1-ttox)*(1-teff)*(exp(c)-1)/(exp(c)+1)
        joint_probs_array[i, j,4] = ttox*teff + ttox*teff*(1-ttox)*(1-teff)*(exp(c)-1)/(exp(c)+1)
      }
    }
    return(joint_probs_array)
  }
  
  ############ get outcome for individual outcome ##########
  generate_multinomial_outcome2 = function(tox_d, eff_d, c, num_simulations){
    
    # Generate joint probabilities using the Gumbel2 function
    joint_probs = Gumbel(tox_d, eff_d, c)
    
    outcome=rmultinom(1,num_simulations,joint_probs)
    
    return(outcome)
  }
  
  
  ###################get phi ###########
  get_phi_boundary<-function(phi_T,rho){
    return(log((1-phi_T)/(1-rho*phi_T))/log(rho*((1-phi_T)/(1-rho*phi_T))))
  }
  
  ################## get path ########
  generate_path_S1_values <- function(mat) {
    nrow_mat <- nrow(mat)
    ncol_mat <- ncol(mat)
    
    # Start from bottom left and go up the first column
    path <- cbind(seq(nrow_mat, 1, by = -1), rep(1, nrow_mat))
    
    # Then move horizontally to the top right from the first row
    if (ncol_mat > 1) {
      path <- rbind(path, cbind(rep(1, ncol_mat - 1), seq(2, ncol_mat)))
    }
    
    # Extract values from matrix using the path indices
    values <- mat[cbind(path[,1], path[,2])]
    values_matrix <- matrix(values, nrow = 1)
    list(values = values_matrix, positions = path)
  }
  
  
  ####### Function to get the indices of S1 path ##############
  get_S1_indices <- function(nrow_mat, ncol_mat) {
    # Start from bottom left and go up the first column
    path <- cbind(seq(nrow_mat, 1, by = -1), rep(1, nrow_mat))
    
    # Then move horizontally to the top right from the first row
    if (ncol_mat > 1) {
      path <- rbind(path, cbind(rep(1, ncol_mat - 1), seq(2, ncol_mat)))
    }
    
    return(path)
  }
  
  ######### Function to exclude toxic doses and return their indices #########
  exclude_toxic_doses_indices <- function(nrow_mat, ncol_mat, toxic_position) {
    # Exclude doses above the toxic position in the same column and columns to the right
    toxic_rows <- 1:toxic_position[1]
    toxic_cols <- toxic_position[2]:ncol_mat
    
    toxic_indices <- expand.grid(toxic_rows, toxic_cols)
    return(as.matrix(toxic_indices))
  }
  
  ######## Function to remove the union of doses from original matrix ########
  remove_doses <- function(mat, remove_indices) {
    for(i in 1:nrow(remove_indices)) {
      mat[remove_indices[i, 1], remove_indices[i, 2]] <- NA
    }
    return(mat)
  }
  ########################################################################
  compress_matrix <- function(mat) {
    # Identify rows and columns without any NA values
    valid_rows <- apply(mat, 1, function(row) !all(is.na(row)))
    valid_cols <- apply(mat, 2, function(col) !all(is.na(col)))
    
    # Subset the matrix to only valid rows and columns
    compressed_mat <- mat[valid_rows, valid_cols,drop=FALSE]
    
    return(compressed_mat)
  }
  ######################################################
  assign_submatrix <- function(original_matrix, sub_matrix) {
    # Determine starting position for the submatrix
    start_row <- nrow(original_matrix) - nrow(sub_matrix) + 1
    start_col <- ncol(original_matrix) - ncol(sub_matrix) + 1
    
    # Assign values of submatrix to original matrix
    original_matrix[start_row:nrow(original_matrix), start_col:ncol(original_matrix)] <- sub_matrix
    
    return(original_matrix)
  }
  #####################################################
  
  assign_arrary <- function(original_array, sub_array) {
    # Determine starting position for the submatrix
    #start_row <- nrow(original_matrix) - nrow(sub_matrix) + 1
    #start_col <- ncol(original_matrix) - ncol(sub_matrix) + 1
    
    start_row <- dim(original_array)[1] - dim(sub_array)[1] + 1
    start_col <- dim(original_array)[2] - dim(sub_array)[2] + 1
    
    # Assign values of submatrix to original matrix
    original_array[start_row:nrow(original_array), start_col:ncol(original_array),] <- sub_array
    
    return(original_array)
  }
  
  ####################################################
  divide_equally <- function(total, n_parts) {
    
    # Check if total can be divided by n_parts without a remainder
    if (total %% n_parts == 0) {
      return(rep(total / n_parts, n_parts))
    }
    
    # If not divisible evenly, then use the previous approach
    base <- total %/% n_parts      # Integer division
    remainder <- total %% n_parts  # Modulo operation
    
    # Create a vector with each part as the base value
    parts <- rep(base, n_parts)
    
    # Distribute the remainder to the first few parts
    parts[1:remainder] <- parts[1:remainder] + 1
    
    return(parts)
  }
  
  ###################################################
  get_stage2_outcome<-function(samplesize, true_tox, 
                               true_eff,location){
    
    YT_outcome_matrix <- matrix(0, nrow(true_tox), ncol(true_tox))
    YE_outcome_matrix <- matrix(0, nrow(true_eff), ncol(true_eff))
    N_outcome_matrix <- matrix(0, nrow(true_tox), ncol(true_tox))
    
    joint_TE_array_S2<-array(0, dim = c(nrow(true_tox), 
                                        ncol(true_tox), 4))
    
    # For each intersection location, simulate outcomes
    for (i in 1:nrow(location)) {
      r <- location[i, "row"]
      c <- location[i, "col"]
      N_outcome_matrix[r, c] <- samplesize[i]
      prob_T <- true_tox[r, c]
      prob_E <- true_eff[r, c]
      
      S2_outcome = generate_multinomial_outcome2(prob_T,prob_E,rho,samplesize[i])
      
      YT_outcome_matrix[r, c] <- sum(S2_outcome[3,],S2_outcome[4,])
      YE_outcome_matrix[r, c] <- sum(S2_outcome[2,],S2_outcome[4,])
      joint_TE_array_S2[r,c,] <- S2_outcome
      
    }
    return(list(N2_matrix=N_outcome_matrix,
                YT2_matirx=YT_outcome_matrix,
                YE2_matrix=YE_outcome_matrix,
                joint_TE_arrary2 =joint_TE_array_S2))
  }
  ####################################################
  library(Iso)
  model_T <- function(nti,yti) {
    pt = yti/nti
    pt[is.na(pt)] = 0
    pt_flipped <- pt[nrow(pt):1, ]
    nti_flipped<- nti[nrow(nti):1, ]
    rate_flipped<-biviso(y=pt_flipped,w=nti_flipped,fatal = F,warn = F)
    T_est=rate_flipped[nrow(rate_flipped):1,]
    return(T_est)
  }
  
  
  get_est_dose<-function(input_arrary,N_matrix){
    rate_array <- array(0, dim = dim(input_arrary))
    
    # Loop through the third dimension and add corresponding matrices
    for (k in 1:dim(input_arrary)[3]) {
      rate_array[,,k] <- input_arrary[,,k]/N_matrix
      rate_array[,,k][is.na(rate_array[,,k])] = 0
    }
    return(rate_array)
  }
  
  #########
  ## (T,E): (0,0), (0,1), (1,0), (1,1)
  #Uti_score<-array(c(40,100,0,60),c(1,1,4))
  Uti_score = array(Utility_score,c(1,1,4))
   
  Uti_EST<-function(Utility,prop_est){
    Uti_est <- matrix(0, nrow = dim(prop_est)[1],ncol = dim(prop_est)[2])
    result_array <- array(0, dim(prop_est))
    
    for (k in 1:dim(prop_est)[3]) {
      result_array[,,k] <- prop_est[,,k] * Utility[,,k]
    }
    Uti_est = result_array[,,1]+result_array[,,2]+
      result_array[,,3]+result_array[,,4]
    
    return(Uti_est)
  }
  ######################
  True_Uti<-Uti_EST(Uti_score,Gumbel2(true_matrix_tox,true_matrix_eff,rho))/100
  #################
  find_high_values_in_intersection <- function(matrix, intersection,gamma) {
    extract_values <- function(matrix, intersection){
      sapply(1:nrow(intersection), function(i) {
        matrix[intersection[i, "row"], intersection[i, "col"]]
      })
    }
    
    # Use the function to extract values
    extracted_values <- extract_values(matrix, intersection2)
    
    max_value <- max(extracted_values)
    
    # Identify locations in Mean_Uti_est with values >= gamma * max_value
    identify_locations <- function(matrix, threshold){
      which(matrix >= threshold, arr.ind = TRUE)
    }
    
    locations <- identify_locations(matrix, gamma * max_value)
    return(locations)
  }
  
  closest_value_and_position <- function(mat, target = 0.3) {
    # Subtract target from each element
    diff_matrix <- abs(mat - target)
    
    # Find the positions of all the minimum values
    positions <- which(diff_matrix == min(diff_matrix), arr.ind = TRUE)
    
    # Sort by column first (to get the most left) and then by row in decreasing order (to get the most bottom)
    sorted_positions <- positions[order(positions[, "col"], -positions[, "row"]), ]
    
    # Extract the most left and most bottom position
    desired_position <- sorted_positions[1, ]
    
    # Retrieve the closest value to target from the original matrix
    closest_value <- mat[desired_position["row"], desired_position["col"]]
    
    list(value = closest_value, position = desired_position)
  }
  
  ####################################################
  
  gen.tite <- function(dist = 1,
                       n,
                       pi,
                       alpha = 0.5,
                       Tobs = 1) {
    ############ subroutines ############
    weib <- function(n, pi, pihalft)
    {
      ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
      
      alpha = log(log(1 - pi) / log(1 - pihalft)) / log(2)
      
      lambda = -log(1 - pi) / (Tobs ^ alpha)
      
      t = (-log(runif(n)) / lambda) ^ (1 / alpha)
      
      return(list(t=t, b=(1/lambda)^(1/alpha),k=alpha))
      
    }
    
    llogit <- function(n, pi, pihalft)
    {
      ## solve parameters for log-logistic given pi=1-S(T) and phalft=1-S(T/2)
      #alpha = log((1 / (1 - pi) - 1) / (1 / (1 - pihalft) - 1)) / log(2)
      
      #lambda = (1 / (1 - pi) - 1) / (Tobs ^ alpha)
      
      #t = ((1 / runif(n) - 1) / lambda) ^ (1 / alpha)
      
      kk = log((pi*(1-pihalft))/(pihalft*(1-pi)))/log(2)
      bb = Tobs*((1-pi)/pi)^(1/kk)
      r = runif(n)
      t = bb * (r / (1 - r)) ^ (1 / kk)
      
      return(list(t=t, b=bb,k=kk))
      
    }
    
    pareto <- function(n, pi, pihalft){
      
      alpha = log((1-pihalft)/(1-pi))/log(2)
      beta= Tobs * (1-pi)^(1/alpha)
      r = runif(n)
      t = beta/r^(1/alpha)
      return(list(t=t, b=beta,a=alpha))
    }
    
    ############ end of subroutines ############
    
    
    tox = rep(0, n)
    
    t.tox = rep(0, n)
    
    #### Weibull
    if (dist == 1)
    {
      pihalft = alpha * pi
      # alpha*100% event in (0, 1/2T)
      t.tox = weib(n, pi, pihalft)$t
      
      tox[t.tox <= Tobs] = 1
      
      ntox.st = sum(tox)
      
      t.tox[tox == 0] = Tobs
      b = weib(n, pi, pihalft)$b
      k = weib(n, pi, pihalft)$k
      S_t = function(t){
        exp(-(t/b)^k)
      }
      RMST = integrate(S_t, lower = 0, upper = Tobs)
      Med = b*(log(2))^(1/k)
    }
    #### log-logistic
    else if (dist == 2)
    {
      pihalft = alpha * pi
      # alpha*100% event in (0, 1/2T)
      t.tox = llogit(n, pi, pihalft)$t
      
      tox[t.tox <= Tobs] = 1
      
      ntox.st = sum(tox)
      
      t.tox[tox == 0] = Tobs
      b = llogit(n, pi, pihalft)$b
      k = llogit(n, pi, pihalft)$k
      S_t = function(t){
        1 - (t^k/(b^k+t^k))
      }
      RMST = integrate(S_t, lower = 0, upper = Tobs)
      Med = b
    }
    else if (dist == 3)
    {
      pihalft = alpha * pi
      # alpha*100% event in (0, 1/2T)
      t.tox = pareto(n, pi, pihalft)$t
      
      tox[t.tox <= Tobs] = 1
      
      ntox.st = sum(tox)
      
      t.tox[tox == 0] = Tobs
      b = pareto(n, pi, pihalft)$b
      k = pareto(n, pi, pihalft)$a
      S_t <- function(t) {
        ifelse(t < b, 1, (b/t)^k)
      }
      
      RMST = ifelse(Tobs < b, Tobs, integrate(S_t, lower = 0, upper = Tobs)$value)
      Med = b*(2^(1/k))
    }
    
    return(list(
      tox = tox,
      t.tox = t.tox,
      b = b,
      k = k,
      RMST_true = RMST,
      Med_true = Med
    ))
    
  }
  
  get_stage3_outcome_TES<-function(samplesize, true_tox, 
                                   true_eff,true_surv,
                                   dist_surv,alpha_S,ass_window,
                                   location, cur_N){
    
    YT_outcome_matrix <- matrix(0, nrow(true_tox), ncol(true_tox))
    YE_outcome_matrix <- matrix(0, nrow(true_eff), ncol(true_eff))
    N_outcome_matrix <- matrix(0, nrow(true_tox), ncol(true_tox))
    
    joint_TE_array_S2<-array(0, dim = c(nrow(true_tox), 
                                        ncol(true_tox), 4))
    
    total_surv_sample<- samplesize+sapply(1:nrow(location), function(i) {
      cur_N[location$row[i], location$col[i]]
    })
    
    surv_out_arary<-lapply(total_surv_sample, function(rows) {
      matrix(0, nrow = rows, ncol = 2)})
    
    # For each intersection location, simulate outcomes
    for (i in 1:nrow(location)) {
      r <- location[i, "row"]
      c <- location[i, "col"]
      N_outcome_matrix[r, c] <- samplesize[i]
      prob_T <- true_tox[r, c]
      prob_E <- true_eff[r, c]
      
      surv_rate<- 1 - true_surv[r, c]
      
      
      surv_data<-gen.tite(dist=dist_surv,n=total_surv_sample[i],
                          pi=surv_rate,alpha = alpha_S,Tobs = ass_window)
      
      
      surv_out_arary[[i]]<-matrix(c(surv_data$t.tox,surv_data$tox),ncol = 2)
      
      
      S2_outcome = generate_multinomial_outcome2(prob_T,prob_E,rho,samplesize[i])
      
      YT_outcome_matrix[r, c] <- sum(S2_outcome[3,],S2_outcome[4,])
      YE_outcome_matrix[r, c] <- sum(S2_outcome[2,],S2_outcome[4,])
      joint_TE_array_S2[r,c,] <- S2_outcome
      
    }
    return(list(N3_matrix=N_outcome_matrix,
                YT3_matirx=YT_outcome_matrix,
                YE3_matrix=YE_outcome_matrix,
                joint_TE_arrary3 =joint_TE_array_S2,
                surv_array = surv_out_arary))
  }
  
  get_surv_est<-function(surv_array, ass_window){
    results<-rep(NA,length(surv_outcome))
    for(i in 1:length(surv_outcome)) {
      # Extract time and status
      time <- surv_outcome[[i]][, 1]
      status <- surv_outcome[[i]][, 2]
      
      # Perform the survival analysis
      fit <- survfit(Surv(time, status) ~ 1)
      
      # Extract summary at time = 90
      result <- summary(fit, times = ass_window,extend = T)$surv
      
      # Store the result
      results[i] <- result
    }
    return(results)
  }
  
  ####################################################
  ntrial<-nrow(true_matrix_tox)
  ###### get phi ######
  
  phi = get_phi_boundary(phiT,rhoT)
  
  for (trial in 1:sim_trial) {
    S1 = 1
    if(S1==1){
      ######## get S1 #########
      tox_values_S1 <- generate_path_S1_values(true_matrix_tox)$values
      
      eff_values_S1 <- generate_path_S1_values(true_matrix_eff)$values
      #########################
      
      ########total data###########
      N<-matrix(0,nrow = nrow(true_matrix_tox),ncol = ncol(true_matrix_tox))
      YTOX<-matrix(0,nrow = nrow(true_matrix_tox),ncol = ncol(true_matrix_tox))
      YEFF<-matrix(0,nrow = nrow(true_matrix_eff),ncol = ncol(true_matrix_eff))
      joint_TE_array = array(0, dim = c(nrow(true_matrix_tox), 
                                        ncol(true_matrix_tox), 4))
      ######subtrial S1##########
      d = 1
      ytox1 = rep(0,length(tox_values_S1))
      yeff1 = rep(0,length(eff_values_S1))
      n1 = rep(0,length(tox_values_S1))
      ytox_eff1 = array(0, dim = c(nrow(tox_values_S1), 
                                   ncol(tox_values_S1), 4))
      
      
      ncorhort1 = length(tox_values_S1)
      tox_dose = 0
      
      for (i_s1 in 1:ncorhort1){
        ind_outcomes<- generate_multinomial_outcome2(tox_values_S1[,d], 
                                                     eff_values_S1[,d],
                                                     rho, cohortsize1)
        ytox1[d] = ytox1[d] + sum(ind_outcomes[3,],ind_outcomes[4,])
        yeff1[d] = yeff1[d] + sum(ind_outcomes[2,],ind_outcomes[4,])
        ytox_eff1[,d,]<-ind_outcomes
        
        n1[d] = n1[d] + cohortsize1
        
        pT_hat = ytox1[d]/n1[d]
        
        if(pT_hat<phi){
          d = d + 1
        }else{
          tox_dose = d
          break
        }
      }
      ##########################
      cohort_num_s1 = sum(n1)/cohortsize1
      ######replace value#######
      S1_result <- generate_path_S1_values(true_matrix_tox)
      path_positions <- S1_result$positions
      for(i in seq_along(n1)) {
        N[path_positions[i, 1], path_positions[i, 2]] <- n1[i]
      }
      
      for(i in seq_along(ytox1)) {
        YTOX[path_positions[i, 1], path_positions[i, 2]] <- ytox1[i]
      }
      
      for(i in seq_along(yeff1)) {
        YEFF[path_positions[i, 1], path_positions[i, 2]] <- yeff1[i]
      }
      
      for (i in seq_along(ytox1)) {
        joint_TE_array[path_positions[i, 1], path_positions[i, 2],]<-ytox_eff1[1,i,]
      }
      
      
      ########################
      tox_location1 <- path_positions[tox_dose,]
      
      
      s1_indices <- get_S1_indices(nrow(true_matrix_tox), 
                                   ncol(true_matrix_tox))
      toxic_pos <- tox_location1
      
      if(length(toxic_pos)==0){
        union_indices<- unique(rbind(s1_indices, s1_indices), by = c("Var1", "Var2"))
      }else{
        toxic_indices <- exclude_toxic_doses_indices(nrow(true_matrix_tox), 
                                                     ncol(true_matrix_tox), toxic_pos)
        # Union of S1 and toxic doses
        union_indices <- unique(rbind(s1_indices, toxic_indices), by = c("Var1", "Var2"))
      }
      
      # Removing the union from original matrix
      resulting_matrix_tox1 <- remove_doses(true_matrix_tox, union_indices)
      resulting_matrix_eff1 <- remove_doses(true_matrix_eff, union_indices)
      
      subtrial_matrix_tox <- compress_matrix(resulting_matrix_tox1)
      subtrial_matrix_eff <- compress_matrix(resulting_matrix_eff1)
    }
    ############################
    
    for (subtrial in 2:ntrial) {
      S1=0
      
      if(length(subtrial_matrix_tox)==0||cohort_num_s1>=max_cor_size_s1){
        break
      }else{
        
        tox_values_sub <- generate_path_S1_values(subtrial_matrix_tox)$values
        eff_values_sub <- generate_path_S1_values(subtrial_matrix_eff)$values
        
        N_sub<-matrix(0,nrow = nrow(subtrial_matrix_tox),ncol = ncol(subtrial_matrix_tox))
        YTOX_sub<-matrix(0,nrow = nrow(subtrial_matrix_tox),ncol = ncol(subtrial_matrix_tox))
        YEFF_sub<-matrix(0,nrow = nrow(subtrial_matrix_eff),ncol = ncol(subtrial_matrix_eff))
        joint_TE_array_sub = array(0, dim = c(nrow(subtrial_matrix_tox), 
                                              ncol(subtrial_matrix_tox), 4))
        
        
        d = 1
        ytox_sub = rep(0,length(tox_values_sub))
        yeff_sub = rep(0,length(eff_values_sub))
        n_sub = rep(0,length(tox_values_sub))
        ytox_eff_sub = array(0, dim = c(nrow(tox_values_sub), 
                                        ncol(tox_values_sub), 4))
        
        
        ncorhort_sub = min(length(tox_values_sub),max_cor_size_s1-cohort_num_s1)
        tox_dose_sub = 0
        
        for (i_sub in 1:ncorhort_sub){
          
          ind_outcomes_sub<- generate_multinomial_outcome2(tox_values_sub[,d], 
                                                           eff_values_sub[,d],
                                                           rho, cohortsize1)
          
          ytox_sub[d] = ytox_sub[d] + sum(ind_outcomes_sub[3,],ind_outcomes_sub[4,])
          yeff_sub[d] = yeff_sub[d] + sum(ind_outcomes_sub[2,],ind_outcomes_sub[4,])
          ytox_eff_sub[,d,]<-ind_outcomes_sub
          n_sub[d] = n_sub[d] + cohortsize1
          
          pT_hat = ytox_sub[d]/n_sub[d]
          cohort_num_s1 = cohort_num_s1 + 1
          
          if(pT_hat<phi){
            d = d + 1
          }else{
            tox_dose_sub = d
            break
          }
        }
        
        
        sub_result <- generate_path_S1_values(subtrial_matrix_tox)
        path_positions <- sub_result$positions
        for(i in seq_along(n_sub)) {
          N_sub[path_positions[i, 1], path_positions[i, 2]] <- n_sub[i]
        }
        
        for(i in seq_along(ytox_sub)) {
          YTOX_sub[path_positions[i, 1], path_positions[i, 2]] <- ytox_sub[i]
        }
        
        
        for(i in seq_along(yeff_sub)) {
          YEFF_sub[path_positions[i, 1], path_positions[i, 2]] <- yeff_sub[i]
        }
        
        for (i in seq_along(ytox_sub)) {
          joint_TE_array_sub[path_positions[i, 1], path_positions[i, 2],]<-ytox_eff_sub[1,i,]
        }
        
        
        ######################
        tox_location_sub <- path_positions[tox_dose_sub,]
        
        ssub_indices <- get_S1_indices(nrow(subtrial_matrix_tox),
                                       ncol(subtrial_matrix_tox))
        toxic_pos_sub <- tox_location_sub
        
        
        if(length(toxic_pos_sub)==0){
          union_indices_sub <- unique(rbind(ssub_indices, ssub_indices), by = c("Var1", "Var2"))
        }else{
          toxic_indices_sub <- exclude_toxic_doses_indices(nrow(subtrial_matrix_tox),
                                                           ncol(subtrial_matrix_tox),
                                                           toxic_pos_sub)
          # Union of S1 and toxic doses
          union_indices_sub <- unique(rbind(ssub_indices, toxic_indices_sub), by = c("Var1", "Var2"))
        }
        
        
        # Find non-NA values in resulting_matrix_tox
        non_na_indices <- which(!is.na(resulting_matrix_tox1), arr.ind = TRUE)
        
        
        
        # Removing the union from original matrix
        resulting_matrix_tox <- remove_doses(subtrial_matrix_tox,
                                             union_indices_sub)
        resulting_matrix_eff <- remove_doses(subtrial_matrix_eff,
                                             union_indices_sub)
        
        
        
        subtrial_matrix_tox <- compress_matrix(resulting_matrix_tox)
        subtrial_matrix_eff <- compress_matrix(resulting_matrix_eff)
        
        resulting_matrix_tox1[non_na_indices] <- resulting_matrix_tox
        
        
        N = assign_submatrix(N, N_sub)
        YEFF = assign_submatrix(YEFF, YEFF_sub)
        YTOX = assign_submatrix(YTOX, YTOX_sub)
        
        joint_TE_array = assign_arrary(joint_TE_array,joint_TE_array_sub)
        
        
      }
    }
    
    
    ####### Construction of admissible set by end of stage I #####
    P_tox_rate_I<-YTOX/N    
    
    location_T<-which(P_tox_rate_I < phi & !is.nan(P_tox_rate_I), arr.ind = TRUE)
    
    Adm_eff<-pbeta(phiE,YEFF+alpha_E,N-YEFF+beta_E,lower.tail = F)
    location_E<-which(Adm_eff > tau_E & !is.nan(Adm_eff), arr.ind = TRUE)
    ### admissible set ###
    intersection<-merge(location_T, location_E, by = c("row", "col"))
    
    if(nrow(intersection)!=0){
      ########## Assign N2 patients in stage II ######
      #S2_ind_size<-divide_equally(maxn2, nrow(intersection))
      S2_ind_size<-divide_equally(maxn2+(maxn1-sum(N)), nrow(intersection))
      S2_outcome<-get_stage2_outcome(S2_ind_size,true_matrix_tox ,
                                     true_matrix_eff ,intersection)
      
      YEFF2 = YEFF+S2_outcome$YE2_matrix
      N2 = N +S2_outcome$N2_matrix
      YTOX2 = YTOX +S2_outcome$YT2_matirx
      
      joint_TE_array_SS2 = joint_TE_array+S2_outcome$joint_TE_arrary2
      
      EST_TOX2<-model_T(N2,YTOX2)
      
      refine_rate<-function(matrix_in){
        for (i in 1:nrow(matrix_in)) {
          matrix_in[i,] <- matrix_in[i,] + 
            seq(0, 0.0000001, length.out = ncol(matrix_in))
        }
        
        # Increment by rows for each column
        for (j in 1:ncol(matrix_in)) {
          matrix_in[,j] <- matrix_in[,j] + 
            seq(0.0000001, 0, length.out = nrow(matrix_in))
        }
        
        return(matrix_in)
      }
      
      EST_TOX2_new<-refine_rate(EST_TOX2)
      
      find_min_value<-function(matrix_in, target){
        diff_matrix <- abs(matrix_in - target)
        
        # Find the position (i, j) of the minimum value
        position <- which(diff_matrix == min(diff_matrix), arr.ind = TRUE)
        
        # Retrieve the closest value to target_T from the original matrix
        closest_value <- matrix_in[position[1, "row"], position[1, "col"]]
        
        return(closest_value)
      }
      
      phi_T2 = find_min_value(EST_TOX2_new,phiT)
      
      
      location_T2<-which(EST_TOX2_new <= phi_T2 & !is.nan(EST_TOX2_new), arr.ind = TRUE)
      
      Adm_eff2<-pbeta(phiE,YEFF2+alpha_E,N2-YEFF2+beta_E,lower.tail = F)
      location_E2<-which(Adm_eff2 >tau_E & !is.nan(Adm_eff2), arr.ind = TRUE)
      
      intersection2<-merge(location_T2, location_E2, by = c("row", "col"))
      
      rate_array<-get_est_dose(joint_TE_array_SS2,N2)
      Mean_Uti_est = Uti_EST(Uti_score,rate_array)/100
      
      if (nrow(intersection2)!=0){
        high_locationS2<-find_high_values_in_intersection(Mean_Uti_est, intersection2,gamma)
        A3<-merge(intersection2, high_locationS2, 
                  by = c("row", "col"))
        S3_ind_size<-divide_equally(maxn3,nrow(A3))
        
        S3_outcome<-get_stage3_outcome_TES(S3_ind_size,true_matrix_tox,true_matrix_eff,true_matrix_surv,
                                           surv_dist,alpha,maxt,A3, N2)
        
        YEFF3 = YEFF2+S3_outcome$YE3_matrix
        N3 = N2 +S3_outcome$N3_matrix
        YTOX3 = YTOX2 +S3_outcome$YT3_matirx
        
        joint_TE_array_SS3 = joint_TE_array_SS2+S3_outcome$joint_TE_arrary3
        
        EST_TOX3<-model_T(N3,YTOX3)
        
        EST_TOX3_new<-refine_rate(EST_TOX3)
        
        phi_T3 = find_min_value(EST_TOX3_new,phiT)
        
        location_T3<-which(EST_TOX3_new <= phi_T3 & !is.nan(EST_TOX3_new), arr.ind = TRUE)
        
        intersection3<-merge(location_T3, A3, by = c("row", "col"))
        
        surv_outcome<-S3_outcome$surv_array
        
        surv_est <- get_surv_est(surv_outcome, maxt)
        
        YS <- matrix(0, nrow = nrow(true_matrix_surv), ncol = ncol(true_matrix_surv))
        
        if(nrow(intersection3)==0){
          dopt = data.frame(row=0,col=0)
        } else{
          
          for(i in 1:nrow(A3)) {
            YS[A3$row[i], A3$col[i]] <- surv_est[i]
          }
          
          extracted_values_surv <- sapply(1:nrow(intersection3), function(i) {
            YS[intersection3$row[i], intersection3$col[i]]
          })
          
          max_value_surv <- max(extracted_values_surv)
          
          if (max_value_surv <= phi_s){
            dopt = data.frame(row=0,col=0)
          }else{
            #max_idx <- which.max(extracted_values_surv)
            # Corresponding location in intersection2
            max_idx<-which(extracted_values_surv==max_value_surv)
            max_location <- intersection3[max_idx,]
            if(nrow(max_location)>1){
              rate_array2<-get_est_dose(joint_TE_array_SS3,N3)
              Mean_Uti_est2 = Uti_EST(Uti_score,rate_array2)/100
              high_locationS3<-find_high_values_in_intersection(Mean_Uti_est2, 
                                                                max_location,1)
              new_location<-merge(max_location, high_locationS3, 
                                  by = c("row", "col"))
              locations = new_location
              
              dopt = locations[1,]
              
              
            }else{
              dopt = max_location
            }
          }
        }
      }else{
        N3 = N2
        YEFF3 = YEFF2
        YTOX3 = YTOX2
        dopt = data.frame(row=0,col=0)
        joint_TE_array_SS3 = joint_TE_array_SS2
      }
    } else{
      
      N3 = N
      YEFF3 = YEFF
      YTOX3 = YTOX
      joint_TE_array_SS3 = joint_TE_array
      Mean_Uti_est = matrix(0,nrow = nrow(true_matrix_tox),ncol = ncol(true_matrix_tox))
      dopt = data.frame(row=0,col=0)
    }
    
    
    ##########
    
    
    
    simulation_N[,,trial] <- N3
    #simulation_Tox[,,trial] <- YTOX3
    #simulation_Eff[,,trial] <- YEFF3
    simulation_Tox[,,trial] = joint_TE_array_SS3[,,3]+joint_TE_array_SS3[,,4]
    simulation_Eff[,,trial] = joint_TE_array_SS3[,,2]+joint_TE_array_SS3[,,4]
    simulation_Uti[,,trial] = Mean_Uti_est
    d_opt[trial,] = dopt
  }
  matrix_sums <- apply(simulation_N, 3, sum)
  max_index <- which.max(matrix_sums)
  max_matrix <- simulation_N[,,max_index]
  
  result<-list(N_ave=apply(simulation_N, c(1,2), mean),
              Tox_ave=apply(simulation_Tox, c(1,2), mean),
              Eff_ave=apply(simulation_Eff, c(1,2), mean),
              max_N=max_matrix,
              Uti_ave=apply(simulation_Uti, c(1,2), function(x) mean(x, na.rm = TRUE)),
              True_Utility = True_Uti,
              dopt=d_opt)
  counts<-result$dopt
  result.list <-list(Sel = round(table(counts)/sum(table(counts))*100,1),
                     Patients_allocation_N=result$N_ave,
                     Total_N=sum(result$N_ave),
                     Patients_allocation_Percentage=round(result$N_ave/sum(result$N_ave)*100,1),
                     True_Utility=result$True_Utility);result.list
  return(result.list)

}

##########################################################################################################
##########################################################################################################
############           Function to get trial implementation of the Great Wall design  

##########################################################################################################
##########################################################################################################
## Function to get dose escalation decision for incoming cohort of the stage I of Great Wall design 
## ytox: data for toxicity
## n: cohortsize
## phiT: highest acceptable toxicity rate
## rhoT: dose deescalation boundary

Decision_stageI<-function(ytox,n,phiT,rhoT){
  get_phi_boundary<-function(phi_T,rho){
    return(log((1-phi_T)/(1-rho*phi_T))/log(rho*((1-phi_T)/(1-rho*phi_T))))
  }
  pT_hat = ytox/n
  phi = get_phi_boundary(phiT,rhoT)
  
  if(pT_hat<phi){
    d = "Escalate"
  }else{
    d = "overly toxic"
  }
  return(d)
}

##########################################################################################################
##########################################################################################################
## Function to get dose finding path at stage I  
## mat: dose combination matrix

generate_path_S1_values <- function(mat) {
  nrow_mat <- nrow(mat)
  ncol_mat <- ncol(mat)
  
  # Start from bottom left and go up the first column
  path <- cbind(seq(nrow_mat, 1, by = -1), rep(1, nrow_mat))
  
  # Then move horizontally to the top right from the first row
  if (ncol_mat > 1) {
    path <- rbind(path, cbind(rep(1, ncol_mat - 1), seq(2, ncol_mat)))
  }
  
  # Extract values from matrix using the path indices
  values <- mat[cbind(path[,1], path[,2])]
  values_matrix <- matrix(values, nrow = 1)
  result<-list(path1 = values_matrix, positions = path)
  return(result$path1)
}

##########################################################################################################
##########################################################################################################
## Function to get decision by the end of stage I  
## YTOX: toxicity data at each combination by the end of stage I
## YEFF: efficacy data at each combination by the end of stage I
## N: enrolled sample size at each combination by the end of stage I
## phiT: highest acceptable toxicity rate
## rhoT: dose deescalation boundary
## phiE: lower limit of efficacy rate 
## tau_E: cutoff probability for futile combination
## alpha_E: hyper parameter for beta prior of efficacy 
## beta_E: hyper parameter for beta prior of efficacy 

decision_end_stageI<-function(YTOX,YEFF,N,phiT,rhoT,phiE,alpha_E,beta_E,tau_E){
  P_tox_rate_I<-YTOX/N    
  get_phi_boundary<-function(phi_T,rho){
    return(log((1-phi_T)/(1-rho*phi_T))/log(rho*((1-phi_T)/(1-rho*phi_T))))
  }
  phi = get_phi_boundary(phiT,rhoT)
  location_T<-which(P_tox_rate_I < phi & !is.nan(P_tox_rate_I), arr.ind = TRUE)
  
  Adm_eff<-pbeta(phiE,YEFF+alpha_E,N-YEFF+beta_E,lower.tail = F)
  location_E<-which(Adm_eff > tau_E & !is.nan(Adm_eff), arr.ind = TRUE)
  ### admissible set ###
  intersection<-merge(location_T, location_E, by = c("row", "col"))
  extracted_values <- mapply(function(r, c) dose_matrix[r, c], intersection$row, intersection$col)
  return(extracted_values)
}



##########################################################################################################
##########################################################################################################
## Function to get decision by the end of stage II  
## YTOX2: toxicity data at each combination by the end of stage II
## YEFF2: efficacy data at each combination by the end of stage II
## N2: enrolled sample size at each combination by the end of stage II
## joint_TE_array_SS2: joint toxicity and efficacy data at each combination by the end of stage II with 
##                     format (T, E): (0,0), (0,1), (1,0), (1,1)
## phiT: highest acceptable toxicity rate
## rhoT: dose deescalation boundary
## phiE: lower limit of efficacy rate 
## tau_E: cutoff probability for futile combination
## alpha_E: hyper parameter for beta prior of efficacy 
## beta_E: hyper parameter for beta prior of efficacy 
## gamma:  flexibility of candidate set selection at end of stage II
## Utility_score: utility score for each possible outcome [(T,E): (0,0), (0,1), (1,0), (1,1)]

decision_end_stageII<-function(YTOX2,N2,YEFF2,joint_TE_array_SS2,phiT,rhoT,phiE,alpha_E,beta_E,tau_E,
                               Utility_score,gamma){
  library(Iso)
  model_T <- function(nti,yti) {
    pt = yti/nti
    pt[is.na(pt)] = 0
    pt_flipped <- pt[nrow(pt):1, ]
    nti_flipped<- nti[nrow(nti):1, ]
    rate_flipped<-biviso(y=pt_flipped,w=nti_flipped,fatal = F,warn = F)
    T_est=rate_flipped[nrow(rate_flipped):1,]
    return(T_est)
  }
  
  get_est_dose<-function(input_arrary,N_matrix){
    rate_array <- array(0, dim = dim(input_arrary))
    
    # Loop through the third dimension and add corresponding matrices
    for (k in 1:dim(input_arrary)[3]) {
      rate_array[,,k] <- input_arrary[,,k]/N_matrix
      rate_array[,,k][is.na(rate_array[,,k])] = 0
    }
    return(rate_array)
  }
  ## (T,E): (0,0), (0,1), (1,0), (1,1)
  Uti_score<-array(Utility_score,c(1,1,4))
  Uti_EST<-function(Utility,prop_est){
    Uti_est <- matrix(0, nrow = dim(prop_est)[1],ncol = dim(prop_est)[2])
    result_array <- array(0, dim(prop_est))
    
    for (k in 1:dim(prop_est)[3]) {
      result_array[,,k] <- prop_est[,,k] * Utility[,,k]
    }
    Uti_est = result_array[,,1]+result_array[,,2]+
      result_array[,,3]+result_array[,,4]
    
    return(Uti_est)
  }
  find_high_values_in_intersection <- function(matrix, intersection,gamma) {
    extract_values <- function(matrix, intersection){
      sapply(1:nrow(intersection), function(i) {
        matrix[intersection[i, "row"], intersection[i, "col"]]
      })
    }
    
    # Use the function to extract values
    extracted_values <- extract_values(matrix, intersection2)
    
    max_value <- max(extracted_values)
    
    # Identify locations in Mean_Uti_est with values >= gamma * max_value
    identify_locations <- function(matrix, threshold){
      which(matrix >= threshold, arr.ind = TRUE)
    }
    
    locations <- identify_locations(matrix, gamma * max_value)
    return(locations)
  }
  
  closest_value_and_position <- function(mat, target = 0.3) {
    # Subtract target from each element
    diff_matrix <- abs(mat - target)
    
    # Find the positions of all the minimum values
    positions <- which(diff_matrix == min(diff_matrix), arr.ind = TRUE)
    
    # Sort by column first (to get the most left) and then by row in decreasing order (to get the most bottom)
    sorted_positions <- positions[order(positions[, "col"], -positions[, "row"]), ]
    
    # Extract the most left and most bottom position
    desired_position <- sorted_positions[1, ]
    
    # Retrieve the closest value to target from the original matrix
    closest_value <- mat[desired_position["row"], desired_position["col"]]
    
    list(value = closest_value, position = desired_position)
  }
  
  refine_rate<-function(matrix_in){
    for (i in 1:nrow(matrix_in)) {
      matrix_in[i,] <- matrix_in[i,] + 
        seq(0, 0.0000001, length.out = ncol(matrix_in))
    }
    
    # Increment by rows for each column
    for (j in 1:ncol(matrix_in)) {
      matrix_in[,j] <- matrix_in[,j] + 
        seq(0.0000001, 0, length.out = nrow(matrix_in))
    }
    
    return(matrix_in)
  }
  
  find_min_value<-function(matrix_in, target){
    diff_matrix <- abs(matrix_in - target)
    
    # Find the position (i, j) of the minimum value
    position <- which(diff_matrix == min(diff_matrix), arr.ind = TRUE)
    
    # Retrieve the closest value to target_T from the original matrix
    closest_value <- matrix_in[position[1, "row"], position[1, "col"]]
    
    return(closest_value)
  }
  
  find_closest_and_less_indices <- function(matrix, target) {
    closest_values <- apply(matrix, 1, function(row) {
      # Find the value closest to the target in each row
      row[which.min(abs(row - target))]
    })
    
    # Matrix to store results of comparison for less than closest value
    less_than_indices <- matrix(, nrow = 0, ncol = 2)
    colnames(less_than_indices) <- c("row", "col")
    
    # Identify indices where elements are less than the closest value
    for (i in 1:nrow(matrix)) {
      for (j in 1:ncol(matrix)) {
        if (matrix[i, j] <= closest_values[i]) {
          less_than_indices <- rbind(less_than_indices, c(i, j))
        }
      }
    }
    
    # Convert the indices matrix to integer
    less_than_indices <- apply(less_than_indices, 2, as.integer)
    
    list(closest_values = closest_values, less_than_indices = less_than_indices)
  }
  
  EST_TOX2<-model_T(N2,YTOX2)
  EST_TOX2_new<-refine_rate(EST_TOX2)
  
  
  
  #phi_T2 = find_min_value(EST_TOX2_new,phiT)
  
  #location_T2<-which(EST_TOX2_new <= phi_T2 & !is.nan(EST_TOX2_new), arr.ind = TRUE)
  result <- find_closest_and_less_indices(EST_TOX2_new, phiT)
  
  
  location_T2<-result$less_than_indices
  
  Adm_eff2<-pbeta(phiE,YEFF2+alpha_E,N2-YEFF2+beta_E,lower.tail = F)
  location_E2<-which(Adm_eff2 >tau_E & !is.nan(Adm_eff2), arr.ind = TRUE)
  
  intersection2<-merge(location_T2, location_E2, by = c("row", "col"))
  
  rate_array<-get_est_dose(joint_TE_array_SS2,N2)
  Mean_Uti_est = Uti_EST(Uti_score,rate_array)/100
  
  high_locationS2<-find_high_values_in_intersection(Mean_Uti_est, intersection2,gamma)
  A3<-merge(intersection2, high_locationS2, 
            by = c("row", "col"))
  
  mapply(function(r, c) dose_matrix[r, c], A3$row, A3$col)
}


##########################################################################################################
##########################################################################################################
## Function to get survival estimates by the end of stage III  
## surv_outcome: survival outcome at each combination by the end of stage III
## ass_window: follow up time for survival outcomes
get_surv_est<-function(surv_outcome, ass_window){

    # Extract time and status
    time <- surv_outcome[, 1]
    status <- surv_outcome[, 2]
    
    # Perform the survival analysis
    fit <- survfit(Surv(time, status) ~ 1)
    
    # Extract summary at time = 90
    result <- summary(fit, times = ass_window,extend = T)$surv

  
  return(result)
}


##########################################################################################################
##########################################################################################################
## Function to get decision by the end of stage III  
## YTOX3: toxicity data at each combination by the end of stage III
## N3: enrolled sample size at each combination by the end of stage III
## YS_est: survival estimates at each combination by the end of stage III
## phiT: highest acceptable toxicity rate
## rhoT: dose deescalation boundary
## phi_s: survival measurement lower limit
decision_end_stageIII<-function(YTOX3,N3,YS_est,phiT,rhoT,phi_s){
  library(Iso)
  model_T <- function(nti,yti) {
    pt = yti/nti
    pt[is.na(pt)] = 0
    pt_flipped <- pt[nrow(pt):1, ]
    nti_flipped<- nti[nrow(nti):1, ]
    rate_flipped<-biviso(y=pt_flipped,w=nti_flipped,fatal = F,warn = F)
    T_est=rate_flipped[nrow(rate_flipped):1,]
    return(T_est)
  }
  
  refine_rate<-function(matrix_in){
    for (i in 1:nrow(matrix_in)) {
      matrix_in[i,] <- matrix_in[i,] + 
        seq(0, 0.0000001, length.out = ncol(matrix_in))
    }
    
    # Increment by rows for each column
    for (j in 1:ncol(matrix_in)) {
      matrix_in[,j] <- matrix_in[,j] + 
        seq(0.0000001, 0, length.out = nrow(matrix_in))
    }
    
    return(matrix_in)
  }
  
  find_closest_and_less_indices <- function(matrix, target) {
    closest_values <- apply(matrix, 1, function(row) {
      # Find the value closest to the target in each row
      row[which.min(abs(row - target))]
    })
    
    # Matrix to store results of comparison for less than closest value
    less_than_indices <- matrix(, nrow = 0, ncol = 2)
    colnames(less_than_indices) <- c("row", "col")
    
    # Identify indices where elements are less than the closest value
    for (i in 1:nrow(matrix)) {
      for (j in 1:ncol(matrix)) {
        if (matrix[i, j] <= closest_values[i]) {
          less_than_indices <- rbind(less_than_indices, c(i, j))
        }
      }
    }
    
    # Convert the indices matrix to integer
    less_than_indices <- apply(less_than_indices, 2, as.integer)
    
    list(closest_values = closest_values, less_than_indices = less_than_indices)
  }
  
 
  EST_TOX3<-model_T(N3,YTOX3)
  
  EST_TOX3_new<-refine_rate(EST_TOX3)
  result <- find_closest_and_less_indices(EST_TOX3_new, phiT)

  
  location_T3<-result$less_than_indices
  
  extracted_values <- apply(location_T3, 1, function(idx) YS_est[idx[1], idx[2]])
  
  max_value_surv <- max(extracted_values)
  if (max_value_surv <= phi_s){
    dopt = data.frame(row=0,col=0)
  }else{
    max_idx<-which(extracted_values==max_value_surv)
    max_location <- location_T3[max_idx,]
  }
  
  return(mapply(function(r, c) dose_matrix[r, c], max_location["row"], max_location["col"]))
}


########################################## Examples #######################################

####### generate operating characteristics of the Great Wall design ################################
## Scenario 3
# Assume that the target toxicity rate is 30%, the target efficacy rate is 25%, and lower limit 
# of survival is 0.3.

dose_matrix_tox <- matrix(c(0.2,0.4,0.5,
                            0.1,0.35,0.45),byrow = T,nrow = 2,ncol = 3)

dose_matrix_eff <- matrix(c(0.53,0.4,0.4,
                            0.55,0.6,0.6),byrow = T,nrow = 2,ncol = 3)

dose_matrix_surv<- matrix(c(0.35,0.25,0.4,
                            0.55,0.35,0.2),byrow = T,nrow = 2,ncol = 3)



get.Great_Wall.oc(18,36,20, 0.3,1.4,dose_matrix_tox,dose_matrix_eff,
                   dose_matrix_surv,1,0.5,180,
                   0.25,0.05,1,1,
                   3,10000,0.5,0.7,0.3, c(40,100,0,60))
# $Sel
# col
# row    0    1    2    3
# 0  3.5  0.0  0.0  0.0
# 1  0.0  6.8  0.1  0.0
# 2  0.0 83.8  5.8  0.1
# 
# $Patients_allocation_N
# [,1]    [,2]   [,3]
# [1,] 18.0864  7.7077 3.8072
# [2,] 22.9295 13.6569 5.8314
# 
# $Total_N
# [1] 72.0191
# 
# $Patients_allocation_Percentage
# [,1] [,2] [,3]
# [1,] 25.1 10.7  5.3
# [2,] 31.8 19.0  8.1
# 
# $True_Utility
# [,1] [,2] [,3]
# [1,] 0.638 0.48 0.44
# [2,] 0.690 0.62 0.58

####### generate decision of the Great Wall design based on data in trial implementation section #####
########## The data is generated based on Scenario 3 #################
dose_matrix<-matrix(c("2,1","2,2","2,3","1,1","1,2","1,3"),2,3,byrow = T)
####### Get first sub-path ########
generate_path_S1_values(dose_matrix)
#[,1]  [,2]  [,3]  [,4] 
#[1,] "1,1" "2,1" "2,2" "2,3"
######### decision at A1B1 #######
Decision_stageI(1,3,0.3,1.4)
######### decision at A2B2 #######
Decision_stageI(2,3,0.3,1.4)
######### Generate sub-path 2 #######
##### First excluding those too toxic combination and get reduced dose space ####
dose_matrix_2<-matrix(c("1,2","1,3"),1,2,byrow = T)
##### get sub-path 2 ####
generate_path_S1_values(dose_matrix_2)
########## Get decision by the end of stage I based on accumluated data #####
N<-matrix(c(3,3,0,3,3,3),2,3,byrow = T)
YTOX<-matrix(c(1,2,0,0,0,2),2,3,byrow = T)
YEFF<-matrix(c(2,1,0,2,3,1),2,3,byrow = T)
decision_end_stageI(YTOX,YEFF,N,0.3,1.4,0.25,1,1,0.05)
# [1] "2,1" "1,1" "1,2"
# those three combinations will be moved forward to stage II ####


######### Get decision by the end of stage II based on accumluated data #####
YEFF2 = matrix(c(6,1,0,5,8,1),2,3,byrow = T)
N2 = matrix(c(10,3,0,10,10,3),2,3,byrow = T)
YTOX2 = matrix(c(2,2,0,0,1,2),2,3,byrow = T)
joint_TE_array_SS2 <- array(c(4, 5, 1, 2, 0, 1, 
                              4, 5, 0, 7, 0, 0, 
                              0, 0, 1, 0, 0, 1,
                              2, 0, 1, 1, 0, 1), 
                            dim = c(2, 3, 4))

decision_end_stageII(YTOX2,N2,YEFF2,joint_TE_array_SS2,0.3,1.4,0.25,1,1,0.05, c(40,100,0,60),
                     0.7)
# [1] "2,1" "1,1" "1,2"
# those three combinations will be moved forward to stage III ####

######### Get decision by the end of stage II based on accumluated data #####
######### first get survival estimates at each dose combination ######
### For example at A2B1 #########
get_surv_est(matrix(c(131, 180, 137,  40, 178, 161, 180, 150, 180, 180,
                      59, 180, 83, 147, 95, 180,
                      1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0), ncol = 2),180)

## 0.375 ##
### For example at A1B1 #########
get_surv_est(matrix(c(44, 180, 180, 180, 128, 94, 180, 92, 9, 180,
                      180, 180, 87, 180, 136, 148,
                      1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1), ncol = 2),180)

## 0.5 ##
##### ...... #####
YS_est<-matrix(c(0.375,0,0,0.5,0.3125,0),2,3,byrow = T)

YEFF3<-matrix(c(5,0,0,4,2,0),2,3,byrow = T)+YEFF2
N3<-matrix(c(6,0,0,6,6,0),2,3,byrow = T)+N2
YTOX3<-matrix(c(1,0,0,1,2,0),2,3,byrow = T)+YTOX2

decision_end_stageIII(YTOX3,N3,YS_est,0.3,1.4,0.3)
# "1,1"
### Therefore, the ODC is combination A1B1
