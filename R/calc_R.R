

fix_beta_large <- function(params, S0, I, R0, use_eig=FALSE, beta=NULL, symp_trans=NULL){
  S0 <- as.numeric(S0)/params$beta_norm
  params$large_mixing_matrix <- get_lmm(params)
  ng <- get_next_gen_large(params, S0, beta, symp_trans=symp_trans)
  if(sum(I)==0) use_eig <- T
  
  if(use_eig){
    max_eigen <- Re(eigen(ng, only.values=T)$values[1])
    max_eigen
    return(R0/max_eigen)
  }else{
    le <- sum(ng %*% as.numeric(I)) / sum(I)#R_from_NGM(ng, I)
    weight <- 3
    return( R0/ le)
  }
}




estimate_Rt_large_new <- function( traj, params,beta){
  colsS <-paste0("S[", 1:(params$n*params$n_vac), "]")
  colsI <- paste0("I[", 1:(params$n*params$n_vac*params$n_strain), "]")
  params$large_mixing_matrix <- get_lmm(params)
  n <- params$n*params$n_vac
  Rts <- c()
  for(t in 1:nrow(traj)){
    S <- unlist(traj[t,] %>% dplyr::select(colsS))
    I <- unlist(traj[t,] %>%dplyr::select(colsI))
    S <- S /params$beta_norm
    new_I <- 0
    beta <- params$beta_day
    for(i_strain in 1:params$n_strain){
      ng <- get_next_gen_large(params, S, beta[unlist(traj[t, "t"])]*
                                          params$beta_strain[i_strain], i_strain)
      new_I <- new_I + sum(ng %*% I[((i_strain - 1)*n):(i_strain*n-1)+1])
    }
    Rts <- c(Rts, new_I/ sum(I))
  }
  return(data.frame(Rt=Rts))
}
  
se1e2iiaR_calculate_beta_duration <- function(
                                             a2 = a2,
                                             gamma = gamma,
                                             presymptomaticRelativeInfectiousness = presymptomatic_relative_infectiousness,
                                             asymptomaticProb = asymptomatic_prob,
                                             asymptomaticRelativeInfectiousness = asymptomatic_relative_infectiousness, symp_trans=1) {
  d <- presymptomaticRelativeInfectiousness * (1 - asymptomaticProb) / a2 +
    (1 - asymptomaticProb) * symp_trans/ gamma +
    (asymptomaticRelativeInfectiousness * asymptomaticProb) / gamma
  return(d)
}


get_next_gen_large <- function(params, S, beta, i_strain, symp_trans=NULL){
  if(is.null(symp_trans)) symp_trans <- as.numeric(params$symp_trans[,,i_strain])
  asymp_frac <- as.numeric(params$asympt_frac[,, i_strain]*params$susceptibility_asymp[,, i_strain]/(params$asympt_frac[,, i_strain]*params$susceptibility_asymp[,, i_strain] + (1- params$asympt_frac[,, i_strain])*params$susceptibility_symp[,, i_strain]))
  vac_susc <- as.numeric(params$susceptibility_asymp[,, i_strain]*params$asympt_frac[,, i_strain] + (1-params$asympt_frac[,, i_strain])*params$susceptibility_symp[,, i_strain])
  D <- se1e2iiaR_calculate_beta_duration(1/params$pre_sympt_period, 1/params$infectious_period, params$pre_sympt_infect, asymp_frac, params$asympt_infect, symp_trans=symp_trans)
  ng <- sweep(params$large_mixing_matrix, MARGIN=1, vac_susc*D*S*rep(as.numeric(beta), params$n_vac), `*`)
  ng <- sweep(ng, MARGIN=2, as.numeric(params$transmisibility[,,i_strain]), `*`)
  return(ng)

}

