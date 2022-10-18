
add_per_n <- function(res, key, params){

  n_list <- seq(1, params$n*params$n_vac*params$n_strain,by=params$n)
  for(i in 0:(params$n-1)){
    vals <- res[, paste0(key, "[",n_list+i,"]")]
    if(!is.null(dim(vals))){
      vals <- rowSums(vals)
    }
    res <- cbind(res, vals)
    colnames(res)[ncol(res)] <- paste0(key, "_", i+1)
  }
  return(res)
}


add_per_age <- function(res, key, params){

  fact <- 1
  if(!is.null(params$merge_half_age)){
    fact <- 2
  }
  n_list <- seq(1, params$n*params$n_vac*params$n_strain,by=params$age_groups/fact)
  for(i in 0:(params$age_groups/fact-1)){
    vals <- res[, paste0(key, "[",n_list+i,"]")]
    if(!is.null(dim(vals))){
      vals <- rowSums(vals)
    }
    res <- cbind(res, vals)
    colnames(res)[ncol(res)] <- paste0(key, "_age_", i +1)
  }
  return(res)
}


add_per_vaccine_prio_group <- function(res, key, params, only_reg=FALSE){
  res <- cbind(res, add_region_group(res, key, params, params$reg_prio, only_reg=only_reg))
  colnames(res)[ncol(res)] <- paste0(key, "_plus")
  res <- cbind(res, add_region_group(res, key, params, params$reg_prio_neutral, only_reg=only_reg))
  colnames(res)[ncol(res)] <- paste0(key, "_neutral")
  res <- cbind(res, add_region_group(res, key, params, params$reg_prio_minus, only_reg=only_reg))
  colnames(res)[ncol(res)] <- paste0(key, "_minus")
  return(res)
}


add_region_group <- function(res, key, params, regions, only_reg=FALSE){
  tmp <- rep(0, nrow(res))
  for(i in regions){

    if(only_reg){

      n_list <- i
    }else{
      n_list <- c()
      for(j in 1:params$n_vac){
        for(k in 1:params$n_strain){
          n_list <- c(n_list, 1:params$age_groups + (i-1)*params$age_groups + (j-1)*params$n + (k-1)*params$n*params$n_vac)
        }
      }
    }
    vals <- res[, paste0(key, "[",n_list,"]")]
    if(!is.null(dim(vals))){
      vals <- rowSums(vals)
    }
    tmp <- tmp + vals
  }
  return(tmp)
}
  

add_per_region <- function(res, key, params){
  n_list <- seq(1, params$n*params$n_vac*params$n_strain,by=params$N_regions)
  for(i in 1:params$N_regions){
    n_list <- c()
    for(j in 1:params$n_vac){
      for(k in 1:params$n_strain){
        n_list <- c(n_list, 1:params$age_groups + (i-1)*params$age_groups + (j-1)*params$n + (k-1)*params$n*params$n_vac)
      }
    }
    vals <- res[, paste0(key, "[",n_list,"]")]
    if(!is.null(dim(vals))){
      vals <- rowSums(vals)
    }
    res <- cbind(res, vals)
    colnames(res)[ncol(res)] <- paste0(key, "_region_", i)
  }
  return(res)
}



add_per_vac <- function(res, key, params){
  for(i in 1:params$n_vac){
    n_list <- c()
    for(k in 1:params$n_strain){
      n_list <- c(n_list, 1:params$n  + (i-1)*params$n + (k-1)*params$n*params$n_vac)
    }
    vals <- res[, paste0(key, "[",n_list,"]")]
    if(!is.null(dim(vals))){
      vals <- rowSums(vals)
    }
    res <- cbind(res, vals)
    colnames(res)[ncol(res)] <- paste0(key, "_vac_", i)
  }
  return(res)
}



#'
#' @export
refine_results_odin <- function(res, params){
  d <- res
                                        #  res <- data.table(res)
  N <- nrow(res)

  n <- params$n_vac*params$n_strain*params$n
  ## vac_list <- (N_reg*N_age + 1):params$n
  ## vac_list_dose1 <- (N_reg*N_age + 1):(N_reg*N_age*2)
  ## vac_list_dose2 <- (N_reg*N_age*2 + 1):params$n
  res <- cbind(res, tot_infected=rowSums(res[, paste0("tot_infected[",1:n,"]")]))
#  res <- cbind(res, tot_imp=rowSums(res[, paste0("tot_imp[",1:n,"]")]))
 # non_vac_list <- 1:(N_reg*N_age)
  age_list <- seq(1, params$n*params$n_vac*params$n_strain,by=params$n)

  res <- add_per_age(res, "tot_infected", params)
  res <- add_per_age(res, "tot_resp", params)
  res <- add_per_age(res, "tot_hosp", params)
  res <- add_per_age(res, "D", params)
  res <- add_per_age(res, "tot_vac", params)

  if(params$n_vac >= 2){
    res <- add_per_vac(res, "tot_vac", params)
  }
  
  if(!is.null(params$N_regions) & params$N_regions > 1){
    res <- add_per_region(res, "tot_infected", params)
    res <- add_per_region(res, "tot_resp", params)
    res <- add_per_region(res, "tot_hosp", params)
    res <- add_per_region(res, "D", params)
  }
  if(params$n_vac > 1){
    res <- add_per_vac(res, "N", params)
    res <- add_per_vac(res, "tot_infected", params)
    res <- add_per_vac(res, "tot_hosp", params)
  }
  
  if(!is.null(params$reg_prio)){
    res <- add_per_vaccine_prio_group(res, "tot_infected", params)
    res <- add_per_vaccine_prio_group(res, "tot_resp", params)
    res <- add_per_vaccine_prio_group(res, "tot_hosp", params)
    res <- add_per_vaccine_prio_group(res, "D", params)
    res <- add_per_vaccine_prio_group(res, "tot_vac", params, only_reg=TRUE)

  }
  
#  res <- cbind(res, ICU_R=rowSums(res[, paste0("ICU_R[",1:n,"]")]))
#  res <- cbind(res, MISC=rowSums(res[, paste0("MISC[",1:n,"]")] +
#                                 res[, paste0("MISC_ICU_H[",1:n,"]")] 
 #                                ))
 # res <- cbind(res, MISC_ICU=rowSums(res[, paste0("MISC_ICU[",1:n,"]")]))
  
  res <- cbind(res, D=rowSums(res[, paste0("D[",1:n,"]")]))

  res <- cbind(res, tot_hosp=rowSums(res[, paste0("tot_hosp[",1:n,"]")]))
#  res <- cbind(res, tot_misc=rowSums(res[, paste0("tot_misc[",1:n,"]")]))
  res <- cbind(res, tot_prev=rowSums(res[, c(paste0("I[",1:n,"]"), paste0("A[",1:n,"]"), paste0("P[",1:n,"]"))]))
  res <- cbind(res, tot_resp=rowSums(res[, paste0("tot_resp[",1:n,"]")]))

  #vac_numbers <- rowSums(res[, paste0("tot_vac[",vac_list,"]")])
  #before_start <- vac_numbers[n_vac_start]
  #vac_from_model <- c(vac_numbers[(params$N_vac_start+1):N], rep(vac_numbers[N], params$N_vac_start)) - before_start


  #vac_numbers <- res[, paste0("tot_vac[",vac_list,"]")]
  #before_start <- vac_numbers[params$N_vac_start,]

 ##  tot_vaccinated_all_time_i <- res[, paste0("S[",vac_list_dose1,"]")] +
 ##    res[, paste0("Ea[",vac_list_dose1,"]")] +
 ##    res[, paste0("Es[",vac_list_dose1,"]")] +
 ##    res[, paste0("I[",vac_list_dose1,"]")] +
 ##    res[, paste0("A[",vac_list_dose1,"]")] +
 ##    res[, paste0("P[",vac_list_dose1,"]")] +
 ##    res[, paste0("H[",vac_list_dose1,"]")] +
 ##    res[, paste0("ICU_H[",vac_list_dose1,"]")] +
 ##    res[, paste0("ICU_R[",vac_list_dose1,"]")] +
 ##    res[, paste0("ICU_P[",vac_list_dose1,"]")] +
 ##    res[, paste0("B_D[",vac_list_dose1,"]")] +
 ##    res[, paste0("B_D_H[",vac_list_dose1,"]")] +
 ##    res[, paste0("B_D_ICU[",vac_list_dose1,"]")] +
 ##    res[, paste0("PRE_MISC[",vac_list_dose1,"]")] +
 ##    res[, paste0("MISC[",vac_list_dose1,"]")] +
 ##    res[, paste0("D[",vac_list_dose1,"]")] +
 ##    res[, paste0("R[",vac_list_dose1,"]")] +
 ##    res[, paste0("S[",vac_list_dose2,"]")] +
 ##    res[, paste0("Ea[",vac_list_dose2,"]")] +
 ##    res[, paste0("Es[",vac_list_dose2,"]")] +
 ##    res[, paste0("I[",vac_list_dose2,"]")] +
 ##    res[, paste0("A[",vac_list_dose2,"]")] +
 ##    res[, paste0("P[",vac_list_dose2,"]")] +
 ##    res[, paste0("H[",vac_list_dose2,"]")] +
 ##    res[, paste0("ICU_H[",vac_list_dose2,"]")] +
 ##    res[, paste0("ICU_R[",vac_list_dose2,"]")] +
 ##    res[, paste0("ICU_P[",vac_list_dose2,"]")] +
 ##    res[, paste0("B_D[",vac_list_dose2,"]")] +
 ##    res[, paste0("B_D_H[",vac_list_dose2,"]")] +
 ##    res[, paste0("B_D_ICU[",vac_list_dose2,"]")] +
 ##    res[, paste0("PRE_MISC[",vac_list_dose2,"]")] +
 ##    res[, paste0("MISC[",vac_list_dose2,"]")] +
 ##    res[, paste0("D[",vac_list_dose2,"]")] +
 ##    res[, paste0("R[",vac_list_dose2,"]")] 
    
 ## # vac_from_model <- rbind(vac_numbers[(params$N_vac_start+1):N,], matrix(rep(vac_numbers[N,], each= params$N_vac_start), ncol=length(vac_list)))

 ##  #res <- cbind(res, vac_from_model)
 ##  res <- cbind(res, tot_vac1=rowSums(res[, paste0("tot_vac[",vac_list,"]")]))
 ##  res <- cbind(res, tot_vac_all_time=rowSums(tot_vaccinated_all_time_i))
 ##                                        #  res <- cbind(res, tot_vac2=rep(0, N))
 ##                                        #  res <- cbind(res, tot_vac3=rep(0, N))

 ##  res <- cbind(res, tot_vac_all_time_1=rowSums(tot_vaccinated_all_time_i[,age_list[1:(N_reg*2)]]))
 ##  res <- cbind(res, tot_vac_all_time_2=rowSums(tot_vaccinated_all_time_i[,age_list[1:(N_reg*2)]+1]))
 ##  res <- cbind(res, tot_vac_all_time_3=rowSums(tot_vaccinated_all_time_i[,age_list[1:(N_reg*2)]+2]))
 ##  res <- cbind(res, tot_vac_all_time_4=rowSums(tot_vaccinated_all_time_i[,age_list[1:(N_reg*2)]+3]))
 ##  res <- cbind(res, tot_vac_all_time_5=rowSums(tot_vaccinated_all_time_i[,age_list[1:(N_reg*2)]+4]))
 ##  res <- cbind(res, tot_vac_all_time_6=rowSums(tot_vaccinated_all_time_i[,age_list[1:(N_reg*2)]+5]))
 ##  res <- cbind(res, tot_vac_all_time_7=rowSums(tot_vaccinated_all_time_i[,age_list[1:(N_reg*2)]+6]))
 ##  res <- cbind(res, tot_vac_all_time_8=rowSums(tot_vaccinated_all_time_i[,age_list[1:(N_reg*2)]+7]))
 ##  res <- cbind(res, tot_vac_all_time_9=rowSums(tot_vaccinated_all_time_i[,age_list[1:(N_reg*2)]+8]))
                                                                   

  res <- cbind(res, "I"=rowSums(res[, c(paste0("I[",1:n,"]"), paste0("A[",1:n,"]"))]))
  res <- cbind(res, "A"=rowSums(res[, paste0("A[",1:n,"]")]))
#  res <- cbind(res, "I_vac"=rowSums(res[, c(paste0("I[",vac_list,"]"), paste0("A[",vac_list,"]"))]))
#  res <- cbind(res, hosp_vac=rowSums(res[, c(paste0("H[",vac_list,"]"), paste0("ICU_P[",vac_list,"]"), paste0("ICU_H[",vac_list,"]"),  paste0("ICU_R[",vac_list,"]"))]))
  
  res <- cbind(res, ward=rowSums(res[, c(paste0("H[",1:n,"]"), paste0("ICU_P[",1:n,"]"), paste0("ICU_H[",1:n,"]"))]))
  res <- cbind(res, hosp=rowSums(res[, c("ward", paste0("ICU_R[",1:n,"]"))]))
  res <- cbind(res, resp=rowSums(res[, paste0("ICU_R[",1:n,"]")]))
  incidence <- c(res[2:N, "tot_infected"] - res[1:(N-1), "tot_infected"],0)
  res <- cbind(res, incidence=incidence)
  hosp_incidence <- c(res[2:N, "tot_hosp"] - res[1:(N-1), "tot_hosp"],0)
  res <- cbind(res, hosp_incidence=hosp_incidence)
  return(data.table::data.table(res))
}


#'
#' @export
refine_results_odin_dust <- function(res, params, N_threads){

                                        #  refined <- lapply(1:dim(res)[2], function(n)refine_one_sim(res[,n,], params, n, N_reg, N_age, N_wax))
  refined <- parallel::mclapply(1:dim(res)[2], function(n)refine_one_sim(res[,n,], params,n),
                                mc.cores=N_threads)
##  refined <- lapply(1:dim(res)[2], function(n)refine_one_sim(res[,n,], params,n))

  
  return(data.table::rbindlist(refined))
}


fix_index <- function(index){
  ind <- c()
  for(n in names(index)){
    if(length(index[[n]]) > 1){
      ind <- c(ind, paste0(n, "[", 1:length(index[[n]]),"]"))
    }else{
      ind <- c(ind, n)
    }
  }
    
  
  return(ind)

}

to_results_dt <- function(res, model, filter_dt=FALSE){
  if(class(model$info()) == "list"){
    dust_index <- model$info()[[1]]$index
  }else{
    dust_index <- model$info()$index
  }
  results <- list()
  for(i in 1:dim(res)[2]){
    res_t <- res[,i,]
    res_t <- t(res_t[1:dim(res_t)[1],])
    colnames(res_t) <- fix_index(dust_index)
    results[[length(results) + 1]] <- as.data.table(res_t) %>% dplyr::mutate(sim=i)
  }
  results <- rbindlist(results)
  if (filter_dt){
    results <- results[time %in% min(results$time):max(results$time)]
  }
  return(results)
}

refine_one_sim <- function(res, params,n){
  res <- t(res[1:dim(res)[1],])
  colnames(res) <- fix_index(params$dust_index)
  dt <- refine_results_odin(res,params) %>% dplyr::mutate(sim=n,
                                                           t=1:dim(res)[1])
  
  return(dt)
}

