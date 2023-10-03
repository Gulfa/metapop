
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


add_per_strain <- function(res, key, params){

  n_list <- seq(1, params$n*params$n_vac,by=1)
  for(i in 0:(params$n_strain-1)){
    vals <- res[, paste0(key, "[",n_list+i*params$n*params$n_vac,"]")]
    if(!is.null(dim(vals))){
      vals <- rowSums(vals)
    }
    res <- cbind(res, vals)
    colnames(res)[ncol(res)] <- paste0(key, "_strain_", i+1)
  }
  return(res)
}


add_per_age <- function(res, key, params, use_strain=TRUE, vac_index=NULL){

  fact <- 1
  if(!is.null(params$merge_half_age)){
    fact <- 2
  }
  if(use_strain){
      n_list <- seq(1, params$n*params$n_vac*params$n_strain,by=params$age_groups/fact)
      
  }else{
    if(is.null(vac_index)){
      n_list <- seq(1, params$n*params$n_vac,by=params$age_groups/fact)
    }else{
      n_list <- seq(1, params$n,by=params$age_groups/fact) + (vac_index-1)*params$n
    }
  }
  
  
  for(i in 0:(params$age_groups/fact-1)){
    vals <- res[, paste0(key, "[",n_list+i,"]")]
    if(!is.null(dim(vals))){
      vals <- rowSums(vals)
    }
    #res <- cbind(res, vals)
    #colnames(res)[ncol(res)] <- paste0(key, "_age_", i +1)
    res[, paste0(key, "_age_", i +1)] <- vals
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



add_per_vac <- function(res, key, params, use_strain=TRUE){
  for(i in 1:params$n_vac){
    n_list <- c()
    if(use_strain){
      for(k in 1:params$n_strain){
        n_list <- c(n_list, 1:params$n  + (i-1)*params$n + (k-1)*params$n*params$n_vac)
      }
    }else{
      n_list <- c(n_list, 1:params$n  + (i-1)*params$n)
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
  res <- data.table(res)
  N <- nrow(res)

  n <- params$n_vac*params$n_strain*params$n
  ## vac_list <- (N_reg*N_age + 1):params$n
  ## vac_list_dose1 <- (N_reg*N_age + 1):(N_reg*N_age*2)
  ## vac_list_dose2 <- (N_reg*N_age*2 + 1):params$n
  res <- cbind(res, tot_infected=rowSums(res[, paste0("tot_infected[",1:n,"]")]))
#  res <- cbind(res, tot_imp=rowSums(res[, paste0("tot_imp[",1:n,"]")]))
 # non_vac_list <- 1:(N_reg*N_age)
  age_list <- seq(1, params$n*params$n_vac*params$n_strain,by=params$n)

  if(!"time" %in% colnames(res)){
    res <- cbind(res, "time"=res[, "t"])
  }
  res <- add_per_age(res, "tot_infected", params)
  res <- add_per_age(res, "tot_resp", params)
  res <- add_per_age(res, "tot_hosp", params)
  res <- add_per_age(res, "D", params)
  res <- add_per_age(res, "tot_vac", params, use_strain=FALSE, vac_index=1)
  res <- add_per_age(res, "tot_vac_adm", params, use_strain=FALSE, vac_index=1)
  res <- add_per_age(res, "S", params, use_strain=FALSE)


  if(params$n_vac >= 2){
    res <- add_per_vac(res, "tot_vac", params, use_strain=FALSE)
 
  }
  if(!is.null(params$N_regions) & params$N_regions > 1){
    res <- add_per_region(res, "tot_infected", params)
    res <- add_per_region(res, "tot_resp", params)
    res <- add_per_region(res, "tot_hosp", params)
    res <- add_per_region(res, "D", params)
  }
  if(params$n_vac > 1){
    res <- add_per_vac(res, "N", params, use_strain = FALSE)
    res <- add_per_vac(res, "tot_infected", params)
    res <- add_per_vac(res, "tot_hosp", params)
    res <- add_per_vac(res, "S", params, use_strain=FALSE)
    res <- add_per_vac(res, "R", params, use_strain=FALSE)

  }
  
  if(!is.null(params$reg_prio)){
    res <- add_per_vaccine_prio_group(res, "tot_infected", params)
    res <- add_per_vaccine_prio_group(res, "tot_resp", params)
    res <- add_per_vaccine_prio_group(res, "tot_hosp", params)
    res <- add_per_vaccine_prio_group(res, "D", params)
    res <- add_per_vaccine_prio_group(res, "tot_vac_adm", params, only_reg=TRUE)

  }
  if(params$n_strain > 1){
    res <- add_per_strain(res, "tot_hosp", params)
    res <- add_per_strain(res, "tot_infected", params)
    for(i in 1:params$n_strain){
      incidence <- c(res[2:N, paste0("tot_infected_strain_",i)] - res[1:(N-1), paste0("tot_infected_strain_",i)],0)
      res <- cbind(res, incidence)
      colnames(res)[ncol(res)] <-paste0("incidence_strain_",i)
      hosp_incidence <- c(res[2:N, paste0("tot_hosp_strain_",i)] - res[1:(N-1), paste0("tot_hosp_strain_",i)],0)
      res <- cbind(res, hosp_incidence)
      colnames(res)[ncol(res)] <-paste0("hosp_incidence_strain_",i)
    }
  }
 
  res <- cbind(res, D=rowSums(res[, paste0("D[",1:n,"]")]))

  res <- cbind(res, tot_hosp=rowSums(res[, paste0("tot_hosp[",1:n,"]")]))
#  res <- cbind(res, tot_misc=rowSums(res[, paste0("tot_misc[",1:n,"]")]))
  res <- cbind(res, tot_prev=rowSums(res[, c(paste0("I[",1:n,"]"), paste0("A[",1:n,"]"), paste0("P[",1:n,"]"))]))
  res <- cbind(res, tot_resp=rowSums(res[, paste0("tot_resp[",1:n,"]")]))


  res <- cbind(res, "I"=rowSums(res[, c(paste0("I[",1:n,"]"), paste0("A[",1:n,"]"))]))
  res <- cbind(res, "A"=rowSums(res[, paste0("A[",1:n,"]")]))
  res <- cbind(res, ward=rowSums(res[, c(paste0("H[",1:n,"]"), paste0("ICU_P[",1:n,"]"), paste0("ICU_H[",1:n,"]"))]))
  res <- cbind(res, hosp=rowSums(res[, c("ward", paste0("ICU_R[",1:n,"]"))]))
  res <- cbind(res, resp=rowSums(res[, paste0("ICU_R[",1:n,"]")]))
  if(! "incidence" %in% colnames(res)){
    incidence <- c(res[2:N, "tot_infected"] - res[1:(N-1), "tot_infected"],0)
    res <- cbind(res, incidence=incidence)
  }
  hosp_incidence <- c(res[2:N, "tot_hosp"] - res[1:(N-1), "tot_hosp"],0)
  res <- cbind(res, hosp_incidence=hosp_incidence)
  return(data.table::data.table(res))
}

#'
#' @export
refine_results_odin_minimal <- function(res, params){
  d <- res
  res <- data.table(res)
  N <- nrow(res)

  n <- params$n_vac*params$n_strain*params$n
  res <- cbind(res, tot_infected=rowSums(res[, paste0("tot_infected[",1:n,"]")]))
  if(!"time" %in% colnames(res)){
    res <- cbind(res, "time"=res[, "t"])
  }
  res <- add_per_age(res, "tot_infected", params)
  res <- add_per_age(res, "tot_resp", params)
  res <- add_per_age(res, "tot_hosp", params)
  res <- add_per_age(res, "D", params)
 
 
  res <- cbind(res, D=rowSums(res[, paste0("D[",1:n,"]")]))

  res <- cbind(res, tot_hosp=rowSums(res[, paste0("tot_hosp[",1:n,"]")]))
#  res <- cbind(res, tot_misc=rowSums(res[, paste0("tot_misc[",1:n,"]")]))
  res <- cbind(res, tot_resp=rowSums(res[, paste0("tot_resp[",1:n,"]")]))
  res <- cbind(res, ward=rowSums(res[, c(paste0("H[",1:n,"]"), paste0("ICU_P[",1:n,"]"), paste0("ICU_H[",1:n,"]"))]))
  res <- cbind(res, hosp=rowSums(res[, c("ward", paste0("ICU_R[",1:n,"]"))]))
  res <- cbind(res, resp=rowSums(res[, paste0("ICU_R[",1:n,"]")]))
  if(! "incidence" %in% colnames(res)){
    incidence <- c(res[2:N, "tot_infected"] - res[1:(N-1), "tot_infected"],0)
    res <- cbind(res, incidence=incidence)
  }
  hosp_incidence <- c(res[2:N, "tot_hosp"] - res[1:(N-1), "tot_hosp"],0)
  res <- cbind(res, hosp_incidence=hosp_incidence)
  return(data.table::data.table(res))
}


#'
#' @export
refine_results_odin_dust <- function(res, params, N_threads, type="normal"){

                                        #  refined <- lapply(1:dim(res)[2], function(n)refine_one_sim(res[,n,], params, n, N_reg, N_age, N_wax))
  if(params$use_determinsitic_model){
    refined <- refine_one_sim(res, params, 1, type=type)
  }else{
  refined <- data.table::rbindlist(parallel::mclapply(1:dim(res)[2], function(n)refine_one_sim(res[,n,], params,n, type=type),
                                mc.cores=N_threads))
  }
##  refined <- lapply(1:dim(res)[2], function(n)refine_one_sim(res[,n,], params,n))

  
  return(refined)
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

fix_index_determinsitic <- function(raw_results){
    tmp <- unlist(lapply(stringr::str_split(colnames(raw_results), "\\["), function(x) x[[1]]))
    freq_count <- table(tmp)
    new_names <- c()
    used_base_names <- list()
    for(t in tmp){
      if(freq_count[t] > 1){
        count <- sum(used_base_names==t) + 1
       
        new_t <- paste(t, "[", count, "]", sep="")
      }else{
        new_t <- t
      }
      used_base_names[[length(used_base_names) + 1]] <- t
      new_names[[length(new_names) + 1]] <- new_t
    }
  return(unlist(new_names))
}

refine_one_sim <- function(res, params,n, type="normal"){
  if(!is.null(params$dust_index)){
    res <- t(res[1:dim(res)[1],])
    colnames(res) <- fix_index(params$dust_index)
  }else{
    colnames(res) <- fix_index_determinsitic(res)

  }
  if(type=="normal"){
  dt <- refine_results_odin(res,params)
  }else if(type=="minimal"){
    dt <- refine_results_odin_minimal(res,params)
  }else{
    stop("Unknown refine type")
  }
  return(dt%>% dplyr::mutate(sim=n,t=1:dim(res)[1])
  )
}

