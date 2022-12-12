#' Exporting main model
#'
#' @export model
NULL

#' run_params
#'
#' Run the metapopulation model with a set of parameters and a given number of particles and threads.
#' The number of particles gives the number of samples and threads the number of cores used to sample
#' @export
run_params <- function(params, L=200, N_particles=1, N_threads=1, run_name="run1", run_params=list(), silent=TRUE, estimate_Rt=FALSE,Rt_per_day=1,
                       return_summary_function=NULL){
  params <- fix_beta_mode_params(params)
  if(!silent){
    print(glue::glue("Running {run_name}"))
  }

  params <- fix_init_params(params)
  dust_model <- model$new(pars = params,
                         step = 1,
                         n_particles = N_particles,
                         n_threads = N_threads,
                         )

  
  raw_results <- dust_model$simulate(1:(L/params$dt))
  if(!silent){
    print(glue::glue("Finished raw run {run_name}"))
  }

  params$dust_index <- dust_model$info()$index

  if(!is.null(return_summary_function)){
    return(return_summary_function(raw_results, params) %>% mutate(name=run_name))
  }
  
  results <- refine_results_odin_dust(raw_results, params, N_threads)
  if(!silent){
    print(glue::glue("Finished refine results {run_name}"))
  }


  beta <- list()
  if(params$beta_mode == 3){
    for(n in 1:N_particles){
      beta[[length(beta) + 1]] <- outer(exp(results[sim==n, log_beta]), params$rand_beta_factors)
    }
  }
  else if(params$beta_mode == 2){
    for(n in 1:N_particles){
      beta[[length(beta) + 1]] <- outer(results[sim==n, beta_thresh], params$rand_beta_factors)
    }
  }
  ## else if(params$beta_mode == 4){
  ##   for(n in 1:N_particles){
  ##     beta[[length(beta) + 1]] <- results[sim==n, beta_spont_be
  ##   }
  ## }
  else{
      for(n in 1:N_particles){
        beta[[length(beta) + 1]] <- params$beta_day
      }
      
  }


  if(estimate_Rt){
    Rts <- parallel::mclapply(unique(results$sim), function(x) estimate_Rt_large_new(results%>%dplyr::filter(sim==x), params, beta[[x]],Rt_per_day=Rt_per_day), mc.cores=N_threads)
    results <- results %>% dplyr::mutate(Rt=rep(unlist(Rts)))
  }
  results <- results  %>% dplyr::mutate(name=run_name) %>% dplyr::mutate(!!!run_params)
  return(results)
}



#' Function to run multiple param sets
#'
#' @export
run_param_sets <- function(paramsets, L=100, N_particles=1, N_threads_internal=1, N_threads_external=1, silent=TRUE, return_summary_function=NULL){

  
  results <- parallel::mclapply(
                         1:length(paramsets),
                         function (i){
                           x <- paramsets[[i]]
                           if(! "name" %in% names(x)) { x$name <- i}
                           run_params(x, L,
                                      N_particles=N_particles,
                                      N_threads=N_threads_internal,
                                      run_params=x$run_params,
                                      run_name=x$name,
                                      silent=silent,
                                      return_summary_function=return_summary_function)
                         }, mc.cores=N_threads_external)

  results <- tryCatch({
    rbindlist(results)
  },error=function(cont){
    print(results)
    return(list())
  }
  )

  if(!"name" %in% paramsets[[1]]) results[, sim:=paste(name,sim)]
  
  return(results)                            
  }
