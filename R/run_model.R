
#' run_params
#'
#' Run the metapopulation model with a set of parameters and a given number of particles and threads.
#' The number of particles gives the number of samples and threads the number of cores used to sample

run_params <- function(params, L=200, N_particles=1, N_threads=1, run_name="run1", run_params=list(), silent=TRUE, estimate_Rt=FALSE){
  if(!silent){
    print(glue::glue("Running {run_name}"))
  }
  
  dust_model <- model$new(pars = params,
                         step = 1,
                         n_particles = N_particles,
                         n_threads = N_threads,
                         )

  raw_results <- dust_model$simulate(1:L)
  params$dust_index <- dust_model$info()$index
  results <- refine_results_odin_dust(raw_results, params, N_threads)
  beta <- list()
  if(params$rand_beta == 1){
    for(n in 1:N_particles){
      beta[[length(beta) + 1]] <- outer(exp(results[sim==n, log_beta]), params$rand_beta_factors)
    }
  }
  else if(params$threshold_beta == 1){
    for(n in 1:N_particles){
      beta[[length(beta) + 1]] <- outer(results[sim==n, beta_thresh], params$rand_beta_factors)
    }
  }
  else{
      for(n in 1:N_particles){
        beta[[length(beta) + 1]] <- params$beta_day
      }
      
  }


  if(estimate_Rt){
    Rts <- parallel::mclapply(unique(results$sim), function(x) estimate_Rt_large_new(results%>%dplyr::filter(sim==x), params, params$beta_day), mc.cores=N_threads)
    results <- results %>% dplyr::mutate(Rt=unlist(Rts))
  }
  results <- results  %>% dplyr::mutate(name=run_name) %>% dplyr::mutate(!!!run_params)
  return(results)
}




run_param_sets <- function(paramsets, L=100, N_particles=1, N_threads_internal=1, N_threads_external=1, silent=TRUE){

  
  results <- parallel::mclapply(
                         paramsets,
                         function (x) run_params(x, L,
                                                 N_particles=N_particles,
                                                 N_threads=N_threads_internal,
                                                 run_params=x$run_params,
                                                 run_name=x$name,
                                                 silent=silent), mc.cores=N_threads_external)
  results <- tryCatch({
    rbindlist(results)
  },error=function(cont){
    print(results)
    return(list())
  }
  )

  return(results)                            
  }
