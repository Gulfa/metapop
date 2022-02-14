
#' run_params
#'
#' Run the metapopulation model with a set of parameters and a given number of particles and threads.
#' The number of particles gives the number of samples and threads the number of cores used to sample

run_params <- function(params, L=200, N_particles=1, N_threads=1, run_name="run1", run_params=list()){

  dust_model <- model$new(pars = params,
                         step = 1,
                         n_particles = N_particles,
                         n_threads = N_threads,
                         )

  raw_results <- dust_model$simulate(1:L)
  params$dust_index <- dust_model$info()$index
  results <- refine_results_odin_dust(raw_results, params, N_threads)
  Rts <- parallel::mclapply(unique(results$sim), function(x) estimate_Rt_large_new(results%>%dplyr::filter(sim==x), params, params$beta_day), mc.cores=N_threads)
  Rts <-lapply(unique(results$sim), function(x) estimate_Rt_large_new(results%>%dplyr::filter(sim==x), params, params$beta_day))
  results <- results %>% dplyr::mutate(Rt=unlist(Rts), name=run_name) %>% dplyr::mutate(!!!run_params)
  return(results)
}




run_param_sets <- function(paramsets, L=100, N_particles=1, N_threads_internal=1, N_threads_external=1){
  
  results <- parallel::mclapply(paramsets, function (x) run_params(x, L, N_particles=N_particles,
                                                                   N_threads=N_threads_internal,
                                                                   run_params=x$run_params,
                                                                   run_name=x$name), mc.cores=N_threads_external)
  results <- tryCatch({
    rbindlist(results)
  },error=function(cont){
    print(results)
    return(list())
  }
  )

  return(results)                            
  }
