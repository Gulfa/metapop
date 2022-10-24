


#' @export
get_filter_tot_hosp <- function(incidence,params, n_particles=10, n_threads=10, data_time_interval="day" ){
  
  dust_model <- model$new(pars = params,
                          step = 1,
                          n_particles = 1,
                          n_threads = 1
                          )
  if(data_time_interval == "day"){
    filter_data <- mcstate::particle_filter_data(data = data.table(incidence=incidence, day=1:length(incidence)),
                                                 time = "day",
                                                 rate = 1/params$dt)
  }else{
    stop("Not implemented")
  }

  hosp_compare <- function(state, observed, pars = NULL) {
    exp_noise <- 1e6
    
    incidence_observed <- observed$incidence
    incidence_modelled <- state[1, ,drop=TRUE]
    lambda <- incidence_modelled +
      rexp(n = length(incidence_modelled), rate = exp_noise)
    dpois(x = incidence_observed, lambda = lambda, log = TRUE) #+ dist_ll
  }
  index_list <- list(run=c(dust_model$info()$index$tot_hosp_inc), state=unlist(dust_model$info()$index))
  
  index_func <- purrr::partial(function(x, index_list=index_list) index_list, index_list=!!index_list)
  
  filter <- mcstate::particle_filter$new(data = filter_data,
                                         model = model,
                                         n_particles = n_particles,
                                         n_threads=n_threads,
                                         index=index_func,
                                         compare = hosp_compare)
  return(filter)
}

