

fix_params <- function(params, N, n_vac, n_strain, vac_pars){

  params$asympt_frac <- array(params$asympt_frac, dim=c(N, n_vac, n_strain))
  params$sympt_frac <- array(params$sympt_frac, dim=c(N, n_vac, n_strain))
  params$length_hosp <- array(params$length_hosp, dim=c(N, n_vac, n_strain))
  params$length_icu <- array(params$length_icu, dim=c(N, n_vac, n_strain))
  params$hosp_prob <- array(outer(params$hosp_prob, vac_pars$rr_hosp), dim=c(N, n_vac, n_strain))
  params$icu_prob <- array(outer(params$icu_prob, vac_pars$rr_icu), dim=c(N, n_vac, n_strain))
  params$pre_icu <- array(params$pre_icu, dim=c(N, n_vac, n_strain))
  params$post_icu <- array(params$post_icu, dim=c(N, n_vac, n_strain))
  params$susceptibility <- array(outer(params$susceptibility, vac_pars$rr_inf), dim=c(N, n_vac, n_strain))
  params$susceptibility_symp <- array(vac_pars$rr_inf,, dim=c(N, n_vac, n_strain))
  params$susceptibility_asymp <- array(vac_pars$rr_inf_asymp, dim=c(N, n_vac, n_strain))

  params$transmisibility <- array(outer(params$transmisibility, vac_pars$rr_trans), dim=c(N, n_vac, n_strain))
  params$symp_trans <- array(1, dim=c(N, n_vac, n_strain))

  
  params$prob_death_non_hosp <- array(outer(params$prob_death_non_hosp, vac_pars$rr_death), dim=c(N, n_vac, n_strain))
  params$prob_death_hosp <- array(outer(params$prob_death_hosp, vac_pars$rr_death), dim=c(N, n_vac, n_strain))
  params$prob_death_icu <- array(outer(params$prob_death_icu, vac_pars$rr_death), dim=c(N, n_vac, n_strain))
  return(params)
}
                     
