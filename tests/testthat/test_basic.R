library(data.table)
basic_params <- function(N=5, n_vac=2, L=100, n_strain=1){
  list(
    N_steps=L,
    n_vac=n_vac,
    n_strain=n_strain,
    dt=1,
    T_waning=array(1e10, dim=c(N, n_vac)),
    vaccinations=array(0,dim=c(L, N, n_vac)),
                                        #    import_vec=rep(0, L),
                                        #    import_age_prio=rep(1,dim=(N,n_vac,1)),
    beta_day=matrix(0.05, ncol=N, nrow=L),
                                        #    vac_time_full_effect=array(14.0, N),
    mixing_matrix=matrix(1.0, nrow=N, ncol=N),
    migration_matrix=matrix(0.0, nrow=N, ncol=N),
    latent_period=3.0,
    beta_strain=rep(1, n_strain),
    cross_protection=matrix(1, ncol=n_strain, nrow=n_strain),
    infectious_period=4.0,
    import_vec=array(0,dim=c(L,N,n_vac,n_strain)),
    length_hosp=array(4,dim=c(N,n_vac,n_strain)),
    length_icu=array(15,dim=c(N,n_vac,n_strain)),
    hosp_prob=array(0.1, dim=c(N,n_vac,n_strain)),
    icu_prob=array(0.5, dim=c(N,n_vac,n_strain)),
    pre_icu=array(2,dim=c(N,n_vac,n_strain)),
    post_icu=array(2,dim=c(N,n_vac,n_strain)),
    time_before_death=4,
    time_before_death_hosp=3,
    time_before_death_icu=2,
    susceptibility=array(1,dim=c(N,n_vac,n_strain)),
    transmisibility=array(1,dim=c(N,n_vac,n_strain)),
    susceptibility_asymp=array(1,dim=c(N,n_vac,n_strain)),
    susceptibility_symp=array(1,dim=c(N,n_vac,n_strain)),
    symp_trans=array(1.1,dim=c(N,n_vac,n_strain)),
    prob_death_non_hosp=array(0.03, dim=c(N,n_vac,n_strain)),
    prob_death_hosp=array(0.1, dim=c(N,n_vac,n_strain)),
    prob_death_icu=array(0.5,dim=c(N,n_vac,n_strain)),
    waning_immunity_vax = array(1000, dim=c(N,n_vac,n_strain)),
    sympt_frac=array(0.9,dim=c(N,n_vac,n_strain)),
    asympt_frac=array(0.1,dim=c(N,n_vac,n_strain)),
    pre_sympt_infect=1.5,
    asympt_infect=0.7,
    pre_sympt_period=2,
    n=N,
    S_ini=matrix(1e5, nrow=N, ncol=n_vac),
    I_ini=array(20, dim=c(N,n_vac,n_strain)),
    I_imp_ini=array(0, dim=c(N,n_vac,n_strain)),
    Ea_ini=array(0, dim=c(N,n_vac,n_strain)),
    Es_ini=array(0, dim=c(N,n_vac,n_strain)),
    A_ini=array(0, dim=c(N,n_vac,n_strain)),
    P_ini=array(0, dim=c(N,n_vac,n_strain)),
    H_ini=array(0, dim=c(N,n_vac,n_strain)),
    ICU_H_ini=array(0, dim=c(N,n_vac,n_strain)),

    ICU_R_ini=array(0, dim=c(N,n_vac,n_strain)),
    ICU_P_ini=array(0, dim=c(N,n_vac,n_strain)),
    B_D_ini=array(0, dim=c(N,n_vac,n_strain)),
    B_D_H_ini=array(0, dim=c(N,n_vac,n_strain)),
    B_D_ICU_ini=array(0, dim=c(N,n_vac,n_strain)),
    R_ini=array(0, dim=c(N,n_vac,n_strain)),
    D_ini=array(0, dim=c(N,n_vac,n_strain)),
#    MISC_ini=array(0, dim=c(N,n_vac,n_strain)),
#    MISC_ICU_ini=array(0, dim=c(N,n_vac,n_strain)),
#    MISC_ICU_H_ini=array(0, dim=c(N,n_vac,n_strain)),
#    PRE_MISC_ini=array(0, dim=c(N,n_vac,n_strain)),
#    tot_misc_ini=array(0, dim=c(N,n_vac,n_strain)),
    tot_infected_ini=array(0, dim=c(N,n_vac,n_strain)),
    tot_hosp_ini=array(0, dim=c(N,n_vac,n_strain)),
    tot_resp_ini=array(0, dim=c(N,n_vac,n_strain)),
 #   tot_vac_ini=array(0, dim=c(N,n_vac,n_strain)),
    beta_norm=rep(1e5,N),
    reg_pop_long=rep(1e5,N),
    waning_inf=1e10,
    N_regions=1,
    age_groups=N

  )
}

test_that("Conserve N", {
  results <- run_params(basic_params(), L=100, 3, 3)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(basic_params()$S_ini) + sum(basic_params()$I_ini))
  expect_true(N_t)
})


test_that("Test vaccinaion implemenation", {
  results_no_vax <- run_params(basic_params(), L=100, 3, 3)
  params <- basic_params()
  vax <- array(0,dim=c(100, 5, 2))
  vax[5,,1] <- rep(-50000, params$n)
  vax[5,,2] <- rep(50000,  params$n)

  params$vaccinations <- vax
  params$susceptibility[,2,] <- 0.5
  results <- run_params(params, L=100, 3, 3)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(basic_params()$S_ini) + sum(basic_params()$I_ini))
  expect_true(N_t)

  expect_lte(mean(results[time==100,tot_infected]), mean(results_no_vax[time==100,tot_infected]))
  expect_true(mean(results[time==100,get("tot_infected[1]")])/50000 > mean(results[time==100,get("tot_infected[6]")/150000]))
  
})


test_that("Test Waning vaccine", {
  params <- basic_params(n_vac=3)
  params$T_waning[,2] <- 3
  params$beta_day[,] <- 0

  results <- run_params(params, L=100, 3, 3)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$S_ini) + sum(params$I_ini))
  expect_true(N_t)
  expect_gte(mean(results[time==100,get("S[11]")]), 190000)
  
})

test_that("Test import", {
  params <- basic_params(n_vac=1, n_strain=2, N=5)
  params$I_ini[,,] <- 0
  params$import_vec[30:100,1,1,1] <- 10
  params$import_vec[40:100,1,1,2] <- 10
  results <- run_params(params, L=100, 3, 3)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$S_ini) + sum(params$I_ini))
  expect_true(N_t)

  expect_lte(max(results[time<30,get("I[1]")]), 1)
  expect_gte(max(results[time>30,get("I[1]")]), 1)
  expect_lte(max(results[time<40,get("I[6]")]), 1)
  expect_gte(max(results[time>40,get("I[6]")]), 1)
  
})


test_that("Test Waning infection", {
  params <- basic_params(n_vac=1, L=300)
  params$waning_inf <- 100

  results <- run_params(params, L=300, 3, 3)
  # Crude test for two peaks
  expect_gte(  max(results[t <150, get("I[5]")]), 5000)
  expect_lte(  max(results[t >150 & t < 200, get("I[5]")]), 1000)
  expect_gte(  max(results[t >200, get("I[5]")]), 2000)
  
})



test_that("Test 2 strains", {
  params <- basic_params(n_vac=1, n_strain=2)
  params$beta_strain <- c(1,2)

  results <- run_params(params, L=100, 3, 3)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$S_ini) + sum(params$I_ini))
  expect_true(N_t)
  expect_gte(mean(results[time==100,get("tot_infected[6]")]), mean(results[time==100,get("tot_infected[1]")]))
  
})



## test_that("Basic model FULL run with N=1", {

##   for(i in 1:5){

##     pars <- list(
##       n_strain=2,
##       n_vac=2,
##       n=1,
##       S_ini=matrix(100000,nrow=1, ncol=2),
##       I_ini=array(c(10000,10000,100,100), dim=c(1,2,2)),
##       R_ini=array(c(200000,0,0,0), dim=c(1,2,2)),
##       beta_strain=c(1,4),
##       mixing_matrix=matrix(1, nrow=1, ncol=1),
##       susceptibility=array(c(1,1.0, 1, 1.0), dim=c(1,2,2)),
##       transmisibility=array(1, dim=c(1,2,2)),
##       beta_day=rep(0.4, 2000),
##       dt=0.1,
##       gamma=1/5,
##       N_steps=2000,
##       beta_norm=c(10000),
##       cros_protection=matrix(c(0,0,0,0), nrow=2, ncol=2)
##     )
    
##     mod <- metapopfull$new(pars=pars,
##                            step=1,
##                            n_particles=500,
##                            n_threads=1)
##     res <- mod$simulate(1:500)
##   }
##   expect_equal(1,1)
## })
