library(data.table)
library(dplyr)
basic_params <- function(N=9, n_vac=2, L=100, n_strain=1){
  list(
    N_steps=L,
    n_vac=n_vac,
    n_strain=n_strain,
    dt=1,## 
    T_waning=array(1e10, dim=c(N, n_vac)),
    vaccinations=array(0,dim=c(L, N, n_vac)),
                                        #    import_vec=rep(0, L),
                                        #    import_age_prio=rep(1,dim=(N,n_vac,1)),
    beta_day=matrix(0.05, ncol=N, nrow=L),
                                        #    vac_time_full_effect=array(14.0, N),
    mixing_matrix=matrix(1.0, nrow=N, ncol=N),
    tot_vac_ini=array(0, dim=c(N,n_vac)),
    migration_matrix=matrix(0.0, nrow=N, ncol=N),
    latent_period=3.0,
    beta_strain=rep(1, n_strain),
    cross_protection=matrix(1, ncol=n_strain, nrow=n_strain),
    infectious_period=4.0,
    import_vec=array(0,dim=c(L,N,n_vac,n_strain)),
    length_hosp=array(4,dim=c(N,n_vac,n_strain)),
    length_icu=array(15,dim=c(N,n_vac,n_strain)),
    hosp_prob=array(0.04, dim=c(N,n_vac,n_strain)),
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
    symp_trans=array(1.0,dim=c(N,n_vac,n_strain)),
    prob_death_non_hosp=array(0.03, dim=c(N,n_vac,n_strain)),
    prob_death_hosp=array(0.1, dim=c(N,n_vac,n_strain)),
    prob_death_icu=array(0.5,dim=c(N,n_vac,n_strain)),
    waning_immunity_vax = array(10000, dim=c(N,n_vac,n_strain)),
    sympt_frac=array(0.7,dim=c(N,n_vac,n_strain)),
    asympt_frac=array(0.3,dim=c(N,n_vac,n_strain)),
    pre_sympt_infect=1.2,
    asympt_infect=0.2,
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
    waning_inf=5e10,
    N_regions=1,
    rand_beta_sd=0.1,
    rand_beta_factors=rep(0.05,N),
    S_div_vac=array(0, dim=c(N, n_vac)),
    beta_cut_peak_param=c(0,0,0,0),
    age_groups=N,
    change_factor=c(0,0,0,0),
    threshold=c(100,0),
    beta_mode=1,
    spont_behav_change_params=c(0,0,0,0,0,0),
    expected_health_loss=array(1, dim=c(N,n_vac))
    

  )
}


test_that("Conserve N", {
  results <- run_params(basic_params(L=500), L=500, 30, 3)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(basic_params()$S_ini) + sum(basic_params()$I_ini))
  expect_true(N_t)
  
})

test_that("Both ways of calculating hosp inc are equivalaent", {
  pars <- basic_params(L=300, N=2)
  pars$beta_day <- pars$beta_day*2
  pars <- change_dt(pars, 1/8)
  results <- run_params(pars, L=100, 1, 1)
  tot_hosps <- results %>% filter(time %% 1 == 0 ) %>%pull(tot_hosp)
  est_incs <- tot_hosps[2:length(tot_hosps)] - tot_hosps[1:(length(tot_hosps)-1)]
  expect_equal(est_incs, results %>% filter(time %% 1 == 0 & time > 0) %>%pull(tot_hosp_inc))

})

test_that("Random beta", {
  params <- basic_params()
  params$beta_mode <- 3
  results1 <- run_params(params, L=100, 6, 3)
  N_t <- all(results1 %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(basic_params()$S_ini) + sum(basic_params()$I_ini))
  expect_true(N_t)
  params$rand_beta <- 0.4
  results2 <- run_params(params, L=100, 3, 3)
  tmp <- results1[, .(peak=max(get("I[1]"))), by=sim]
  tmp2 <- results2[, .(peak=max(get("I[1]"))), by=sim]
  expect_gte(max(tmp2) - min(tmp2), max(tmp) - max(tmp))
}
)
test_that("Test vaccinaion implemenation", {
  results_no_vax <- run_params(basic_params(N=2), L=100, 3, 3)
  params <- basic_params(N=2)
  vax <- array(0,dim=c(100, 2, 2))
  vax[5,,1] <- rep(-50000, params$n)
  vax[5,,2] <- rep(50000,  params$n)

  params$vaccinations <- vax
  params$susceptibility_symp[,2,] <- 0.5
  params$susceptibility_asymp[,2,] <- 0.5
  results <- run_params(params, L=100, 3, 3)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$S_ini) + sum(params$I_ini))
  expect_true(N_t)

  expect_lte(mean(results[t==100,tot_infected]), mean(results_no_vax[t==100,tot_infected]))
  expect_true(mean(results[t==100,get("tot_infected[1]")])/50000 > mean(results[t==100,get("tot_infected[2]")/150000]))
  
})

test_that("Test vaccinaion implemenation, vax_type=2", {
  pars <- basic_params(N=1, n_vac=3)
  pars$S_div_vac <- array(c(1,1,1), dim=c(1,3))
  pars$vax_type <- 2
  pars$beta_day <- pars$beta_day*0

  vax <- array(0,dim=c(100, 1, 3))
  vax[1:30,,1] <- 3e5/30
  pars$vaccinations <- vax
  results <- run_params(pars, L=100, 1, 1)

  results[1:35, "S[1]"]
  results[1:35, "S[2]"]
  results[1:35, "S[3]"]

  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(pars$S_ini) + sum(pars$I_ini))
  expect_true(N_t)

  expect_lte(results[t==31, get("S[1]")], 500)
  expect_lte(results[t==31, get("S[2]")], 500)
  expect_gte(results[t==31, get("S[3]")],3e5-500)

})

test_that("Test simple waning from R", {
  params <- basic_params(N=2, n_vac=1)
  params$R_ini[1:2,1,1] <- 1e5
  params$S_ini[1:2,1] <- 0
  params$waning_inf <- 5
  results <- run_params(params, L=100, 1, 1)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$R_ini) + sum(params$I_ini))
  expect_true(N_t)
  expect_gte(mean(results[t==100,get("S[1]")]), 90000)
  expect_lte(mean(results[t==100,get("R[1]")]), 10000)
})


test_that("Test complicated waning from R with strains", {
  params <- basic_params(N=1, n_vac=4, n_strain=2)
  params$R_ini[1,1:4,1] <- 1e5
  params$R_ini[1,1:4,2] <- 2e5
  params$S_ini[1,1:4] <- 0
  params$I_ini <- params$I_ini*1
  params$waning_inf <- 5
  params$beta_day <- params$beta_day*0
  params$vac_struct_length <- 2
  results <- run_params(params, L=100, 1, 1)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$R_ini) + sum(params$I_ini))
  expect_true(N_t)

  expect_equal(mean(results[t==100,get("S[1]")]), 0)
  expect_gte(mean(results[t==100,get("S[2]")]), 4e5)
  expect_equal(mean(results[t==100,get("S[3]")]), 0)
  expect_gte(mean(results[t==100,get("S[4]")]), 8e5)

  params$T_waning <- matrix(c(2e10, 5, 1e10, 5), nrow=1, ncol=4)
  params$include_waning <- 2
  results <- run_params(params, L=100, 1, 1)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$R_ini) + sum(params$I_ini))
  expect_true(N_t)

  expect_equal(mean(results[t==100,get("S[2]")]), 0)
  expect_gte(mean(results[t==100,get("S[1]")]), 4e5)
  expect_equal(mean(results[t==100,get("S[4]")]), 0)
  expect_gte(mean(results[t==100,get("S[3]")]), 8e5)
})


test_that("Test Waning vaccine", {
  params <- basic_params(N=1, n_vac=3)
  params$T_waning[,2] <- 3
  params$beta_day[,] <- 0
  params$include_waning <- 1

  results <- run_params(params, L=100, 3, 3)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$S_ini) + sum(params$I_ini))
  expect_true(N_t)
  expect_gte(mean(results[t==100,get("S[3]")]), 180000)
  
})


test_that("Test Waning vaccine type 2", {
  params <- basic_params(N=1, n_vac=3)
  params$T_waning[,] <- c(1e10,10,3)
  params$beta_day[,] <- 0
  params$include_waning <- 2

  results <- run_params(params, L=100, 3, 3)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$S_ini) + sum(params$I_ini))
  expect_true(N_t)
  expect_gte(mean(results[t==100,get("S[1]")]), 280000)
  expect_gte(mean(results[t==25,get("S[2]")]),mean(results[t==25,get("S[2]")]))
  
})

test_that("Test import", {
  params <- basic_params(n_vac=1, n_strain=2, N=9, L=500)
  params$I_ini[,,] <- 0
  params$import_vec[30:100,1,1,1] <- 10
  params$import_vec[40:100,1,1,2] <- 10
  results <- run_params(params, L=205, 3, 3)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$S_ini) + sum(params$I_ini))
  expect_true(N_t)

  expect_lte(max(results[t<30,get("I[1]")]), 1)
  expect_gte(max(results[t>30,get("I[1]")]), 1)
  expect_lte(max(results[t<40,get("I[10]")]), 1)
  expect_gte(max(results[t>40,get("I[10]")]), 1)
  
})


test_that("Test Waning infection", {
  params <- basic_params(N=9, n_vac=2, L=300)
  params$waning_inf <- 100
  params$I_ini[,2,] <- 0
  params$beta_day=matrix(0.03, ncol=9, nrow=300)
  beta_1 <- fix_beta_large(params, params$S_ini, params$I_ini, 1, beta=params$beta_day[1,], use_eig=TRUE)
  results <- run_params(params, L=300, 3, 3)
  # Crude test for two peaks
  expect_gte(  max(results[t <150, get("I[5]")]), 5000)
  expect_lte(  max(results[t >150 & t < 200, get("I[5]")]), 3000)
  expect_gte(  max(results[t >200, get("I[5]")]), 2000)
  
})

test_that("Test Rt calculation", {
  params <- basic_params(N=5, n_vac=2, L=500)
  params$I_ini[,,] <- 1000
  beta_1 <- fix_beta_large(params, params$S_ini, params$I_ini, 1, beta=params$beta_day[1,], use_eig=TRUE)
  params$beta_day=params$beta_day*beta_1
  params$dt <- 0.1
  results <- run_params(params, L=50, 1, 1)
  inc <- results[, mean(get("incidence")), by=t]

  expect_lte(abs(inc[t==200, V1] - inc[t==400, V1]), 100)

  # complicated contact matrix, susceptibility and transmissibility
  params <- basic_params(N=5, n_vac=2, L=168*10)
  params$S_ini <- array((1:5)*1e5, dim=c(5,2))
  params$S_ini[,2] <- 0
  params$contact_matrix <- matrix(1:25, ncol=5, nrow=5)
  params$susceptibility_asymp <- array(c(0.1, 0.4, 0.7, 0.9, 1.0), dim=c(5,2,1))
  params$susceptibility_symp <- array(c(0.5, 0.8, 0.8, 0.9, 0.5), dim=c(5,2,1))
  params$transmisibility <- array(c(0.4, 1.0, 0.5, 0, 0), dim=c(5,2,1))
  params$I_ini[,,] <- 200
  beta_1 <- fix_beta_large(params, params$S_ini, params$I_ini, 1, beta=params$beta_day[1,], use_eig=TRUE)
  params$beta_day <- params$beta_day*beta_1
  params$dt <- 0.1
  results <- run_params(params, L=168, 3, 3)
  inc <- results[, mean(get("incidence")), by=t]
  expect_lte(abs(inc[t==200, V1] - inc[t==400, V1]), 100)





  ## plot(results[, mean(get("I[1]")), by=t])
  ## plot(results[, mean(get("P[1]")), by=t])
  ## plot(results[, mean(incidence), by=t])
  ## plot(results[, mean(Rt), by=t])
  ## # Crude test for two peaks
  ## expect_gte(  max(results[t <150, get("I[5]")]), 5000)
  ## expect_lte(  max(results[t >150 & t < 200, get("I[5]")]), 3000)
  ## expect_gte(  max(results[t >200, get("I[5]")]), 2000)
  
})



test_that("Test 2 strains", {
  params <- basic_params(n_vac=1, n_strain=2, L=150)
  params$beta_strain <- c(1,2)

  results <- run_params(params, L=100, 3, 3)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$S_ini) + sum(params$I_ini))
  expect_true(N_t)
  expect_gte(mean(results[t==100,get("tot_infected_strain_2")]), mean(results[t==100,get("tot_infected_strain_1")]))
  
})

test_that("Test cut_peak", {
  params <- basic_params(n_vac=1, n_strain=1, L=200)

  beta_1 <- fix_beta_large(params, params$S_ini, params$I_ini, 1, beta=params$beta_day[1,], use_eig=TRUE)
  bv <- 0.05*beta_1*2
  params$beta_day <- params$beta_day*beta_1*2
  results1 <- run_params(params, L=200, 3, 3)
  params$beta_mode <- 5
  params$beta_cut_peak_param <- c(bv, bv*0.2, 0.5,14)
  params$rand_beta_factors <- rep(1,9)
  results2 <- run_params(params, L=200, 3, 3)
  
  params$beta_cut_peak_param <- c(bv, bv*0.2, 0.5,21)
  params$rand_beta_factors <- rep(1,9)
  results3 <- run_params(params, L=200, 3, 3)
  
  ## plot(results1[, mean(tot_infected), by=t])
  ## lines(results2[, mean(tot_infected), by=t], col="red")
  ## lines(results3[, mean(tot_infected), by=t], col="blue")
  
  N_t <- all(results2 %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$S_ini) + sum(params$I_ini))
  expect_true(N_t)
  
  expect_gte(mean(results1[t==200, get("tot_infected")]), mean(results2[t==200,get("tot_infected")])+100000)
  expect_gte(mean(results2[t==200, get("tot_infected")]), mean(results3[t==200,get("tot_infected")])+45000)
  
})

## test_that("Test cut_peak + dynamic change", {
##   params <- basic_params(n_vac=1, n_strain=1, L=200)
##   beta_1 <- fix_beta_large(params, params$S_ini, params$I_ini, 1, beta=params$beta_day[1,], use_eig=TRUE)
##   bv <- 0.05*beta_1*2

##   params$beta_dynamic_change <- 1
##   params$dyn_change <- c(bv, 0.5, 0.1,1)
##   params$rand_beta_factors <- rep(1,9)
##   results1 <- run_params(params, L=200, 3, 3)
##   params$beta_cut_peak <- 1
  
##   params$beta_cut_peak_param <- c(bv, bv*0.2, 1/2,14)
##   params$rand_beta_factors <- rep(1,9)
##   results2 <- run_params(params, L=200, 3, 3)
##   N_t <- all(results2 %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$S_ini) + sum(params$I_ini))
##   expect_true(N_t)
  
##   expect_gte(mean(results1[t==200,get("tot_infected")]), mean(results2[t==200,get("tot_infected")])+100000)
  
## })


## test_that("Test dynamic change", {
##   params <- basic_params(n_vac=1, n_strain=1, L=200)

##   results1 <- run_params(params, L=200, 3, 3)
##   params$beta_dynamic_change <- 1
##   params$dyn_change <- c(0.05, 0.5, 0.1,1)
##   params$rand_beta_factors <- rep(1,9)
##   results2 <- run_params(params, L=200, 3, 3)
##   N_t <- all(results2 %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$S_ini) + sum(params$I_ini))
##   expect_true(N_t)
  
##   expect_gte(max(results1[, "incidence"]), max(results2[,get("incidence")]))
  
## })

test_that("Test spontaneous behviour change", {
  params <- basic_params(n_vac=1, n_strain=1, L=200)
  params$beta_mode <- 4
  params$spont_behav_change_params <- c(0.05, 1, 2000, 4/(0.9*2000), 0.8, 0.7)
  params$expected_health_loss <- array(0, dim=c(9,1))
  results1 <- run_params(params, L=200, 3, 3)

  params$expected_health_loss <- array(50000, dim=c(9,1))
  results2 <- run_params(params, L=200, 3, 3)

  params$expected_health_loss <-  array(100000, dim=c(9,1))
  results3 <- run_params(params, L=200, 3, 3)

  N_t <- all(results2 %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$S_ini) + sum(params$I_ini))
  expect_true(N_t)
  
  expect_gte(max(results1[, "incidence"]), max(results2[,get("incidence")]))
  expect_gte(max(results2[, "incidence"]), max(results3[,get("incidence")]))
  
})



test_that("Test dynamic threshold", {
  params <- basic_params(L=600)
  params$beta_mode <- 2
  params$threshold <- c(1000,0)
  params$threshold_ini <- 0.04
  params$threshold_max <- 0.04
  params$threshold_min <- 0.04*0.1
  params$rand_beta_factors <- rep(1,9)
  params$change_factor <- c(0,2.5, 0.03, 30)
  results <- run_params(params, L=500, 3, 3, estimate_Rt=FALSE)
  N_t <- all(results %>% dplyr::filter(t!=1) %>% dplyr::pull(tot_N) == sum(params$S_ini) + sum(params$I_ini))
  expect_true(N_t)

  expect_gte(2200, max(results[, hosp]))
  
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
