dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt
N_steps <- user()

# Vaccination
vax_time_step[,] <- vaccinations[step, i,j] #interpolate(interpolation_time, vaccinated_cont)
dim(vax_time_step) <- c(n, n_vac)
dim(vaccinations) <- c(N_steps, n, n_vac)
vaccinations[,,] <- user()


N_imp[,,] <- import_vec[step,i,j,k] #interpolate(interpolation_time, N_imp_input)
dim(import_vec) <- c(N_steps, n, n_vac, n_strain)
dim(N_imp) <- c( n, n_vac, n_strain)
import_vec[,,,] <- user()
                                        # import

# Daily varying beta
beta[] <- beta_day[step,i]#*exp(log_beta) #interpolate(interpolation_time, beta_day)
dim(beta) <- c(n)
dim(beta_day) <- c(N_steps, n)
beta_day[,] <- user()


                                        #Varrying severity
#severity[] <- severity_day[step,i]#interpolate(interpolation_time, severity_day)
#dim(severity) <- n
#dim(severity_day) <- c(length(interpolation_time), n)
#severity_day[,] <- user()

          
#vac_time_full_effect[] <- user()
#interpolate(interpolation_time, vac_time_full_effect_day)
#dim(vac_time_full_effect) <- n
#dim(vac_time_full_effect_day) <- c(length(interpolation_time),n)
#vac_time_full_effect_day[,] <- user()


#ng <- n/3
#ng2 <- ng*2


#S_waning[] <- 0
#S_waning[1:ng] <- S[i+ng] + S[i+ng2]
#S_waning[(ng+1):(ng2)] <- -S[i]
#S_waning[(ng2+1):(n)] <- -S[i]
#output(S_waning) <- TRUE
#sum_S_waning <- sum(S_waning[])
#output(sum_S_waning) <- TRUE

#dim(S_waning) <- n
## Core equations for transitions between compartments:
update(S[,]) <-  S[i,j] - sum(n_SE[i,j,]) +  vax_time_step[i,j] - n_waning[i,j] + sum(n_RS[i,j,]) -sum(N_imp[i,j,])# +n_S_2[i]  #+ (sum(mig_S[i,])- S[i]/reg_pop_long[i] * sum(migration_matrix[1:n,i]))*dt + dt*S_waning[i]/waning_immunity_vax[i]*0   - n_imp[i] 
#S[] <- if(S[i] <0 ) 0 S[i]

update(Ea[,,]) <- Ea[i,j,k] + n_SEa[i,j,k] - n_EaA[i,j,k] + n_RIA[i,j,k]# +  dt*(sum(mig_Ea[i,])- Ea[i]/reg_pop_long[i] * sum(migration_matrix[1:n,i]))

update(Es[,,]) <- Es[i,j,k] + n_SEi[i,j,k] - n_EsI[i,j,k] + n_RIS[i,j,k]# +   dt*( sum(mig_Es[i,])- Es[i]/reg_pop_long[i] * sum(migration_matrix[1:n,i]))

update(P[,,]) <-  P[i,j,k] + n_EsI[i,j,k] - n_PI[i,j,k]# + dt*(sum(mig_P[i,])- P[i]/reg_pop_long[i] * sum(migration_matrix[1:n,i]))
update(I[,,]) <-  I[i,j,k] + n_PI[i,j,k] - n_I[i,j,k] + N_imp[i,j,k]# + dt*(sum(mig_I[i,])- I[i]/reg_pop_long[i] * sum(migration_matrix[1:n,i]))

#update(I_imp[]) <-  I_imp[i] -n_I_imp[i] +  n_imp[i]

update(A[,,]) <-  A[i,j,k] + n_EaA[i,j,k] - n_AR[i,j,k]# + dt*(sum(mig_A[i,])- A[i]/reg_pop_long[i] * sum(migration_matrix[1:n,i]))
update(H[,,]) <-  H[i,j,k] + n_IH[i,j,k]  - n_H[i,j,k]
update(ICU_H[,,]) <- ICU_H[i,j,k] + n_IICU[i,j,k] - n_ICU_HR[i,j,k]
update(ICU_R[,,]) <- ICU_R[i,j,k] + n_ICU_HR[i,j,k] - n_ICU_RP[i,j,k]

update(ICU_P[,,]) <- ICU_P[i,j,k] -n_ICU_P[i,j,k] +n_ICU_RP[i,j,k]

update(B_D[,,]) <- B_D[i,j,k] + n_ID[i,j,k] - n_B_D_D[i,j,k] 

update(B_D_H[,,]) <- B_D_H[i,j,k] + n_HD[i,j,k] - n_B_H_D[i,j,k]

update(B_D_ICU[,,]) <- B_D_ICU[i,j,k] + n_ICU_D[i,j,k] - n_B_ICU_D[i,j,k]

#update(PRE_MISC[]) <- PRE_MISC[i] + n_I_PRE_MISC[i] - n_PRE_MISC[i]
#update(MISC[]) <- MISC[i] + n_PRE_MISC_H[i] - n_MISCR[i]

#update(MISC_ICU_H[]) <- MISC_ICU_H[i] + n_PRE_MISC_ICU_H[i] - n_MISC_ICU[i]
#update(MISC_ICU[]) <- MISC_ICU[i] + n_MISC_ICU[i] - n_MISC_ICU_R[i]

update(R[,,]) <-  R[i,j,k] + n_AR[i,j,k] + n_IR[i,j,k]  + n_ICU_R[i,j,k] + n_HR[i,j,k] - n_RI[i,j,k]- n_RS[i,j,k] #+ dt*(sum(mig_R[i,])- R[i]/reg_pop_long[i] * sum(migration_matrix[1:n,i])) + n_MISCR[i]- n_I_PRE_MISC[i] + n_I_imp[i]

update(D[,,]) <- D[i,j,k] + n_B_D_D[i,j,k] + n_B_H_D[i,j,k] + n_B_ICU_D[i,j,k]

#fraction_of_vaccines_given_to_S[1:ng] <- S[i] /(S[i] + Ea[i]+ Es[i] + P[i] + I[i] + A[i] + H[i] + ICU_H[i] + ICU_R[i]+ICU_P[i]+ R[i] + B_D[i] + B_D_H[i] + B_D_ICU[i] + 1e-10)
#fraction_of_vaccines_given_to_S[(ng+1):(ng2)] <- fraction_of_vaccines_given_to_S[i - ng]
#fraction_of_vaccines_given_to_S[(ng2+1):(n)] <- fraction_of_vaccines_given_to_S[i - ng2]
#updated(W[,]) <- W[i,j] + n_waning[i,j] - sum(n_WE[i,j,])


update(tot_infected[,,]) <-  tot_infected[i,j,k] + n_SE[i,j,k]
  

update(tot_hosp[,,]) <-  tot_hosp[i,j,k] + n_IH[i,j,k] + n_IICU[i,j,k]

update(tot_resp[,,]) <- tot_resp[i,j,k] + n_ICU_HR[i,j,k]
#update(tot_vac[]) <- tot_vac[i] + as.integer(round(N_vac[i]))
#update(tot_misc[]) <- tot_misc[i] + n_PRE_MISC[i]
#update(tot_imp[]) <- n_imp[i]
#update(hosp_inc) <- sum(n_IH[]) + sum(n_IICU[])
#update(log_beta) <- log_beta + rnorm(0, 0.05)


                                        #Transition probs:
                                        #S - > Ea , Ei
symp_asymp_effect[,,] <- (asympt_frac[i,j,k]*susceptibility_asymp[i,j,k] + sympt_frac[i,j,k]*susceptibility_symp[i,j,k])
p_SE[,,] <-1- exp(-sum(lambda_ij[i,j,,,k])*symp_asymp_effect[i,j,k]*dt)

p_EI <- 1 - exp(-dt/latent_period)
p_PI <- 1 - exp(-dt/pre_sympt_period)
p_IR <- 1 - exp(-dt/infectious_period)
p_HR[,,] <- 1 - exp(-dt/length_hosp[i,j,k])
p_ICU_R[,,] <- 1 - exp(-dt/pre_icu[i,j,k])
p_ICU_P[,,] <- 1 - exp(-dt/length_icu[i,j,k])
p_ICU_PR[,,] <- 1 - exp(-dt/post_icu[i,j,k])
p_waning[,] <- 1 - exp(-dt/T_waning[i,j])
p_waning_inf <- 1 - exp(-dt/waning_inf)


#p_MISCR[] <- 1 - exp(-dt/misc_length[i])

p_D <- 1 - exp(-dt/time_before_death)
p_D_H <- 1 - exp(-dt/time_before_death_hosp)
p_D_ICU <- 1 - exp(-dt/time_before_death_icu)

#p_full_effect[] <- 1 - exp(-dt/vac_time_full_effect[i,j,k])
#p_waning <- 1 - exp(-dt/waning_immunity)
dim(p_SE) <- c(n,n_vac,n_strain)
dim(p_HR) <- c(n,n_vac,n_strain)
dim(p_ICU_R) <- c(n,n_vac,n_strain)
dim(p_ICU_P) <- c(n,n_vac,n_strain)
dim(p_ICU_PR) <- c(n,n_vac,n_strain)
dim(symp_asymp_effect) <- c(n,n_vac,n_strain)
dim(n_RIA) <- c(n,n_vac,n_strain)
dim(n_RIS) <- c(n,n_vac,n_strain)
dim(n_RI) <- c(n,n_vac,n_strain)
dim(n_RI_op) <- c(n,n_vac,n_strain)
dim(p_waning) <- c(n, n_vac)
#dim(p_full_effect) <- c(n,n_vac,n_strain)
#dim(p_MISCR) <- (n,n_vac,n_strain)

## HARD CODED 2 strains
n_SE_tot[,] <- rbinom(S[i,j], sum(p_SE[i,j,]))
rel_strain[,,] <- p_SE[i,j,k]/sum(p_SE[i,j,])
n_SE[,,] <- if(k==1 || n_strain==1) rbinom(n_SE_tot[i,j],rel_strain[i,j,k]) else
              (if (k==2) n_SE_tot[i,j] - n_SE[i,j,1] else 0)

n_RI[,,] <- if( n_strain==1) 0 else
             (if (k==1) rbinom(R[i,j,1], 1 - exp(-sum(lambda_ij[i,j,,,1])*cross_protection[1,2]*symp_asymp_effect[i,j,k]* dt)) else (if (k==2)  rbinom(R[i,j,2], 1 - exp(-sum(lambda_ij[i,j,,,1])*cross_protection[2,1]*symp_asymp_effect[i,j,k] * dt)) else 0))
n_RI_op[,,] <- if(n_strain==1) 0 else
                                   ( if(k==1) n_RI[i,j,2] else n_RI[i,j,1])


n_RIA[,,] <- rbinom(n_RI_op[i,j,k], pa[i,j,k])
n_RIS[,,] <- n_RI_op[i,j,k] - n_RIA[i,j,k]



n_IR[,,] <- rbinom(I[i,j,k], p_IR)


n_SEa[,,] <- rbinom(n_SE[i,j,k], pa[i,j,k])



dim(n_SE_tot) <- c(n,n_vac)
dim(rel_strain) <- c(n,n_vac,n_strain)
n_waning[,] <- if(j != n_vac) rbinom(S[i,j], p_waning[i,j]) else(
                                                              -sum(n_waning[i, 1:(j-1)]))

#tot_waned[] <- sum(n_waning[i,])

#n_waning[,n_vac] <- -tot_waned[i]
dim(n_waning) <- c(n,n_vac)
#dim(tot_waned) <- n

#n_S_2[] <- rbinom(S[i] - n_SE[i]+n_vac[i], p_full_effect[i])
#n_S_2[1:ng] <- 0
#n_S_2[(ng+1):(ng2)] <- - n_S_2[i]
#n_S_2[(ng2+1):(n)] <- -n_S_2[i-ng]

#n_va[] <- as.integer(round((N_vac[i]*fraction_of_vaccines_given_to_S[i]))*dt)

#n_vac[1:ng] <- if (-n_va[i] < S[i]-n_SE[i]) n_va[i] else -(S[i]-n_SE[i])+1
#n_vac[(ng+1):(ng2)] <- -n_va[i-ng]
#n_vac[(ng2+1):(n)] <- 0


#dim(Simp) <- n
#Simp[] <- (S[i]-n_SE[i]+n_S_2[i]+ n_vac[i])*import_age_prio[i]
#S_I <- sum(Simp)
#n_imp[] <- rbinom(S[i] - n_SE[i]+ n_S_2[i]+n_vac[i], N_imp/S_I*import_age_prio[i])

  
#output(n_S_2) <- TRUE

pa[,,] <- asympt_frac[i,j,k]*susceptibility_asymp[i,j,k]/(asympt_frac[i,j,k]*susceptibility_asymp[i,j,k] + sympt_frac[i,j,k]*susceptibility_symp[i,j,k])

dim(pa) <- c(n, n_vac, n_strain)
n_SEi[,,] <- n_SE[i,j,k] - n_SEa[i,j,k]

n_EaA[,,] <- rbinom(Ea[i,j,k], p_EI)

n_EsI[,,] <- rbinom(Es[i,j,k], p_EI)
n_PI[,,] <- rbinom(P[i,j,k], p_PI)

n_AR[,,] <- rbinom(A[i,j,k], p_IR)
n_I[,,] <- rbinom(I[i,j,k], p_IR)
#n_I_imp[] <- rbinom(I_imp[i], p_IR)


# nested binomial
n_IR[,,] <- rbinom(n_I[i,j,k], p_hosp_icu_R[i,j,k])
p_hosp_icu_R[,,] <- 1 - hosp_prob[i,j,k] - (1-hosp_prob[i,j,k])*prob_death_non_hosp[i,j,k]

n_IH[,,] <- rbinom(n_I[i,j,k] - n_IR[i,j,k], p_hosp[i,j,k])
p_hosp[,,] <- hosp_prob[i,j,k]*(1-icu_prob[i,j,k])/(hosp_prob[i,j,k] + (1-hosp_prob[i,j,k])*prob_death_non_hosp[i,j,k])

n_IICU[,,] <- rbinom(n_I[i,j,k] - n_IR[i,j,k] - n_IH[i,j,k], p_icu[i,j,k])
p_icu[,,] <- hosp_prob[i,j,k]*icu_prob[i,j,k]/(hosp_prob[i,j,k]*icu_prob[i,j,k] + (1-hosp_prob[i,j,k])*prob_death_non_hosp[i,j,k])

n_ID[,,] <- n_I[i,j,k] - n_IR[i,j,k] - n_IH[i,j,k] - n_IICU[i,j,k]


## n_I_PRE_MISC[] <- rbinom(n_IR[i,j,k], prob_misc[i,j,k])
## n_PRE_MISC[] <- rbinom(PRE_MISC[i,j,k], 1 - exp(-dt/time_to_misc))
## n_PRE_MISC_H[] <- rbinom(n_PRE_MISC[i,j,k], 1 - misc_icu[i,j,k])
## n_PRE_MISC_ICU_H[] <- n_PRE_MISC[i] - n_PRE_MISC_H[i]
## n_MISC_ICU[] <- rbinom(MISC_ICU_H[i], 1 - exp(-dt/misc_icu_pre))
## n_MISC_ICU_R[] <- rbinom(MISC_ICU[i], 1 - exp(-dt/misc_icu_length))
## n_MISCR[] <- rbinom(MISC[i], p_MISCR[i])

dim(p_hosp_icu_R)<- c(n, n_vac, n_strain)
dim(p_hosp)<- c(n, n_vac, n_strain)
dim(p_icu)<- c(n, n_vac, n_strain)



n_H[,,] <- rbinom(H[i,j,k], p_HR[i,j,k])

n_HD[,,] <- rbinom(n_H[i,j,k], p_hosp_d[i,j,k])
p_hosp_d[,,] <- prob_death_hosp[i,j,k]
dim(p_hosp_d)<- c(n, n_vac, n_strain)
n_HR[,,] <- n_H[i,j,k] - n_HD[i,j,k]

n_ICU_HR[,,] <- rbinom(ICU_H[i,j,k], p_ICU_R[i,j,k])
n_ICU_RP[,,] <- rbinom(ICU_R[i,j,k], p_ICU_P[i,j,k])
n_ICU_P[,,] <- rbinom(ICU_P[i,j,k], p_ICU_PR[i,j,k])

n_ICU_D[,,] <- rbinom(n_ICU_P[i,j,k], p_icu_d[i,j,k])
p_icu_d[,,] <- prob_death_icu[i,j,k]
n_ICU_R[,,] <- n_ICU_P[i,j,k] - n_ICU_D[i,j,k]

dim(p_icu_d)<- c(n, n_vac, n_strain)

n_B_D_D[,,] <- rbinom(B_D[i,j,k], p_D)
n_B_H_D[,,] <- rbinom(B_D_H[i,j,k], p_D_H)
n_B_ICU_D[,,] <- rbinom(B_D_ICU[i,j,k], p_D_ICU)

n_RS[,,] <- rbinom(R[i,j,k], p_waning_inf)


dim(n_SE)<- c(n, n_vac, n_strain)
#dim(n_S_2)<- c(n, n_vac, n_strain)
#dim(n_imp)<- c(n, n_vac, n_strain)
#dim(n_vac)<- c(n, n_vac, n_strain)
#dim(n_va)<- c(n, n_vac, n_strain)
dim(n_SEi)<- c(n, n_vac, n_strain)
dim(n_SEa)<- c(n, n_vac, n_strain)
dim(n_EsI)<- c(n, n_vac, n_strain)
dim(n_EaA)<- c(n, n_vac, n_strain)
dim(n_PI)<- c(n, n_vac, n_strain)
dim(n_AR)<- c(n, n_vac, n_strain)
dim(n_I)<- c(n, n_vac, n_strain)
#dim(n_I_imp)<- c(n, n_vac, n_strain)
dim(n_IH)<- c(n, n_vac, n_strain)
dim(n_IICU)<- c(n, n_vac, n_strain)
dim(n_IR)<- c(n, n_vac, n_strain)
dim(n_ID)<- c(n, n_vac, n_strain)
dim(n_H)<- c(n, n_vac, n_strain)
dim(n_HD)<- c(n, n_vac, n_strain)
dim(n_HR)<- c(n, n_vac, n_strain)
dim(n_ICU_RP)<- c(n, n_vac, n_strain)
dim(n_ICU_P)<- c(n, n_vac, n_strain)
dim(n_ICU_HR)<- c(n, n_vac, n_strain)
dim(n_ICU_D)<- c(n, n_vac, n_strain)
dim(n_ICU_R)<- c(n, n_vac, n_strain)
dim(n_B_D_D)<- c(n, n_vac, n_strain)
dim(n_B_H_D)<- c(n, n_vac, n_strain)
dim(n_B_ICU_D)<- c(n, n_vac, n_strain)
dim(n_RS)<- c(n, n_vac, n_strain)
#dim(n_I_PRE_MISC)<- c(n, n_vac, n_strain)
#dim(n_PRE_MISC)<- c(n, n_vac, n_strain)
#dim(n_MISCR)<- c(n, n_vac, n_strain)
#dim(n_PRE_MISC_H)<- c(n, n_vac, n_strain)
#dim(n_PRE_MISC_ICU_H)<- c(n, n_vac, n_strain)
#dim(n_MISC_ICU)<- c(n, n_vac, n_strain)
#dim(n_MISC_ICU_R)<- c(n, n_vac, n_strain)



## Total population size (odin will recompute this at each timestep:
## automatically)
N[,] <- S[i,j]  + sum(Ea[i,j,])+ sum(Es[i,j,]) + sum(P[i,j,]) +  sum(I[i,j,]) + sum(A[i,j,]) + sum(H[i,j,]) + sum(ICU_H[i,j,]) + sum(ICU_R[i,j,]) + sum(ICU_P[i,j,]) + sum(B_D[i,j,]) + sum(B_D_ICU[i,j,]) + sum(B_D_H[i,j,]) + sum(R[i,j,]) + sum(D[i,j,])

update(tot_N) <- sum(N[,])

#output(N) <- TRUE
#output(tot_N) <- TRUE

                                        # Transitions

lambda_ij[,,,,] <- beta[i]*beta_strain[i5] * mixing_matrix[i,k]/beta_norm[i]*susceptibility[i,j,i5]*transmisibility[k,l,i5]*(pre_sympt_infect*P[k,l,i5] + symp_trans[k,l,i5]*(I[k,l,i5] ) + asympt_infect*A[k,l,i5]) #I_imp[j]
  
#output(lambda_ij) <- TRUE

## # Migrations
## mig_S[,] <- migration_matrix[i,j]*S[j]/reg_pop_long[j]
## mig_Ea[,] <- migration_matrix[i,j]*Ea[j]/reg_pop_long[j]
## mig_Es[,] <- migration_matrix[i,j]*Es[j]/reg_pop_long[j]
## mig_P[,] <- migration_matrix[i,j]*P[j]/reg_pop_long[j]
## mig_I[,] <- migration_matrix[i,j]*I[j]/reg_pop_long[j]
## mig_A[,] <- migration_matrix[i,j]*A[j]/reg_pop_long[j]
## mig_R[,] <- migration_matrix[i,j]*R[j]/reg_pop_long[j]

## Initial states:
initial(S[,]) <- S_ini[i,j] # will be user-defined
initial(Ea[,,]) <-  Ea_ini[i,j,k]
initial(Es[,,]) <-  Es_ini[i,j,k]
initial(P[,,]) <-  P_ini[i,j,k]
initial(I[,,]) <-  I_ini[i,j,k]
#initial(I_imp[,,]) <-  I_imp_ini[i,j,k]
initial(A[,,]) <-  A_ini[i,j,k]
initial(H[,,]) <-  H_ini[i,j,k]
initial(ICU_H[,,]) <-  ICU_H_ini[i,j,k]
initial(ICU_R[,,]) <- ICU_R_ini[i,j,k]
initial(ICU_P[,,]) <- ICU_P_ini[i,j,k]
initial(B_D[,,]) <- B_D_ini[i,j,k]
initial(B_D_H[,,]) <- B_D_H_ini[i,j,k]
initial(B_D_ICU[,,]) <- B_D_ICU_ini[i,j,k]
#initial(PRE_MISC[,,]) <- PRE_MISC_ini[i,j,k]
#initial(MISC[,,]) <- MISC_ini[i,j,k]
#initial(MISC_ICU_H[,,]) <- MISC_ICU_H_ini[i,j,k]
#initial(MISC_ICU[,,]) <- MISC_ICU_ini[i,j,k]
initial(R[,,]) <-  R_ini[i,j,k]
initial(D[,,]) <-  D_ini[i,j,k]
#initial(Ni[,,]) <-N[i,j,k] 
initial(tot_N) <- 0

initial(tot_infected[,,]) <- tot_infected_ini[i,j,k]
initial(tot_hosp[,,]) <- tot_hosp_ini[i,j,k]
initial(tot_resp[,,]) <- tot_resp_ini[i,j,k]
#initial(tot_vac[,,]) <- tot_vac_ini[i,j,k]
#initial(tot_misc[,,]) <- tot_misc_ini[i,j,k]
#initial(tot_imp[,,]) <- 0
#initial(log_beta) <- 0

## User defined parameters - default in parentheses:


dim(S)<- c(n, n_vac)
dim(Ea)<- c(n, n_vac, n_strain)
dim(Es)<- c(n, n_vac, n_strain)
dim(P)<- c(n, n_vac, n_strain)
dim(I)<- c(n, n_vac, n_strain)
#dim(I_imp)<- c(n, n_vac, n_strain)
dim(A)<- c(n, n_vac, n_strain)
dim(R)<- c(n, n_vac, n_strain)
dim(H)<- c(n, n_vac, n_strain)
dim(ICU_H)<- c(n, n_vac, n_strain)
dim(ICU_R)<- c(n, n_vac, n_strain)
dim(ICU_P)<- c(n, n_vac, n_strain)
dim(B_D)<- c(n, n_vac, n_strain)
dim(B_D_H)<- c(n, n_vac, n_strain)
dim(B_D_ICU)<- c(n, n_vac, n_strain)
#dim(MISC)<- c(n, n_vac, n_strain)
#dim(MISC_ICU)<- c(n, n_vac, n_strain)
#dim(MISC_ICU_H)<- c(n, n_vac, n_strain)
#dim(PRE_MISC)<- c(n, n_vac, n_strain)
dim(D)<- c(n, n_vac, n_strain)


dim(N)<- c(n, n_vac)
dim(tot_infected)<- c(n, n_vac, n_strain)
dim(tot_hosp)<- c(n, n_vac, n_strain)
dim(tot_resp)<- c(n, n_vac, n_strain)
#dim(tot_vac)<- c(n, n_vac, n_strain)
#dim(tot_misc)<- c(n, n_vac, n_strain)
#dim(tot_imp)<- c(n, n_vac, n_strain)
#dim(fraction_of_vaccines_given_to_S)<- c(n, n_vac, n_strain)


## dim(mig_S) <- c(n,n)
## dim(mig_Ea) <- c(n,n)
## dim(mig_Es) <- c(n,n)
## dim(mig_P) <- c(n,n)
## dim(mig_I) <- c(n,n)
## dim(mig_A) <- c(n,n)
## dim(mig_R) <- c(n,n)




#dim(N)<- c(n, n_vac, n_strain)
dim(S_ini)<- c(n, n_vac)
dim(Ea_ini)<- c(n, n_vac, n_strain)
dim(Es_ini)<- c(n, n_vac, n_strain)
dim(P_ini)<- c(n, n_vac, n_strain)
dim(I_ini)<- c(n, n_vac, n_strain)
dim(I_imp_ini)<- c(n, n_vac, n_strain)
dim(A_ini)<- c(n, n_vac, n_strain)
dim(R_ini)<- c(n, n_vac, n_strain)
dim(H_ini)<- c(n, n_vac, n_strain)
dim(ICU_H_ini)<- c(n, n_vac, n_strain)
dim(ICU_R_ini)<- c(n, n_vac, n_strain)
dim(ICU_P_ini)<- c(n, n_vac, n_strain)
dim(B_D_ini)<- c(n, n_vac, n_strain)
dim(B_D_H_ini)<- c(n, n_vac, n_strain)
dim(B_D_ICU_ini)<- c(n, n_vac, n_strain)
#dim(MISC_ini)<- c(n, n_vac, n_strain)
#dim(MISC_ICU_H_ini)<- c(n, n_vac, n_strain)
#dim(MISC_ICU_ini)<- c(n, n_vac, n_strain)
#dim(PRE_MISC_ini)<- c(n, n_vac, n_strain)
dim(D_ini)<- c(n, n_vac, n_strain)

dim(tot_infected_ini)<- c(n, n_vac, n_strain)
dim(tot_hosp_ini)<- c(n, n_vac, n_strain)
dim(tot_resp_ini)<- c(n, n_vac, n_strain)
#dim(tot_vac_ini)<- c(n, n_vac, n_strain)
#dim(tot_misc_ini)<- c(n, n_vac, n_strain)

#dim(reg_pop_long)<- c(n, n_vac, n_strain)
dim(beta_norm)<- n
dim(mixing_matrix) <- c(n, n)
dim(lambda_ij) <- c(n, n_vac, n, n_vac, n_strain)
#dim(migration_matrix) <- c(n,n)
dim(length_hosp)<- c(n, n_vac, n_strain)
dim(hosp_prob)<- c(n, n_vac, n_strain)
dim(length_icu)<- c(n, n_vac, n_strain)
dim(pre_icu)<- c(n, n_vac, n_strain)
dim(post_icu)<- c(n, n_vac, n_strain)
dim(icu_prob)<- c(n, n_vac, n_strain)
dim(prob_death_non_hosp)<- c(n, n_vac, n_strain)
dim(prob_death_hosp)<- c(n, n_vac, n_strain)
dim(prob_death_icu)<- c(n, n_vac, n_strain)
dim(susceptibility)<- c(n, n_vac, n_strain)
dim(transmisibility)<- c(n, n_vac, n_strain)
dim(waning_immunity_vax)<- c(n, n_vac, n_strain)
#dim(prob_misc)<- c(n, n_vac, n_strain)
#dim(misc_length)<- c(n, n_vac, n_strain)
dim(susceptibility_asymp)<- c(n, n_vac, n_strain)
dim(susceptibility_symp)<- c(n, n_vac, n_strain)
#dim(misc_icu)<- c(n, n_vac, n_strain)
dim(asympt_frac)<- c(n, n_vac, n_strain)
dim(sympt_frac)<- c(n, n_vac, n_strain)
dim(symp_trans)<- c(n, n_vac, n_strain)
dim(T_waning)<- c(n, n_vac)
dim(beta_strain) <- n_strain
mixing_matrix[,] <- user()
#migration_matrix[,] <- user()
latent_period <- user()
infectious_period <- user()
length_hosp[,,] <- user()
length_icu[,,] <- user()
hosp_prob[,,] <- user()
icu_prob[,,] <- user()
pre_icu[,,] <- user()
post_icu[,,] <- user()
#time_to_misc <- user()
#prob_misc[,,] <- user()


time_before_death <- user()
time_before_death_hosp <- user()
time_before_death_icu <- user()



susceptibility[,,] <- user()
transmisibility[,,] <- user()
susceptibility_asymp[,,] <- user()
susceptibility_symp[,,] <- user()
symp_trans[,,] <- user()

n_vac <- user()
n_strain <- user()
T_waning[,] <- user()
waning_inf <- user()
prob_death_non_hosp[,,] <- user()
prob_death_hosp[,,] <- user()
prob_death_icu[,,] <- user()
waning_immunity_vax[,,] <- user()
waning_immunity <- user(10000)

sympt_frac[,,] <- user()
asympt_frac[,,] <- user()
pre_sympt_infect <- user()
asympt_infect <- user()

#vac_time <- user()
#vac_time_full_effect_1 <- user()


pre_sympt_period <- user()


dim(cross_protection) <- c(n_strain, n_strain)
cross_protection[,] <- user()
n <- user(198)
S_ini[,] <- user()
I_ini[,,] <- user()
I_imp_ini[,,] <- user(0)
Ea_ini[,,] <- user(0)
Es_ini[,,] <- user(0)
A_ini[,,] <- user(0)
P_ini[,,] <- user(0)
beta_strain[] <- user()
H_ini[,,] <- user(0)
ICU_H_ini[,,] <- user(0)
ICU_R_ini[,,] <- user(0)
ICU_P_ini[,,] <- user(0)
B_D_ini[,,] <- user(0)
B_D_H_ini[,,] <- user(0)
B_D_ICU_ini[,,] <- user(0)
R_ini[,,] <- user(0)
D_ini[,,] <- user(0)
tot_infected_ini[,,] <- user(0)
tot_hosp_ini[,,] <- user(0)
tot_resp_ini[,,] <- user(0)
#tot_vac_ini[,,] <- user(0)
#reg_pop_long[,,] <- user()
beta_norm[] <- user(0)
#MISC_ini[,,] <- user(0)
#MISC_ICU_ini[,,] <- user(0)
#MISC_ICU_H_ini[,,] <- user(0)
#PRE_MISC_ini[,,] <- user(0)
#tot_misc_ini[,,] <- user(0)

