dt <- user(1)
steps_per_day <- 1/dt
initial(time) <- 0
update(time) <- (step + 1) * dt
N_steps <- user()

beta <- beta_day[step]
dim(beta_day) <- N_steps
beta_day[] <- user()
import <- user(0)

## i : metapop
## j : vaccine status
## k: variant

## S [i,j]
## I [i,j,k]
## R [i,j,k]


update(S[,]) <- S[i,j] - sum(n_SI[i,j,]) #- n_import[i]
update(I[,,]) <- I[i,j,k] + n_SI[i,j,k] - n_IR[i,j,k] + n_RI_op[i,j,k]#+ n_import[i]
update(R[,,]) <- R[i,j,k] + n_IR[i,j,k] - n_RI[i,j,k]
update(inc[,,]) <- if(step %% steps_per_day == 0) n_SI[i,j,k] else inc[i,j,k] + n_SI[i,j,k]
## Individual probabilities of transition:

p_SI[,,] <- 1 - exp(-sum(lambda_ij[i,j,,,k])* dt) # S to I
p_IR <- 1 - exp(-gamma * dt) # I to R

## Draws from binomial distributions for numbers changing between
## compartments:
n_SI_tot[,] <- rbinom(S[i,j], sum(p_SI[i,j,]))
rel_strain[,,] <- p_SI[i,j,k]/sum(p_SI[i,j,])

# HARD CODED 2 strains
n_SI[,,] <- if(k==1 || n_strain==1) rbinom(n_SI_tot[i,j],rel_strain[i,j,k]) else
              (if (k==2) n_SI_tot[i,j] - n_SI[i,j,1] else 0)
n_RI[,,2] <- rbinom(R[i,j,2], 1 - exp(-sum(lambda_ij[i,j,,,1])*cross_protection[2,1]* dt))
n_RI_op[,,1] <- n_RI[i,j,2]
n_RI[,,1] <- rbinom(R[i,j,1], 1 - exp(-sum(lambda_ij[i,j,,,1])*cross_protection[1,2]* dt))
n_RI_op[,,2] <- n_RI[i,j,1]

n_IR[,,] <- rbinom(I[i,j,k], p_IR)


lambda_ij[,,,,] <- beta*beta_strain[i5] * mixing_matrix[i,l]*I[k,l,i5]/sum(N)*susceptibility[i,j,i5]*transmisibility[k,l,i5]

#n_import[] <- rbinom(S[i], import/sum(S))

## Total population size
N[,] <- S[i,j] + sum(I[i,j,]) + sum(R[i,j,])

## Initial states:
initial(S[,]) <- S_ini[i,j]
initial(I[,,]) <- I_ini[i,j,k]
initial(R[,,]) <- R_ini[i,j,k]
initial(inc[,,]) <- 0

dim(S) <- c(n,n_vac)
dim(I) <- c(n,n_vac,n_strain)
dim(R) <- c(n,n_vac,n_strain)
dim(N) <- c(n,n_vac)
dim(inc) <- c(n,n_vac,n_strain)
dim(p_SI) <- c(n,n_vac,n_strain)
dim(n_SI) <- c(n, n_vac, n_strain)
dim(n_SI_tot) <- c(n, n_vac)
dim(n_IR) <- c(n,n_vac,n_strain)
dim(n_RI) <- c(n,n_vac,n_strain)
dim(n_RI_op) <- c(n,n_vac,n_strain)

dim(rel_strain) <- c(n,n_vac,n_strain)
#dim(n_import) <- n
dim(lambda_ij) <- c(n,n_vac, n,n_vac, n_strain)
dim(beta_norm) <- n
dim(susceptibility) <- c(n, n_vac, n_strain)
dim(transmisibility) <- c(n, n_vac, n_strain)
dim(cross_protection) <- c(n_strain, n_strain)


## User defined parameters - default in parentheses:
gamma <- user(0.1)

n <- user(4)
S_ini[,] <- user()
I_ini[,,] <- user()
R_ini[,,] <- user()
n_vac <- user()
n_strain <- user()
beta_strain[] <- user()
dim(beta_strain) <- n_strain
mixing_matrix[,] <- user()
dim(mixing_matrix) <- c(n,n)
dim(S_ini) <- c(n, n_vac)
dim(I_ini) <- c(n, n_vac, n_strain)
dim(R_ini) <- c(n, n_vac, n_strain)
beta_norm[] <- user()
susceptibility[,,] <- user()
transmisibility[,,] <- user()
cross_protection[,] <- user()
