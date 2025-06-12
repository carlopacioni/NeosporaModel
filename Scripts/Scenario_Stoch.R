library(data.table)
library(ggplot2)

source(file.path("Scripts","Proc_Neo_IBM.R"))

alpha <- 0.39
InitPrev <- 0.32
zeta_col <- 0.2
zeta_env <- 0.03/alpha
rhov <- 0.43
betas <- 0.05
betaI <- 0.14
eps <- 0.095
g <- 1
K <- 150
Se <- 0.98
Sp <- 0.99

intro <- 0; nsim <- 100; tot.time<- 10; dir.in <- "Results_Stochastic"
#### Baseline ####
parms_BLA <- expand.grid(list(maxAge=12,
                          alpha=alpha,
                          betas=betas,
                          betaI=betaI,
                          rhov=rhov,
                          delta=0.02,
                          eps=eps,
                          sigma=0.07,   
                          zeta_env=zeta_env,
                          zeta_col=zeta_col,
                          p=0.3,
                          g=g,
                          InitPrev=InitPrev,
                          K=K,
                          c=0,        
                          Se=Se,
                          Sp=Sp))

system.time(
  BLA<- proc_IBM(dir.in, intro, nsim, tot.time, parms_BLA, ageI=2, 
                  root_name="BLA", ncore="auto")
)

parms_BLB <- expand.grid(list(maxAge=12,
                              alpha=alpha,
                              betas=betas,
                              betaI=betaI,
                              rhov=rhov,
                              delta=0.3,
                              eps=eps,
                              sigma=0.07,   
                              zeta_env=zeta_env,
                              zeta_col=zeta_col,
                              p=0.3,
                              g=g,
                              InitPrev=InitPrev,
                              K=K,
                              c=0,        
                              Se=Se,
                              Sp=Sp))

system.time(
  BLB<- proc_IBM(dir.in, intro, nsim, tot.time, parms_BLB, ageI=2, 
                 root_name="BLB", ncore="auto")
)

parms_BioS <- expand.grid(list(maxAge=12,
                              alpha=alpha,
                              betas=betas,
                              betaI=betaI,
                              rhov=rhov,
                              delta=0.3,
                              eps=eps,
                              sigma=0.03,   
                              zeta_env=zeta_env,
                              zeta_col=zeta_col,
                              p=0.3,
                              g=g,
                              InitPrev=InitPrev,
                              K=K,
                              c=0,        
                              Se=Se,
                              Sp=Sp))

system.time(
  BioS<- proc_IBM(dir.in, intro, nsim, tot.time, parms_BioS, ageI=2, 
                 root_name="BioS", ncore="auto")
)

#### Test New Animals ####
parms_TNA <- expand.grid(list(maxAge=12,
                          alpha=alpha,
                          betas=betas,
                          betaI=betaI,
                          rhov=rhov,
                          delta=0.02,
                          eps=eps,
                          sigma=0.07,   
                          zeta_env=zeta_env,
                          zeta_col=zeta_col,
                          p=0.002,
                          g=g,
                          InitPrev=InitPrev,
                          K=K,
                          c=0,        
                          Se=Se,
                          Sp=Sp))

system.time(
  TNA <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                    root_name="TNA", ncore="auto")
)

parms_TNB <- expand.grid(list(maxAge=12,
                              alpha=alpha,
                              betas=betas,
                              betaI=betaI,
                              rhov=rhov,
                              delta=0.3,
                              eps=eps,
                              sigma=0.07,   
                              zeta_env=zeta_env,
                              zeta_col=zeta_col,
                              p=0.002,
                              g=g,
                              InitPrev=InitPrev,
                              K=K,
                              c=0,        
                              Se=Se,
                              Sp=Sp))

system.time(
  TNB <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                  root_name="TNB", ncore="auto")
)

parms_TNBioS <- expand.grid(list(maxAge=12,
                              alpha=alpha,
                              betas=betas,
                              betaI=betaI,
                              rhov=rhov,
                              delta=0.02,
                              eps=eps,
                              sigma=0.03,   
                              zeta_env=zeta_env,
                              zeta_col=zeta_col,
                              p=0.002,
                              g=g,
                              InitPrev=InitPrev,
                              K=K,
                              c=0,        
                              Se=Se,
                              Sp=Sp))

system.time(
  TNBioS <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                  root_name="TNBioS", ncore="auto")
)

#### Test and Cull ####
parms_TCA <- expand.grid(list(maxAge=12,
                              alpha=alpha,
                              betas=betas,
                              betaI=betaI,
                              rhov=rhov,
                              delta=0.02,
                              eps=eps,
                              sigma=0.07,   
                              zeta_env=zeta_env,
                              zeta_col=zeta_col,
                              p=0.002,
                              g=g,
                              InitPrev=InitPrev,
                              K=K,
                              c=0.50,        
                              Se=Se,
                              Sp=Sp))

system.time(
  TCA <- proc_IBM(dir.in, intro, nsim, tot.time, parms_TCA, ageI=2, 
                    root_name="TCA", ncore="auto")
)

parms_TCB <- expand.grid(list(maxAge=12,
                              alpha=alpha,
                              betas=betas,
                              betaI=betaI,
                              rhov=rhov,
                              delta=0.3,
                              eps=eps,
                              sigma=0.07,   
                              zeta_env=zeta_env,
                              zeta_col=zeta_col,
                              p=0.002,
                              g=g,
                              InitPrev=InitPrev,
                              K=K,
                              c=0.50,        
                              Se=Se,
                              Sp=Sp))

system.time(
  TCB <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                  root_name="TCB", ncore="auto")
)

parms_TCBioS <- expand.grid(list(maxAge=12,
                              alpha=alpha,
                              betas=betas,
                              betaI=betaI,
                              rhov=rhov,
                              delta=0.02,
                              eps=eps,
                              sigma=0.03,   
                              zeta_env=zeta_env,
                              zeta_col=zeta_col,
                              p=0.002,
                              g=g,
                              InitPrev=InitPrev,
                              K=K,
                              c=0.50,        
                              Se=Se,
                              Sp=Sp))

system.time(
  TCBioS <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                  root_name="TCBioS", ncore="auto")
)


debug(proc_IBM)
