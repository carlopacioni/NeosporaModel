library(data.table)
library(ggplot2)

source(file.path("Scripts","Proc_Neo_IBM.R"))

intro <- 0; nsim <- 100; tot.time<- 20; dir.in <- "Results_Stochastic" 
alpha <- 0.43

#### eps ####
parms <- expand.grid(list(maxAge=12,
                          alpha=alpha,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=seq(0.1, 0.4, by=0.1), # 0.095
                          sigma=0,
                          zeta_env=0.03/alpha,
                          zeta_col=0.03/alpha,
                          p=0.3,
                          g=1,
                          InitPrev=0.3,
                          K=150,              
                          c=0,        
                          Se=0.98,
                          Sp=0.99))

system.time(
eps <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                   root_name="eps", ncore="auto")
)

#### sigma ####
parms <- expand.grid(list(maxAge=12,
                          alpha=alpha,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=0.1,
                          sigma=c(0, 0.005, 0.07),
                          zeta_env=0.03/alpha,
                          zeta_col=0.03/alpha,
                          p=0.3,
                          g=1,
                          InitPrev=0.3,
                          K=150,
                          c=0, 
                          Se=0.98, 
                          Sp=0.99))

system.time(
  sigma <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                  root_name="sigma", ncore="auto")
)

#### sigma vs p ####
parms <- expand.grid(list(maxAge=12,
                          alpha=alpha,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=0.1,
                          sigma=c(0, 0.005, 0.07),
                          zeta_env=0.03/alpha,
                          zeta_col=0.03/alpha,
                          p=c(0, 0.3),
                          g=1,
                          InitPrev=0.3,
                          K=150,
                          c=0,
                          Se=0.98,
                          Sp=0.99))

system.time(
  sigma <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                    root_name="sigmaVsP", ncore="auto")
)


#### No Neo ####
parms <- expand.grid(list(maxAge=12,
                          alpha=alpha,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=0.1,
                          sigma=0,
                          zeta_env=0.03/alpha,
                          zeta_col=0.03/alpha,
                          p=0.,
                          g=1,
                          InitPrev=0.0,
                          K=150,
                          c=0,        
                          Se=0.98,
                          Sp=0.99))

system.time(
  NoNeo <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                    root_name="NoNeo", ncore="auto")
)


#### No Neo with intro ####
parms <- expand.grid(list(maxAge=12,
                          alpha=alpha,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=0.1,
                          sigma=0,
                          zeta_env=0.03/alpha,
                          zeta_col=0.03/alpha,
                          p=0.,
                          g=1,
                          InitPrev=0.0,
                          K=150,
                          c=0,        
                          Se=0.98,
                          Sp=0.99))

system.time(
  NoNeoIntro <- proc_IBM(dir.in, intro=2, nsim, tot.time, parms, ageI=2, 
                    root_name="NoNeoIntro", ncore="auto")
)

#### French et al ####
parms <- expand.grid(list(maxAge=12,
                          alpha=0.2/(1-0.16),
                          betas=1-(0.2157/(0.2/(1-0.16))),
                          betaI=0.16,
                          rhov=0.95,
                          delta=0.3,
                          eps=0., ####
                          sigma=c(0, 0.025), 
                          zeta_env=c(0, 0.025/0.2), # to match publish value beta=0.03=alpha*zeta_env      
                          zeta_col=0.0/alpha,
                          p=0,
                          g=1,
                          InitPrev=0.02,
                          K=150, ###
                          c=0,        
                          Se=0.98,
                          Sp=0.99))


system.time(
  res_French <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                         root_name="res_French", ncore="auto")
)
res_French
parms





debug(proc_IBM)
debug(Neo.ibm)
debug(schedule)
debug(age.animals)
debug(keepNconstant)
