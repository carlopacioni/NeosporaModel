library(data.table)
library(ggplot2)

source(file.path("Scripts","Proc_Neo_IBM.R"))

Se=0.98
Sp=0.99

intro <- 0; nsim <- 100; tot.time<- 10; dir.in <- "Results_Stochastic"
#### Do Nothing ####
parms <- expand.grid(list(maxAge=12,
                          alpha=0.43,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=0.095,
                          sigma=0.07,   
                          zeta=0.03,      
                          p=0.3,
                          g=1,
                          InitPrev=0.32,
                          K=450,
                          c=0,        
                          Se=0.98,
                          Sp=0.99))

system.time(
  DN <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                  root_name="DN", ncore="auto")
)

#### Test New Animals ####
parms <- expand.grid(list(maxAge=12,
                          alpha=0.43,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=0.095,
                          sigma=0.03, #   
                          zeta=0.03,      
                          p=1-Se,         #  
                          g=1,
                          InitPrev=0.32,
                          K=450,
                          c=0,        
                          Se=0.98,
                          Sp=0.99))

system.time(
  TNA <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                    root_name="TNA", ncore="auto")
)

#### Test and Cull ####
parms <- expand.grid(list(maxAge=12,
                          alpha=0.43,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59, #
                          delta=0.02, #
                          eps=0.095,
                          sigma=0.03, #   
                          zeta=0.03,      
                          p=0.3, #
                          g=1,
                          InitPrev=0.32,
                          K=450,
                          c=0.5,        
                          Se=0.98,
                          Sp=0.99))

system.time(
  TC <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                    root_name="TC", ncore="auto")
)



#### Test NA + T & C ####
parms <- expand.grid(list(maxAge=12,
                          alpha=0.43,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59, #
                          delta=0.02, #
                          eps=0.095,
                          sigma=0.03, #   
                          zeta=0.03,      
                          p=1-Se, #
                          g=1,
                          InitPrev=0.32,
                          K=450,
                          c=0.5,        
                          Se=0.98,
                          Sp=0.99))

system.time(
  TNTC <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                  root_name="TNTC", ncore="auto")
)





debug(proc_IBM)

