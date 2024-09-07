library(data.table)
library(ggplot2)

source(file.path("Scripts","Proc_Neo_IBM.R"))

intro <- 0; nsim <- 5; tot.time<- 20; dir.in <- "Results_Stochastic"
parms <- expand.grid(list(maxAge=5,
                          alpha=0.3,
                          betas=0.02,
                          betaI=0.08,
                          rhov=0.9,
                          delta=0.1,
                          eps=0.1,
                          sigma=0,
                          zeta=0.028,
                          p=0.3,
                          g=1,
                          InitPrev=0.3,
                          K=1000))

proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, root_name="Sigma0")

debug(proc_IBM)
