library(data.table)
library(ggplot2)

source(file.path("Scripts","Proc_Neo_IBM.R"))

intro <- 0; nsim <- 100; tot.time<- 20; dir.in <- "Results_Stochastic"
#### eps ####
parms <- expand.grid(list(maxAge=5,
                          alpha=0.3,
                          betas=0.02,
                          betaI=0.08,
                          rhov=0.9,
                          delta=0.1,
                          eps=seq(0.1, 0.4, by=0.1), # 0.095
                          sigma=0,
                          zeta=0.028,
                          p=0.3,
                          g=1,
                          InitPrev=0.3,
                          K=1000))

system.time(
eps <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                   root_name="eps", ncore="auto")
)

#### sigma ####
parms <- expand.grid(list(maxAge=5,
                          alpha=0.3,
                          betas=0.02,
                          betaI=0.08,
                          rhov=0.9,
                          delta=0.1,
                          eps=0.1,
                          sigma=c(0, 0.005, 0.01, 0.02),
                          zeta=0.028,
                          p=0.3,
                          g=1,
                          InitPrev=0.3,
                          K=1000))

system.time(
  sigma <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                  root_name="sigma", ncore="auto")
)

#### No Neo ####
parms <- expand.grid(list(maxAge=5,
                          alpha=0.3,
                          betas=0.02,
                          betaI=0.08,
                          rhov=0.,
                          delta=0.1,
                          eps=0.1,
                          sigma=0,
                          zeta=0.0,
                          p=0.,
                          g=1,
                          InitPrev=0.0,
                          K=c(1000)))

system.time(
  NoNeo <- proc_IBM(dir.in, intro, nsim, tot.time, parms, ageI=2, 
                    root_name="NoNeo", ncore="auto")
)
res_NoNeo <- proc_res("PremResultsAge", parms, plot_name = "NoNeo_model_plot.png")
res_NoNeo






debug(proc_IBM)
