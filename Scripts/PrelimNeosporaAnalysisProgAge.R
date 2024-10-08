source(file.path("Scripts", "NeoModelProgAge.R"))
dir.create("PremResultsAge", showWarnings = FALSE)
  
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
                           K=1000,
                           c=0,
                           Se=0.99,
                           Sp=0.99))

res_eps <- proc_res("PremResultsAge", parms, plot_name = "eps_model_plot.png")

# return a list where the first element is a list with the PremResultsAge and the second 
# is a plot where all the plots of each parameter combinations are combined
# the list of PremResultsAge has 3 elements, the first are the parameter values, the second
# is the result from the model and the third is the plot 
res_eps


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
                          K=1000,
                          c=0,
                          Se=0.99,
                          Sp=0.99))

res_sigma <- proc_res("PremResultsAge", parms, plot_name = "sigma_model_plot.png")
res_sigma

#### initial prevalence ####
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
                          InitPrev=c(0.3, 0.5,0.7),
                          K=1000,
                          c=0,
                          Se=0.99,
                          Sp=0.99))

res_initPrev <- proc_res("PremResultsAge", parms, plot_name = "initPrev_model_plot.png")
res_initPrev

#### initial N ####
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
                          K=c(100, 300, 1000),
                          c=0,
                          Se=0.99,
                          Sp=0.99))

res_InitN <- proc_res("PremResultsAge", parms, plot_name = "initN_model_plot.png")
res_InitN

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
                          K=c(1000),
                          c=0,
                          Se=0.99,
                          Sp=0.99))

res_NoNeo <- proc_res("PremResultsAge", parms, plot_name = "NoNeo_model_plot.png")
res_NoNeo

#### one param at the time ####
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
                          InitPrev=1,
                          K=c(1000),
                          c=0,
                          Se=0.99,
                          Sp=0.99))

proc_res("PremResultsAge", parms, plot_name = "OneP_model_plot.png")

#### test culling ####
parms <- expand.grid(list(maxAge=5,
                          alpha=0.3,
                          betas=0.02,
                          betaI=0.08,
                          rhov=0.9,
                          delta=0.1,
                          eps=0.1,
                          sigma=0,
                          zeta=0.16,
                          p=0.2,
                          g=1,
                          InitPrev=0.3,
                          K=c(1000),
                          c=0.9,
                          Se=0.99,
                          Sp=0.99))

proc_res("PremResultsAge", parms, plot_name = "Test_culling_model_plot.png")
