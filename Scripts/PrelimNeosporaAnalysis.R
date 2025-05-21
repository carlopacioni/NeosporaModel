source(file.path("Scripts", "NeoModelProgAge.R"))
dir.create("PremPremResults", showWarnings = FALSE)

alpha <- 0.43

# Prelim estimates from the farms indicate that delta + rhoH = 0.0145
# assuming that alpha*zeta_env=0.03, with a Prev=0.3 --> rhoH = 0.009 
# hence sigma = 0.0145 - 0.009 = 0.00435

  
  #### eps ####
parms <- expand.grid(list(maxAge=12,
                          alpha=alpha,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=seq(0.1, 0.4, by=0.1), ####
                          sigma= 0.07,   
                          zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env 
                          zeta_col=0.03/alpha,
                          p=0.3,
                          g=1,
                          InitPrev=0.32,
                          K=150, ###
                          c=0,        
                          Se=0.98,
                          Sp=0.99))

res_eps <- proc_res("PremResults", parms, plot_name = "eps_model_plot.png")
res_eps[[1]][[1]][[2]][, N]
res_eps[[1]][[2]][[2]][, N]
res_eps[[1]][[3]][[2]][, N]
res_eps[[1]][[4]][[2]][, N]
# return a list where the first element is a list with the PremResults and the second 
# is a plot where all the plots of each parameter combinations are combined
# the list of PremResults has 3 elements, the first are the parameter values, the second
# is the result from the model and the third is the plot 
res_eps

#### delta ####
parms <- expand.grid(list(maxAge=12,
                          alpha=alpha,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=seq(0.02, 0.4, by=0.1),
                          eps=0.1, ####
                          sigma= 0.07,   
                          zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env      
                          zeta_col=0.03/alpha,
                          p=0.3,
                          g=1,
                          InitPrev=0.32,
                          K=150, ###
                          c=0,        
                          Se=0.98,
                          Sp=0.99))

res_delta <- proc_res("PremResults", parms, plot_name = "delta_model_plot.png")
res_delta[[1]][[1]][[2]][, N]
res_delta[[1]][[4]][[2]][, N]

#### Culling ####
parms <- expand.grid(list(maxAge=12,
                          alpha=alpha,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=c(0.02, 0.4),
                          eps=0.1, ####
                          sigma= 0.07,   
                          zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env      
                          zeta_col=0.03/alpha,
                          p=0.3,
                          g=1,
                          InitPrev=0.32,
                          K=150, ###
                          c=0.5,        
                          Se=0.98,
                          Sp=0.99))

res_culling <- proc_res("PremResults", parms, plot_name = "culling_model_plot.png")
res_culling[[1]][[1]][[2]][, N]
res_culling[[1]][[2]][[2]][, N]

#### sigma ####
parms <- expand.grid(list(maxAge=12,
                          alpha=alpha,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=0.1, ####
                          sigma=c(0, 0.005, 0.01, 0.02), 
                          zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env      
                          zeta_col=0.03/alpha,
                          p=0.3,
                          g=1,
                          InitPrev=0.32,
                          K=150, ###
                          c=0,        
                          Se=0.98,
                          Sp=0.99))

res_sigma <- proc_res("PremResults", parms, plot_name = "sigma_model_plot.png")
res_sigma

#### initial prevalence ####
parms <- expand.grid(list(maxAge=12,
                          alpha=alpha,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=0.1, ####
                          sigma=0.0, #
                          zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env      
                          zeta_col=0.03/alpha,
                          p=0.3,
                          g=1,
                          InitPrev=c(0, 0.3, 0.5),
                          K=150, ###
                          c=0,        
                          Se=0.98,
                          Sp=0.99))

res_initPrev <- proc_res("PremResults", parms, plot_name = "initPrev_model_plot.png")
res_initPrev

#### test zeta_col ####
parms <- expand.grid(list(maxAge=12,
                          alpha=alpha,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=0.1, ####
                          sigma=0.0, #
                          zeta_env=c(0, 0.03/alpha), # to match publish value beta=0.03=alpha*zeta_env      
                          zeta_col=c(0, 0.03/alpha),
                          p=0.,
                          g=1,
                          InitPrev=0.3,
                          K=150, ###
                          c=0,        
                          Se=0.98,
                          Sp=0.99))

res_zeta_col <- proc_res("PremResults", parms, times=0:50, plot_name = "zeta_col_model_plot.png")
res_zeta_col
res_zeta_col[[1]][[2]][[2]][, Prev]
res_zeta_col[[1]][[3]][[2]][, Prev]
res_zeta_col[[1]][[4]][[2]][, Prev]

#### French et al ####
# rho1=alpha(1-betaI)=0.2 from the paper, so derive alpha keeping our abortion rate
parms <- expand.grid(list(maxAge=12,
                          alpha=0.2/(1-0.16),
                          betas=1-(0.2157/(0.2/(1-0.16))),
                          betaI=0.16,
                          rhov=0.95,
                          delta=0.3,
                          eps=0., ####
                          sigma=c(0, 0.025), 
                          zeta_env=c(0, 0.025/0.2), # to match publish value beta=0.03=alpha*zeta_env      
                          zeta_col=0.03/alpha,
                          p=0,
                          g=1,
                          InitPrev=0.02,
                          K=150, ###
                          c=0,        
                          Se=0.98,
                          Sp=0.99))


res_French <- proc_res("PremResults", parms, times=0:50, plot_name = "French_plot.png")
res_French
parms
res_French[[1]][[1]][[2]][, N]
res_French[[1]][[4]][[2]][, Prev]
