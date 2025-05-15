source(file.path("Scripts", "NeoModelProgAge.R"))
dir.create("PremPremResults", showWarnings = FALSE)
  
  #### eps ####
parms <- expand.grid(list(maxAge=12,
                          alpha=0.43,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=seq(0.1, 0.4, by=0.1), ####
                          sigma= 0.07,   
                          zeta=0.03/0.43, # to match publish value beta=0.03=alpha*zeta      
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
                          alpha=0.43,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=seq(0.02, 0.4, by=0.1),
                          eps=0.1, ####
                          sigma= 0.07,   
                          zeta=0.03/0.43, # to match publish value beta=0.03=alpha*zeta      
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

#### sigma ####
parms <- expand.grid(list(maxAge=12,
                          alpha=0.43,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=0.1, ####
                          sigma=c(0, 0.005, 0.01, 0.02), 
                          zeta=0.03/0.43, # to match publish value beta=0.03=alpha*zeta      
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
                          alpha=0.43,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=0.1, ####
                          sigma=0.0, #
                          zeta=0.03/0.43, # to match publish value beta=0.03=alpha*zeta      
                          p=0.3,
                          g=1,
                          InitPrev=c(0, 0.3, 0.5),
                          K=150, ###
                          c=0,        
                          Se=0.98,
                          Sp=0.99))

res_initPrev <- proc_res("PremResults", parms, plot_name = "initPrev_model_plot.png")
res_initPrev

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
                          zeta=c(0, 0.025/0.2), # to match publish value beta=0.03=alpha*zeta      
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
