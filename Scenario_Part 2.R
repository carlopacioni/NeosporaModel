#### Scenarios  ####

### Do Nothing  
parms_DoN <- expand.grid(list(maxAge=12,
                              alpha=0.43,
                              betas=0.05,
                              betaI=0.16,
                              rhov=0.59,
                              delta=0.02,
                              eps=0.095,
                              sigma= 0.07,   
                              zeta=0.03,      
                              p=0.3,
                              g=1,
                              InitPrev=0.32,
                              K=450,
                              c=0.20,        
                              Se=0.98,
                              Sp=0.99))

res_sigma_DoN <- proc_res("PremResultsAge", parms_DoN, plot_name = "DoN_model_A_plot.png")

parms <- expand.grid(list(maxAge=12,
                          alpha=0.43,
                          betas=0.05,
                          betaI=0.16,
                          rhov=0.59,
                          delta=0.02,
                          eps=0.095,
                          sigma= 0.03,
                          zeta=0.03,
                          p=0,
                          g=1,
                          InitPrev=0.32,
                          K=450,
                          c=0.36,
                          Se=0.98,
                          Sp=0.99))

res_sigma <- proc_res("PremResultsAge", parms, plot_name = "DoN_model_B_plot.png")
res_sigma

### Test New Animals  
parms_TNA <- expand.grid(list(maxAge=12,
                              alpha=0.43,
                              betas=0.05,
                              betaI=0.16,
                              rhov=0.59,
                              delta=0.02,
                              eps=0.095,
                              sigma= 0.05,   
                              zeta=0.03,      
                              p=0.03,           
                              g=1,
                              InitPrev=0.32,
                              K=450,
                              c=0.20,        
                              Se=0.98,
                              Sp=0.99))

res_sigma_TNA <- proc_res("PremResultsAge", parms_TNA, plot_name = "NewAni_A_model_plot.png")

parms_TNB <- expand.grid(list(maxAge=12,
                              alpha=0.43,
                              betas=0.05,
                              betaI=0.16,
                              rhov=0.59,
                              delta=0.02,
                              eps=0.095,
                              sigma= 0.03,   
                              zeta=0.03,      
                              p=0,           
                              g=1,
                              InitPrev=0.32,
                              K=450,
                              c=0.36,        
                              Se=0.98,
                              Sp=0.99))

res_sigma_TNB <- proc_res("PremResultsAge", parms_TNB, plot_name = "NewAni_B_model_plot.png")

### Test and Cull  
parms_TNC_A <- expand.grid(list(maxAge=12,
                                alpha=0.43,
                                betas=0.05,
                                betaI=0.16,
                                rhov=0,
                                delta=0.02,
                                eps=0.095,
                                sigma= 0.03,   
                                zeta=0.03,      
                                p=0.30,
                                g=1,
                                InitPrev=0.32,
                                K=450,
                                c=0.50,        
                                Se=0.98,
                                Sp=0.99))

res_sigma_TNC_A <- proc_res("PremResultsAge", parms_TNC_A, plot_name = "TNC_A_model_plot.png")

parms_TNC_B <- expand.grid(list(maxAge=12,
                                alpha=0.43,
                                betas=0.05,
                                betaI=0,
                                rhov=0,        
                                delta=0.02,
                                eps=0.095,
                                sigma= 0,  
                                zeta=0.03,      
                                p=0,
                                g=1,
                                InitPrev=0.32,
                                K=450,
                                c=0.50,        
                                Se=0.98,
                                Sp=0.99))

res_sigma_TNC_B <- proc_res("PremResultsAge", parms_TNC_B, plot_name = "TNC_B_model_plot.png")



####EXPLORING

# Adjusted Parameters for Infection Control without Extreme Culling
parms_adjusted <- expand.grid(list(
  maxAge = 12,
  alpha = 0.43,
  betas = 0.05,
  betaI = 0.16,
  rhov = 0.59,
  delta = 0.02,
  eps = 0.095,
  sigma = 0.01,    # Reduced environmental transmission
  zeta = 0.01,     # Reduced within-herd transmission
  p = 0,           # No infected introductions
  g = 1,
  InitPrev = 0.32,
  K = 1000,
  c = 0.20,        # Moderate culling
  Se = 0.98,
  Sp = 0.99
))

res_adjusted <- proc_res("PremResultsAge", parms_adjusted, plot_name = "Adjusted_Infection_Control.png")




# parameter grid varying sigma, zeta, and p values
parms_grid <- expand.grid(
  sigma = c(0.01, 0.03, 0.07),
  zeta = c(0.01, 0.03, 0.10),
  p = c(0, 0.20, 0.50),
  maxAge = 12,
  alpha = 0.43,
  betas = 0.05,
  betaI = 0.16,
  rhov = 0.59,
  delta = 0.02,
  eps = 0.095,
  InitPrev = 0.32,
  K = 1000,
  c = 0.20,   
  g = 1,
  Se = 0.98,
  Sp = 0.99
)

run_scenario <- function(parms, scenario_name) {
  res <- proc_res("PremResultsAge", parms, plot_name = paste(scenario_name, "_plot.png"))
  return(res)
}

