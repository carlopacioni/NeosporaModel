source(file.path("Scripts", "NeoModelProgAge.R"))
# dir.create("PremResultsAge", showWarnings = FALSE)
#### Scenarios  ####

alpha <- 0.43
InitPrev <- 0.32
zeta_col <- 0.2

### Do Nothing  
parms_BLA <- expand.grid(list(maxAge=12,
                              alpha=alpha,
                              betas=0.05,
                              betaI=0.14,
                              rhov=0.43,
                              delta=0.02,
                              eps=0.095,
                              sigma= 0.07,   
                              zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env 
                              zeta_col=zeta_col,     
                              p=0.3,
                              g=1,
                              InitPrev=InitPrev,
                              K=150,
                              c=0,        
                              Se=0.98,
                              Sp=0.99))

res_BLA <- proc_res("PremResultsAge", parms_BLA, plot_name = "BLA_model_plot.png")

parms_BLB <- expand.grid(list(maxAge=12,
                              alpha=alpha,
                              betas=0.05,
                              betaI=0.14,
                              rhov=0.43,
                              delta=0.1,
                              eps=0.095,
                              sigma= 0.07,   
                              zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env 
                              zeta_col=zeta_col,     
                              p=0.3,
                              g=1,
                              InitPrev=InitPrev,
                              K=150,
                              c=0,        
                              Se=0.98,
                              Sp=0.99))

res_BLB <- proc_res("PremResultsAge", parms_BLB, plot_name = "BLB_model_plot.png")

parms_BioS <- expand.grid(list(maxAge=12,
                          alpha=0.39,
                          betas=0.05,
                          betaI=0.14,
                          rhov=0.43,
                          delta=0.02,
                          eps=0.095,
                          sigma= 0.03,
                          zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env 
                          zeta_col=zeta_col,
                          p=0.3,
                          g=1,
                          InitPrev=0.32,
                          K=150,
                          c=0,
                          Se=0.98,
                          Sp=0.99))

res_BioS <- proc_res("PremResultsAge", parms_BioS, plot_name = "BioS_model_plot.png")

### Test New Animals  
parms_TNA <- expand.grid(list(maxAge=12,
                              alpha=0.39,
                              betas=0.05,
                              betaI=0.14,
                              rhov=0.43,
                              delta=0.02,
                              eps=0.095,
                              sigma= 0.07,   
                              zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env 
                              zeta_col=zeta_col,      
                              p=0.002,           
                              g=1,
                              InitPrev=0.32,
                              K=150,
                              c=0,        
                              Se=0.98,
                              Sp=0.99))

res_TNA <- proc_res("PremResultsAge", parms_TNA, plot_name = "TNA_model_plot.png")

parms_TNB <- expand.grid(list(maxAge=12,
                              alpha=0.39,
                              betas=0.05,
                              betaI=0.14,
                              rhov=0.43,
                              delta=0.1,
                              eps=0.095,
                              sigma= 0.07,   
                              zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env 
                              zeta_col=zeta_col,      
                              p=0.002,           
                              g=1,
                              InitPrev=0.32,
                              K=150,
                              c=0,        
                              Se=0.98,
                              Sp=0.99))

res_TNB <- proc_res("PremResultsAge", parms_TNB, plot_name = "TNB_model_plot.png")


parms_TNBioS <- expand.grid(list(maxAge=12,
                              alpha=0.39,
                              betas=0.05,
                              betaI=0.14,
                              rhov=0.43,
                              delta=0.02,
                              eps=0.095,
                              sigma= 0.03,   
                              zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env 
                              zeta_col=zeta_col,       
                              p=0.002,           
                              g=1,
                              InitPrev=0.32,
                              K=150,
                              c=0,        
                              Se=0.98,
                              Sp=0.99))

res_TNBioS <- proc_res("PremResultsAge", parms_TNBioS, plot_name = "TNBioS_model_plot.png")

### Test and Cull  
parms_TNC_A <- expand.grid(list(maxAge=12,
                                alpha=0.39,
                                betas=0.05,
                                betaI=0.14,
                                rhov=0.43,
                                delta=0.02,
                                eps=0.095,
                                sigma= 0.07,   
                                zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env 
                                zeta_col=zeta_col,     
                                p=0.002,
                                g=1,
                                InitPrev=0.32,
                                K=150,
                                c=0.50,        
                                Se=0.98,
                                Sp=0.99))

res_TNC_A <- proc_res("PremResultsAge", parms_TNC_A, plot_name = "TNC_A_model_plot.png")

parms_TNC_B <- expand.grid(list(maxAge=12,
                                alpha=0.39,
                                betas=0.05,
                                betaI=0.14,
                                rhov=0.43,
                                delta=0.1,
                                eps=0.095,
                                sigma= 0.07,   
                                zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env 
                                zeta_col=zeta_col,     
                                p=0.002,
                                g=1,
                                InitPrev=0.32,
                                K=150,
                                c=0.50,        
                                Se=0.98,
                                Sp=0.99))

res_TNC_B <- proc_res("PremResultsAge", parms_TNC_B, plot_name = "TNC_B_model_plot.png")


parms_TNC_BioS <- expand.grid(list(maxAge=12,
                                alpha=0.39,
                                betas=0.05,
                                betaI=0.14,
                                rhov=0.43,        
                                delta=0.02,
                                eps=0.095,
                                sigma= 0.03,  
                                zeta_env=0.03/alpha, # to match publish value beta=0.03=alpha*zeta_env 
                                zeta_col=zeta_col,      
                                p=0.002,
                                g=1,
                                InitPrev=0.32,
                                K=150,
                                c=0.50,        
                                Se=0.98,
                                Sp=0.99))

res_TNC_BioS <- proc_res("PremResultsAge", parms_TNC_BioS, plot_name = "TNC_BioS_model_plot.png")

extract_prev <- function(l) {
  dl <- vector("list", length = length(l))
  for(i in seq_along(l)) {
    temp <- l[[i]][[1]][[1]][[2]][, .(time, Prev)]
    name <- names(l)[i]
    temp[, Scenario := name]
    temp[, HerdManagement := ifelse(substr(name, 1, 3) == "Int", "Internal", "External")]
    temp[, Treatment := substr(name, 4, 5)]
    temp[, BioSecurity := ifelse(substr(name, 6, 6) == "Y", "Yes", "No")]
    dl[[i]] <- temp
  }
  return(rbindlist(dl))
}

data <- extract_prev(l)

l <- list(
  IntBLN = res_BLA,
  ExtBLN = res_BLB,
  IntBLY = res_BioS,
  IntTNN = res_TNA,
  ExtTNN = res_TNB,
  IntTNY = res_TNBioS,
  IntTCN = res_TNC_A,
  ExtTCN = res_TNC_B,
  IntTCY = res_TNC_BioS
)

# Plot
ggplot(data, aes(x = time, y = Prev,
                 color = Treatment,
                 linetype = HerdManagement,
                 shape = BioSecurity)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Seroprevalence Over Time by Management, Treatment, and Biosecurity",
    x = "Years",
    y = "Seroprevalence (%)"
  ) +
  theme_minimal()
