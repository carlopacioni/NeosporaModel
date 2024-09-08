proc_IBM <- function(dir.in, intro, nsim, tot.time, params, ageI, root_name, ncore="auto") {
  dir.create(dir.in, showWarnings = FALSE)
  if(ncore == "auto") ncore <- detectCores()
  if(ncore>nsim) ncore <- nsim
  res_out<-vector("list", length = nrow(params))
  for(rn in seq_len(nrow(params))) {
    
     parms <- list(maxAge=params[rn, "maxAge"],
                                alpha=params[rn, "alpha"],
                                betas=params[rn, "betas"],
                                betaI=params[rn, "betaI"],
                                rhov=params[rn, "rhov"],
                                delta=params[rn, "delta"],
                                eps=params[rn, "eps"], 
                                sigma=params[rn, "sigma"],
                                zeta=params[rn, "zeta"],
                                p=params[rn, "p"],
                                g=params[rn, "g"],
                                InitPrev=params[rn, "InitPrev"],
                                K=params[rn, "K"])
    
     # S, I
     init.pop<- c(parms$K * (1 - parms$InitPrev), parms$K * parms$InitPrev)
     
     if(ncore == 1) {
       source(file.path("Scripts","Neo IBM.R"))
       out <- lapply(1:nsim, function(z){Neo.ibm(popsize, init.pop, tot.time, 
                                                 intro, ageI=ageI, parms)})
     } else {
     
     library(parallel)
     cl <- makeCluster(5)
     on.exit(stopCluster(cl))
     clusterExport(cl, varlist = c("init.pop", "tot.time", "intro", "parms"), 
                   envir = environment())
     out <- parLapply(cl=cl, 1:nsim, function(z){
       source(file.path("Scripts","Neo IBM.R"))
       Neo.ibm(popsize, init.pop, tot.time, intro, ageI=ageI, parms)
     })
     }
     res <- rbindlist(out, idcol = "Iter")
     res[, N:= S + I]
     res[, Prev := I/N]
     
     p <- ggplot(res, aes(time, Prev) ) + 
       geom_line(aes(group=Iter), col="red", alpha=0.6) + 
       stat_summary(fun = "median", geom = "line", color = "black", size = 1.2) +
       stat_summary(fun = quantile, fun.args = list(probs=0.0275), geom = "line", 
                    color = "black", size = 1.2, linetype=2) +
       stat_summary(fun = quantile, fun.args = list(probs=0.975), geom = "line", 
                    color = "black", size = 1.2, linetype=2) +
       theme(legend.position = "none")
     
     # summary stats - needs clean up
     suppressWarnings(
     res_long <- melt.data.table(res, id.vars = "time", variable.name = "Parameter")
     )
     
     res_summary <- res_long[, .(Mean=mean(value), SD=sd(value), Median=median(value), 
                                 lcl=quantile(value, 0.0275, na.rm=TRUE),
                                 ucl=quantile(value,0.975, na.rm=TRUE)), by=list(time, Parameter)]
     
     new_parm <- c(parms, intro=intro, nsim=nsim, tot.time=tot.time, ageI=ageI)
     res_out[[rn]] <- list(new_parm, res, res_long, res_summary, p)
  }
  
  library("cowplot")
  pgrid <- plot_grid(plotlist = lapply(res_out, "[[", 5),
                     labels = seq_len(nrow(params)),
                     ncol = 1)
  #pgrid
  save_plot(file.path(dir.in, paste0(root_name, "_plot.png")), plot = pgrid, nrow = 2,
            base_height = 7, base_asp = 3)
  list_out <- list(res_out, pgrid) 
  save(list_out, file = file.path(dir.in, paste0(root_name, "_data.Rda")))
  return(list_out)
}



