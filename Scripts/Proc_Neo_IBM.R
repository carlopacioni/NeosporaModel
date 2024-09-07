proc_IBM <- function(dir.in, intro, nsim, tot.time, params, ageI, root_name) {
  dir.create(dir.in, showWarnings = FALSE)
  if(any(parms[, "maxAge"]< 4)) stop("maxAge is assumed to be 4 or more years")
  res_out<-vector("list", length = nrow(parms))
  for(rn in seq_len(nrow(parms))) {
    
     parms <- list(maxAge=parms[rn, "maxAge"],
                                alpha=parms[rn, "alpha"],
                                betas=parms[rn, "betas"],
                                betaI=parms[rn, "betaI"],
                                rhov=parms[rn, "rhov"],
                                delta=parms[rn, "delta"],
                                eps=parms[rn, "eps"], 
                                sigma=parms[rn, "sigma"],
                                zeta=parms[rn, "zeta"],
                                p=parms[rn, "p"],
                                g=parms[rn, "g"],
                                InitPrev=parms[rn, "InitPrev"],
                                K=parms[rn, "K"])
    
     # S, I
     init.pop<- c(parms$K * (1 - parms$InitPrev), parms$K * parms$InitPrev)
     
     #can parallelise the following
     #out <- lapply(1:nsim, function(z){Neo.ibm(popsize, init.pop, tot.time, intro, ageI=2, parms)})
     
     library(parallel)
     cl <- makeCluster(5)
     on.exit(stopCluster(cl))
     clusterExport(cl, varlist = c("init.pop", "tot.time", "intro", "parms"), 
                   envir = environment())
     out <- parLapply(cl=cl, 1:nsim, function(z){
       source(file.path("Scripts","Neo IBM.R"))
       Neo.ibm(popsize, init.pop, tot.time, intro, ageI=ageI, parms)
     })
     
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
     
     res_out[[rn]] <- list(c(parms, intro=intro, nsim=nsim, tot.time=tot.time, ageI=ageI), 
                           res, res_long, res_summary, p)
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



