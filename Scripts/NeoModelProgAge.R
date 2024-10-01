library(deSolve)
library(data.table)
library(ggplot2)

#### fit model ####
fitDetNeospora <- function(dir.in,
                           maxAge,
                           times,
                           alpha,
                           betas,
                           betaI,
                           rhov,
                           delta,
                           eps,
                           sigma,
                           zeta,
                           p,
                           g,
                           InitPrev,
                           K,
                           c,
                           Se,
                           Sp){


  SI5yrEq <- function(time, state, params){

    with(as.list(c(state, params)),{
      
      # index of adult cows
      cindex <- 3:(maxAge)
      Sc <- state[cindex]
      tSc <- sum(Sc) # total susceptible cows
      
      Ic <- state[(maxAge + 3):length(state)]
      tIc <- sum(Ic) # total infected cows
      N <- So+Sh+tSc+Io+Ih+tIc
      theta <- K - N
      
      # create empty d* vector so that I can use indeces
      dSc <- Sc
      dIc <- Ic

      if(theta<0) {
        theta_neg <- theta
        theta_pos <- 0
      } else {
        theta_neg <- 0
        theta_pos <- theta
      }

      rhoh <- alpha*zeta*(tIc/(tSc + tIc))
      po <- Io/(Io+So)

      dSo <- (1-po)*theta_neg + alpha*(1-betas)*tSc+alpha*(1-betaI)*(1-rhov)*tIc-
        (delta+rhoh+sigma)*So-g*So - 
        c*(1-Sp)*So # test and culling false positives
      dSh <- (1-p)*theta_pos + g*So-(delta+rhoh+sigma)*Sh-g*Sh - 
        c*(1-Sp)*Sh
      dSc[1] <- g*Sh -(delta+rhoh+sigma)*Sc[1]-g*Sc[1] - c*(1-Sp)*Sc[1]
      
      dIo <- po*theta_neg + alpha*rhov*(1-betaI)*tIc + (rhoh+sigma)*So-delta*Io-
        g*Io - # aging
        c*Se*Io # test & culling detected positives
      dIh <- p*theta_pos+g*Io+(rhoh+sigma)*Sh-delta*Ih-g*Ih- c*Se*Ih
      dIc[1] <- g*Ih+(rhoh+sigma)*Sc[1]-delta*Ic[1]-g*Ic[1] -  c*Se*Ic[1]
      
      for(i in 2:length(dSc)) {
        if(i < length(dSc)) {
          dSc[i] <- g*Sc[i - 1]-(delta+rhoh+sigma)*Sc[i]-g*Sc[i] - c*(1-Sp)*Sc[i]
          dIc[i] <- g*Ic[i - 1]+(rhoh+sigma)*Sc[i]-delta*Ic[i]-g*Ic[i]-c*Se*Ic[i]
        } else {
          dSc[i] <- g*Sc[i - 1]-(delta+rhoh+sigma+eps)*Sc[i] - c*(1-Sp)*Sc[i]
          dIc[i] <- g*Ic[i - 1]+(rhoh+sigma)*Sc[i]-(delta+eps)*Ic[i] -  c*Se*Ic[i]
        }
      }

      return(list(c(dSo, dSh, dSc, dIo, dIh, dIc)))
    })
  }

  params <- c(
    maxAge=maxAge,
    alpha=alpha,
    betas=betas,
    betaI=betaI,
    rhov=rhov,
    delta=delta,
    eps=eps, # 0.095
    sigma=sigma,
    zeta=zeta,
    p=p,
    K=K,
    g=g,
    c=c,
    Se=Se,
    Sp=Sp
  )

  # work out the initial states
  dg <- dgeom(1:maxAge, 0.2)
  Sstate <- (dg/sum(dg))*K*(1-InitPrev)
  So=Sstate[1]; Sh=Sstate[2]; Sc=tail(Sstate, -2)
  
  Istate <- (dg/sum(dg))*K*(InitPrev)
  Io=Istate[1]; Ih=Istate[2]; Ic=tail(Istate, -2)

  initial_state <- c(So=So, Sh=Sh, Sc=Sc, 
                     Io=Io, Ih=Ih, Ic=Ic)


  #sum(initial_state)
  #times <- 0:10

  model <- ode(initial_state, times, SI5yrEq, params)
  #model[1:10,]
  #rowSums(model[,-1])
  res<-data.table(model)
  Scols <- names(res)[grep(pattern = "^S", names(res))]
  Icols <- names(res)[grep(pattern = "^I", names(res))]
  res[, S:= rowSums(.SD), .SDcols=Scols]
  res[, I:= rowSums(.SD), .SDcols=Icols]
  res[, Sc:= rowSums(.SD), .SDcols=tail(Scols, -2)]
  res[, Ic:= rowSums(.SD), .SDcols=tail(Icols, -2)]
  res[, N:=S + I]
  res[, Prev:=I/N]
  res[, rhoh:=alpha*zeta*((Ic)/(Sc + Ic))]
  #res

  res_long <- melt(res, id.vars = "time", variable.name = "Compartment", value.name = "Values")
  res_long[, Group:=substr(Compartment, start = 1, stop = 1)]
  #res_long[, Group := as.character(Group)]
  res_long_tmp <- res_long[Group!="P" & Group!="N",]
  # ggplot(res_long_tmp, aes(time, Values, col=Compartment)) +
  #   geom_line(linewidth=1) +
  #   facet_grid(.~Group)
  # ggsave(file.path(dir.in, "Comparts_plot.png"))

  p<-ggplot(res_long[Compartment=="Prev",], aes(time, Values)) +
    geom_line(linewidth=1)
  #ggsave(file.path(dir.in, "Prev_plot.png"))

  n<-ggplot(res_long[Compartment %in% c("S", "I", "N"),], 
            aes(time, Values, col=Compartment)) +
    geom_line(linewidth=1) + theme(legend.position="none") +
    scale_color_manual(values=c("green", "red", "blue")) # S, I, Total

  #ggsave(file.path(dir.in, "S_I_N_plot.png"))

  library("cowplot")
  pgrid <- plot_grid(p, n,
                     #labels = c("A", "B"),
                     ncol = 2)
  pgrid
  #save_plot(file.path(dir.in, "res_plot.png"), plot = pgrid, nrow = 2)

  return(list(c(params, InitPrev=InitPrev,N=K, Times=max(times)), res, pgrid))

}

#### process results ####
# return a list where the first element is a list with the results and the second
# is a plot where all the plots of each parameter combinations are combined
# the list of results has 3 elements, the first are the parameter values, the second
# is the result from the model and the third is the plot
proc_res <- function(dir.in, parms, plot_name) {
  if(any(parms[, "maxAge"]< 4)) stop("maxAge is assumed to be 4 or more years")
  res<-vector("list", length = nrow(parms))
  for(rn in seq_len(nrow(parms))) {
    res[[rn]] <- fitDetNeospora(dir.in,
                                times=0:20,
                                maxAge=parms[rn, "maxAge"],
                                alpha=parms[rn, "alpha"],
                                betas=parms[rn, "betas"],
                                betaI=parms[rn, "betaI"],
                                rhov=parms[rn, "rhov"],
                                delta=parms[rn, "delta"],
                                eps=parms[rn, "eps"], # 0.095
                                sigma=parms[rn, "sigma"],
                                zeta=parms[rn, "zeta"],
                                p=parms[rn, "p"],
                                g=parms[rn, "g"],
                                InitPrev=parms[rn, "InitPrev"],
                                K=parms[rn, "K"],
                                c=parms[rn, "c"],
                                Se=parms[rn, "Se"],
                                Sp=parms[rn, "Sp"])

  }

  library("cowplot")
  pgrid <- plot_grid(plotlist = lapply(res, "[[", 3),
                     labels = seq_len(nrow(parms)),
                     ncol = 1)
  #pgrid
  save_plot(file.path(dir.in, plot_name), plot = pgrid, nrow = 2,
            base_height = 7, base_asp = 3)
  return(list(res,pgrid))
}
