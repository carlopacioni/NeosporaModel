library(deSolve)
library(data.table)
library(ggplot2)

#### fit model ####
fitDetNeospora <- function(dir.in,
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
                           K){


  SI5yrEq <- function(time, state, params){

    with(as.list(c(state, params)),{

      Sc <- Sc2+Sc3+Sc4
      Ic <- Ic2+Ic3+Ic4
      N<- So+Sh+Sc+Io+Ih+Ic
      theta <- K - N

      if(theta<0) {
        theta_neg <- theta
        theta_pos <- 0
      } else {
        theta_neg <- 0
        theta_pos <- theta
      }

      rhoh <- alpha*zeta*(Ic/(Sc + Ic))
      po<- Io/(Io+So)

      dSo <- (1-po)*theta_neg + alpha*(1-betas)*Sc+alpha*(1-betaI)*(1-rhov)*Ic-
        (delta+rhoh+sigma)*So-g*So
      dSh <- (1-p)*theta_pos + g*So-(delta+rhoh+sigma)*Sh-g*Sh
      dSc2 <- g*Sh -(delta+rhoh+sigma)*Sc2-g*Sc2
      dSc3 <- g*Sc2-(delta+rhoh+sigma)*Sc3-g*Sc3
      dSc4 <- g*Sc3-(delta+rhoh+sigma+eps)*Sc4

      dIo <- po*theta_neg + alpha*rhov*(1-betaI)*Ic + (rhoh+sigma)*So-delta*Io-g*Io
      dIh <- p*theta_pos+g*Io+(rhoh+sigma)*Sh-delta*Ih-g*Ih
      dIc2 <- g*Ih+(rhoh+sigma)*Sc2-delta*Ic2-g*Ic2
      dIc3 <- g*Ic2+(rhoh+sigma)*Sc3-delta*Ic3-g*Ic3
      dIc4 <- g*Ic3+(rhoh+sigma)*Sc4-(delta+eps)*Ic4


      return(list(c(dSo, dSh, dSc2, dSc3, dSc4, dIo, dIh, dIc2, dIc3, dIc4)))
    })
  }

  params <- c(
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
    g=g
  )

  # work out the initial states
  So=0.25*K*(1-InitPrev); Sh=0.25*K*(1-InitPrev); Sc2=0.5*(1/3)*K*(1-InitPrev);
  Sc3=0.5*(1/3)*K*(1-InitPrev); Sc4=0.5*(1/3)*K*(1-InitPrev);
  Io=0.25*K*InitPrev; Ih=0.25*K*InitPrev; Ic2=0.5*(1/3)*K*InitPrev; Ic3=0.5*(1/3)*K*InitPrev;
  Ic4=0.5*(1/3)*K*InitPrev

  initial_state <- c(So=So, Sh=Sh, Sc2=Sc2, Sc3=Sc3, Sc4=Sc4,
                     Io=Io, Ih=Ih, Ic2=Ic2, Ic3=Ic3, Ic4=Ic4)


  #sum(initial_state)
  #times <- 0:10

  model <- ode(initial_state, times, SI5yrEq, params)
  #model[1:10,]
  #rowSums(model[,-1])
  res<-data.table(model)
  res[, S:= So+Sh+Sc2+Sc3+Sc4]
  res[, I:=Io+Ih+Ic2+Ic3+Ic4]
  res[, N:=S + I]
  res[, Prev:=I/N]
  res[, rhoh:=alpha*zeta*((Ic2+Ic3+Ic4)/(Sc2+Sc3+Sc4 + Ic2+Ic3+Ic4))]
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
    scale_color_manual(values=c("green", "red", "blue")) 

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
  res<-vector("list", length = nrow(parms))
  for(rn in seq_len(nrow(parms))) {
    res[[rn]] <- fitDetNeospora(dir.in,
                                times=0:20,
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
                                K=parms[rn, "K"])

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
