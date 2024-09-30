
options(stringsAsFactors=FALSE)

fill.pop<- function(cat, inits, maxAge) {
  pop.list<- list()
  ind<- 1
  for(i in 1:length(inits)){
    n<- inits[i]
    if(n > 0) {
      for(j in 1:n) {
        pop.list[[ind]]<- data.frame(cat=cat[i],
                                     age=sample(1:maxAge,1, 
                                                prob = dgeom(1:maxAge, 0.2)))
        ind <- ind+1
      }
    }
  }
  pop.list
}

#-----------------------------------------------------------------------
schedule<- function(poplist, maxtime, parms, rhoh) {

  # death is combined with removal rate for the maxAge animals
  death<- function(ii, maxtime, age, parms) {
    et<- NULL
    rate<- parms$delta
    if(age == parms$maxAge)  rate <- rate + parms$eps
    etime<- -log(runif(1))/rate
    if(etime < maxtime) et<- data.frame(ID=ii,type="death",time=etime)
    et
  }

  infected<- function(ii, maxtime, parms) {
    et<- NULL
    rate<- rhoh + parms$sigma
    etime<- -log(runif(1))/rate
    if(etime < maxtime) et<- data.frame(ID=ii,type="infected",time=etime)
    et
  }

  elist<- list()
  N<- length(poplist)
  if(length(N) == 0) return(NULL)
  #df<- do.call('rbind',poplist)
  #nI<- length(df$cat[df$cat %in% c("I")])

  for(i in 1:N) {
    ind<- poplist[[i]]
    switch(ind$cat,
           S = {
             etype<- infected(i, maxtime, parms)
             if(!is.null(etype)) elist[[length(elist)+1]]<- etype
             etype<- death(i, maxtime, ind$age, parms)
             if(!is.null(etype)) elist[[length(elist)+1]]<- etype
           },
           
           I = {
             etype<- death(i,maxtime, ind$age, parms)
             if(!is.null(etype)) elist[[length(elist)+1]]<- etype
           },
           D = {
             etype<- death(i,maxtime, ind$age, parms)
             if(!is.null(etype)) elist[[length(elist)+1]]<- etype
           }
           )
  }

  elist
}
#--------------------------------------------------------------------
advance<- function(poplist, elist) {

  mylistsort <- function(xx) {
    xx[order(sapply(xx,function(x) x$time))]
  }

  elist<- mylistsort(elist)
  kill.list<- NULL
  n <- length(elist)
  # nbirths <- 0
  if(n > 0) {
    for(i in 1:n) {
      etype<- elist[[i]]
      switch(etype$type,
             death = {
               id<- etype$ID
               poplist[[id]]$cat<- "D"
             },
             infected = {
               id<- etype$ID
               if(poplist[[id]]$cat %in% "S") poplist[[id]]$cat<- "I"
             }
             )
    }
  }
  # Update the kill.list - DAVE DID NOT HAVE THIS, CHECK THIS IS RIGHT
  df<- do.call('rbind', poplist)
  kill.list <- df$cat == "D"
  if(!is.null(kill.list)) poplist<- poplist[-seq_along(poplist)[kill.list]]
  poplist
}
#---------------------------------------------------------

pop.census<- function(poplist, category) {
  pop<- do.call('rbind',poplist)
  pop.sum<- table(pop$cat)
  pop.line<- pop.sum[match(category,names(pop.sum))]
  pop.line[is.na(pop.line)]<- 0
  pop.line
}
#--------------------------------------------------------
age.animals<- function(poplist, parms) {
  kill.list<- NULL
  if(length(poplist) > 0) {
    poplist <- lapply(poplist, function(x) {
      x$age <- x$age + 1
      return(x)
    })
    
    df <- do.call('rbind', poplist)
    kill.list <- df$age > parms$maxAge
    if(!is.null(kill.list)) poplist <- poplist[-seq_along(poplist)[kill.list]]
  }
  poplist
}
#-------------------------------------------------------
add.infected<- function(poplist, ageI) {
  n<- length(poplist)
  poplist[[n+1]]<- data.frame(cat="I", age=ageI)
  poplist
}

#-----------------------------------------------------
keepNconstant <- function(poplist, K, p){
  npop<- length(poplist)
  theta <- K - npop
  df<- do.call('rbind', poplist)
  # work out the proportion of offspring that are infected
  po <- length(df$age[df$age == 1 & df$cat == "I"]) / length(df$age[df$age == 1])
  kill.list<- NULL
  
  if(theta > 0) { # if need to add animals
    for(j in seq_len(theta)) poplist[[npop+j]] <- 
        data.frame(cat=sample(c("S", "I"), size = 1, prob = c(1-p, p)), age=2)
  } else {
    if(theta < 0) { # if need to remove animals, rm calves
      clicker <- 0
      for(i in seq_along(poplist)) {
        if(poplist[[i]]$age > 1) {
          next # skip if not a calf
      } else {
        if(clicker == abs(theta)) break  # quit the loop if we reach the target #
        if(poplist[[i]]$cat == "S" & rbinom(1, 1, prob = 1 - po)) {
          poplist[[i]]$cat <- "D" # Flag to remove if is a calf, with prob 1-po if it is S
        } else {
          if(poplist[[i]]$cat == "I" & rbinom(1, 1, prob = po)) { # or with prob po if it is I
            poplist[[i]]$cat <- "D"
          }
        }
      }
        
    }
    }
  }
  if(!is.null(kill.list)) poplist <- poplist[-seq_along(poplist)[kill.list]]
  poplist
}
#-----------------------------------------------------

birth <- function(poplist, parms, maxtime) {
  poplist_temp <- list()
  counter <- 1
  for(i in seq_along(poplist)) {
    if(poplist[[i]]$age < 3) next
    if(poplist[[i]]$cat == "S") {
      rate <- parms$alpha * (1 - parms$betas)
      etime <- -log(runif(1))/rate
      if(etime < maxtime) {
        poplist_temp[[counter]] <- data.frame(cat="S", age=1)
        counter <- counter + 1
      } 
    } else {
      rate <- parms$alpha * (1 - parms$betaI)
      etime <- -log(runif(1))/rate
      if(etime < maxtime) {
        poplist_temp[[counter]] <- 
          data.frame(cat=sample(c("S", "I"), size=1, prob=c(1-parms$rhov, parms$rhov)), 
                     age=1)
        counter <- counter + 1
    }
  }
  }
 poplist_temp
}
#-----------------------------------------------------

update_rhoh <- function(poplist, parms) {
  df<- do.call('rbind', poplist)
  nIc <- sum(df$cat == "I" & df$age > 2)
  nSc <- sum(df$cat == "S" & df$age > 2)
  rhoh <-  parms$alpha*parms$zeta*(nIc/(nSc + nIc))
  return(rhoh)
}
# ageI=the age of the infected animal if intro=1

Neo.ibm<- function(popsize, init.pop, tot.time, intro=0, ageI=2, parms) {
  category<- c("S", "I")
  
  pop<- fill.pop(category, init.pop, parms$maxAge)
  rhoh <- update_rhoh(pop, parms)
  pop.sum<- matrix(0,nrow=tot.time,ncol=length(category))
  pop.sum<- data.frame(pop.sum)
  names(pop.sum)<- category

  if(intro==1) pop<- add.infected(pop, ageI)
  pop.sum[1,]<- pop.census(pop, category)

  for(i in 2:tot.time) {
    if(length(pop) == 0) break
    if(i == intro) pop<- add.infected(pop, ageI)
    event.list<- schedule(pop, 1, parms, rhoh)
    pop<- advance(pop, event.list)
    offspring_pop <- birth(pop, parms, maxtime=1)
    pop <- age.animals(pop, parms)
    # combined new offspring with existing pop after aging
    pop <- c(pop, offspring_pop)
    pop <- keepNconstant(pop, parms$K, parms$p)
    pop.sum[i,] <- pop.census(pop, category)
    
    # Needs to update rhoh with new prevalence
    rhoh <- update_rhoh(pop, parms)
  }
  time<- 1:tot.time
  pop.sum<- cbind(time,pop.sum)
  pop.sum
}

