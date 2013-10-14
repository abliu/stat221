######################################
### Stat 221 Hw 3
######################################

## Data
Impala<-c(15,20,21,23,26)
Waterbucks<-c(53,57,66,67,72)

## Functions to evaluate posterior

## marginal posterior of N
margN<-function(N,y){ 
  if(N>=max(y)) {
    lgamma(length(y)*N-sum(y)+1)-lgamma(length(y)*N+2)+sum(lchoose(N,y))-log(N)
  }else{
    -Inf
  }
}

####MCMC Sampler

margN.MCMC<-function(init,data,simsize=1000){
  
  #vectors to store values
  Nval<-rep(NA,simsize)
  thetaval <- rep(NA, simsize)
  
  #set initial value
  Nold<-init[1]
  thetaold <- init[2]
  
  #find the lower bound for N
  lower<-max(data)
  
  #generate uniform random variables
  flip<-log(runif(simsize))
  
  size = lower+15
  for(i in 1:simsize){
    p = runif(1,0.02,0.98)
    #draws from a negative binomial conditional on previous draw of n
    # eold<-(Nold-lower+1/2)/lower; pold<-eold/(1+eold)
    # Nnew<-rnbinom(1,size=lower,prob=pold)+lower
    
    #proposal distribution if use new N
    Nnew <- rnbinom(1, size=size, prob=p) + lower
    thetanew <- rbeta(1, 1+sum(data), 1+Nnew*length(data)-sum(data))
    
    #decide whether to accept or reject
    score<-margN(N=(Nnew-lower),y=data)+dbeta(thetanew, 1+sum(data), 
                                              1+Nnew*length(data)-sum(data), log=T)
    -dnbinom(x=(Nnew-lower),size=size,prob=p,log=T)-
      margN(N=(Nold-lower),y=data)+dnbinom(x=(Nold-lower),size=size,prob=p,log=T)
    -dbeta(thetaold, 1+sum(data), 1+Nold*length(data)-sum(data), log=T)
    
    if(score>=flip[i]){Nval[i]<-Nnew; thetaval[i] <- thetanew}else{
      Nval[i]<-Nold
      thetaval[i]<- thetaold
    }
    
    Nold<-Nnew
    thetaold <- thetanew
  } 
  
  data.frame(nval=Nval, thetaval=thetaval)
}

#MCMC draws
runMCMC <- function(ID, data) {
  if (missing(ID)) { args <- commandArgs(TRUE) ; ID <- as.numeric(args[1]) }
  outdir <- sprintf("out-%i", ID) ; dir.create(outdir)
  p_init = runif(1,0.02,0.98)
  n_init = rnbinom(1, size=max(data), prob=p_init) + max(data)
  results <- margN.MCMC(init=c(n_init, p_init),data=data,simsize=20000)
  save(results, file=sprintf("%s/output.rda", outdir))
}
