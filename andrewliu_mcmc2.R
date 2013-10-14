######################################
### Stat 221 Hw 3
######################################

##Data
  Impala<-c(15,20,21,23,26)
  Waterbucks<-c(53,57,66,67,72)
  n<-length(Impala)

#functions to evaluate posterior
  
  #marginal posterior of N
  margN<-function(N,y){ 
    if(N>=max(y)){
    lgamma(n*N-sum(y)+1)-lgamma(n*N+2)+sum(lchoose(N,y))-log(N)
    }else{
      print("N must be greater than max(y)")
    }
  }

####MCMC Sampler

  #samples from the marginal posterior of N
  margN.MCMC<-function(init,data,simsize=1000){
    #init is initial value for N
#     init<-50; p=0.5; data<-Impala; simsize=1000
    
    #vector to store values
    Nval<-rep(NA,simsize)
    
    #set initial value
    Nold<-init
    
    #find the lower bound for N
    lower<-max(data)
    
    #generate uniform random variables
    flip<-log(runif(simsize))
    
    for(i in 1:simsize){
    #draws from a negative binomial conditional on previous draw of n
    eold<-(Nold-lower+1/2)/lower; pold<-eold/(1+eold)
    Nnew<-rnbinom(1,size=lower,prob=pold)+lower
    
    #proposal distribution if use new N
    enew<-(Nnew-lower+1/2)/lower; pnew<-enew/(1+enew)
    
    #decide whether to accept or reject
    score<-margN(N=Nnew,y=data)-dnbinom(x=(Nnew-lower),size=lower,prob=pold,log=T)-
          margN(N=Nold,y=data)+dnbinom(x=(Nold-lower),size=lower,prob=pnew,log=T)

    if(score>=flip[i]){Nval[i]<-Nnew}else{
      Nval[i]<-Nold
    }
    
    Nold<-Nnew
    } 
    
    Nval
  }

  #MCMC draws
  Ndraw<-margN.MCMC(init=100,data=Impala,simsize=20000)
 
  #gridding procedure
  margNprob<-sapply(seq(max(Impala),20000,1),margN,y=Impala)
  margNprob<-margNprob-max(margNprob)
  Ndrawgrid<-sample(seq(max(Impala),20000,1),size=10000,replace=T,prob=exp(margNprob))
 
  #comparison of quantiles
  quantile(Ndrawgrid,prob=seq(0.05,0.95,0.1))
  quantile(Ndraw,prob=seq(0.05,0.95,0.1))
