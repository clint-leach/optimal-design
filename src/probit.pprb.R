probit.pprb <- function(y,X,beta.prop){

#
#  PPRB for Optimal Design Problem
#

####
####  Setup Variables
####

n=length(y)
q=dim(beta.prop)[2]
n.mcmc=dim(beta.prop)[1]
beta.save=matrix(0,n.mcmc,q)
beta=beta.prop[10,]
beta.save[1,]=beta
p=pnorm(X%*%beta)

####
####  Begin PPRB Loop
####

for(k in 2:n.mcmc){
#  if(k%%10000==0){cat(k," ")}

  ####
  #### Update Beta
  ####

  idx.star=sample(1:n.mcmc,1)
  beta.star=beta.prop[idx.star,]
  p.star=pnorm(X%*%beta.star)

  mh.1=sum(dbinom(y,1,p.star,log=TRUE))
  mh.2=sum(dbinom(y,1,p,log=TRUE))
  mh=exp(mh.1-mh.2)
  if(mh > runif(1)){
    beta=beta.star
    p=p.star
  }

  beta.save[k,]=beta

}
#cat("\n")

####
####  Write Output
####

list(beta.save=beta.save)

}
