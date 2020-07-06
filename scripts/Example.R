#########################################################################
### Basic Optimal Design Example - Proof of Concept
### Perry Williams & Mevin Hooten
### Computation time with iMac Pro (2017) with 3.0 GHz, 128 GB Mem, 20 cores:
###               18.6 secs : 15 potential designs,   100,000 MCMC samples
###               2.50 mins : 15 potential designs, 1,000,000 MCMC samples
#########################################################################


###
### Load/Install required packages and functions
###

required.packages=c("ContourFunctions",
                    "doParallel",
                    "laGP",
                    "lattice",
                    "msm",
                    "RColorBrewer",
                    "scatterplot3d")

lapply(required.packages,install.packages,character.only=TRUE)
lapply(required.packages,library,character.only=TRUE)

###
### Functions
###

## GLM - probit pprb
source("probit.pprb.R")

## GLM - probit link
source("probit.reg.mcmc.R")

###
### Dimensions
###

nc=5
nr=5
N=nc*nr

###
### Parameters
###

beta0=rnorm(1)
beta1=rnorm(1)
beta=beta.truth=matrix(c(beta0,beta1),2)

###
### Covariates
###

X.all=cbind(1,rnorm(N,0,3))
mu=rnorm(N,X.all%*%beta,1)

###
### Response
###

Y.all=ifelse(mu>=0,1,0)

###
### Sample data
###

n=10
sample1=sort(sample(1:N,n))
y.orig=Y.all[sample1]
X1=X.all[sample1,]

#########################################################################
### Fit model
#########################################################################

beta.mn=matrix(0,2,1)
beta.var=2.25
n.iter=100000
model.output=probit.reg.mcmc(y=y.orig,
                             X=X1,
                             beta.mn,
                             beta.var,
                             n.iter)
beta.samp=model.output$beta.save

###
### Plot MCMC output
###

par(mfrow=c(2,2))
plot(density(beta.samp[,1]),type='l',xlab='Iteration',ylab="beta0",main="")
abline(v=beta0,col=4)
plot(density(beta.samp[,2]),type='l',xlab='Iteration',ylab="beta1",main="")
abline(v=beta1,col=4)
plot((beta.samp[,1]),type='l',xlab='Iteration',ylab="beta0")
abline(h=beta0,col=4)
plot((beta.samp[,2]),type='l',xlab='Iteration',ylab="beta1")
abline(h=beta1,col=4)

###
### Calculate posterior predictive distribution (ppd) of mu
###

mu.ppd=matrix(rnorm(n.iter*N,mean=X.all%*%t(beta.samp),sd=1),N,n.iter)

###
### Calculate ppd of y
###

y.ppd=ifelse(mu.ppd>=0,1,0)
y.ppd.var1=apply(y.ppd,1,var)

#########################################################################
### Where to sample new data?
#########################################################################

###
### Number of new samples
###

new.samples=1

###
### Unsampled sites
###

unsampled=(1:N)[!(1:N%in%sample1)]

###
### All possible combinations of unsampled sites
###

D=t(combn(unsampled,new.samples))
M=dim(D)[1]

###
### Vector to store scores of each sample
###

score=rep(NA,M)

#########################################################################
### Loop over all designs
#########################################################################

##
## Parallel computation
##

registerDoParallel(cores=detectCores())
score=foreach(d=1:M,.combine='c',.export=c('probit.pprb')) %dopar%{

    cat(round(d/M*100),"% - ")

    ##
    ## y.ppd of new sample
    ##

    p.pred=apply(matrix(y.ppd[D[d,],],new.samples,n.iter),1,mean)
    y.pred=ifelse(p.pred>0.5,1,0)
    X.pred=X.all[D[d,],]

    ##
    ## Update model with posterior preds as new data
    ##

    model.output3=probit.pprb(y.pred,X.pred,beta.samp)

    ##
    ## Calculate new score
    ##

    beta.samp.new=model.output3$beta.save
    p.ppd=pnorm(matrix(X.all%*%t(beta.samp.new), N, n.iter))
    var.p.ppd=p.ppd*(1-p.ppd)
    sum.var.p=apply(var.p.ppd,2,sum)
    score[d]=mean(sum.var.p)
}

tick=Sys.time()
tick-tock

###
### Color palatte
###

SpectralColors=colorRampPalette(brewer.pal(11, "Spectral"))
color=SpectralColors(length(score))

if(new.samples==2){
    summary=cbind(X.all[D[,1],2],X.all[D[,2],2],score)
    summary=summary[order(summary[,3]),]

    ##
    ## View results graphically
    ##

    par(mfrow=c(1,1))
    plot(x=summary[,1],y=summary[,2],col=color,pch=16)
    cf_data(x=summary[,1],y=summary[,2],z=summary[,3],
            bar=TRUE,
            color.palette=topo.colors,
            with_lines=TRUE)
}

if(new.samples==1){
    summary=cbind(X.all[D[,1],2],score)
    summary=summary[order(summary[,2]),]

    ##
    ## View results graphically
    ##

    quartz()
    par(mfrow=c(2,1))

    ##
    ## Plot 1
    ##

    X1=seq(-10,10,0.1)
    X=cbind(1,X1)
    plot(x=summary[,1],
         y=summary[,2],
         xlim=range(X1),
         col=color,
         pch=16,
         xlab="X-value",
         ylab="Score"
         )
    tmp=summary[order(summary[,1]),]
    lines(tmp[,1],tmp[,2],col=2)
    low.score=summary[,1][which(summary[,2]==min(summary[,2]))]
    abline(v=low.score,col=4)

    ##
    ## Plot 2
    ##

    beta=apply(beta.samp,2,quantile,0.5)
    plot(X1,pnorm(X%*%beta),xlim=range(X1),ylim=c(0,1),type='l')
    ind=round(seq(1,n.iter,length.out=300))
    for(i in ind){
        beta=beta.samp[i,]
        lines(X1,pnorm(X%*%beta),xlim=range(X1),ylim=c(0,1),type='l',
              col=rgb(0,0,0,.1))
        abline(v=low.score,col=4)
    }
}


