probit.reg.mcmc=function(y,X,beta.mn,beta.var,n.iter){

    ##
    ##  Subroutines
    ##

    rtn=function(n,mu,sig2,low,high){
        flow=pnorm(low,mu,sqrt(sig2))
        fhigh=pnorm(high,mu,sqrt(sig2))
        u=runif(n)
        tmp=flow+u*(fhigh-flow)
        x=qnorm(tmp,mu,sqrt(sig2))
        x
    }

    ##
    ##  Preliminary Variables
    ##

    X=as.matrix(X)
    y=as.vector(y)
    n=length(y)
    l=dim(X)[2]

    beta.save=matrix(0,n.iter,l)
    mu.save=matrix(0,n.iter,n)
    p.save=matrix(0,n.iter,n)

    y1=(y==1)
    y0=(y==0)

    ##
    ##  Hyperparameters and Starting Values
    ##

    Sig.beta.inv=solve(beta.var*diag(l))
    beta=beta.mn

    ##
    ##  Gibbs Loop
    ##

    for(k in 1:n.iter){
        if(k%%1000==0) cat(k," ")

        ##
        ## Sample mu
        ##

        mu=matrix(0,n,1)
        mu1=rtnorm(sum(y),(X%*%beta)[y1],1,0,Inf)
        mu0=rtnorm(sum(1-y),(X%*%beta)[y0],1,-Inf,0)
        mu[y1,]=mu1
        mu[y0,]=mu0

        ##
        ## Sample beta
        ##

        tmp.chol=chol(t(X)%*%X+Sig.beta.inv)
        beta=backsolve(tmp.chol,backsolve(
                                    tmp.chol,t(X)%*%mu+
                                             Sig.beta.inv%*%
                                             beta.mn,transpose=TRUE)+
                                rnorm(l))

        ##
        ## Calculate p
        ##

        p=pnorm(X%*%beta)

        ##
        ## Save Samples
        ##

        beta.save[k,]=beta
        mu.save[k,]=mu
        p.save[k,]=p


    }
    cat("\n")

    ##
    ##  Write output
    ##

    list(y=y,
         X=X,
         n.iter=n.iter,
         mu.save=mu.save,
         beta.save=beta.save,
         p.save=p.save
         )

}
