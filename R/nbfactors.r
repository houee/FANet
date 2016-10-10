nbfactors = function(dta, nb.iter=50, sv.sig=0.10) {
    dta = t(dta)
    n = nrow(dta)
    m = ncol(dta)
    H = (1/n)*matrix(1,n,n)
    res = scale(dta)
    uu = fast.svd(t(res),tol=0)
    ndf = n - 1
    dstat = uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
    dstat0 = matrix(0,nrow=nb.iter,ncol=ndf)
    
    for(i in 1:nb.iter){
        res0 = t(apply(t(res), 1, sample, replace=FALSE))
        res0 = t(scale(t(res0)))
        uu0 = fast.svd(res0,tol=0)
        dstat0[i,] = uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
    }
    psv = rep(1,n)
    for(i in 1:ndf){
        psv[i] = mean(dstat0[,i] >= dstat[i])
    }
    for(i in 2:ndf){
        psv[i] = max(psv[(i-1)],psv[i])
    }
    nsv = sum(psv <= sv.sig)
    return(as.numeric(list(n.sv = nsv)))
}