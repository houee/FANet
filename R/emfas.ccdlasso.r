emfas.ccdlasso = function (dta, nbf, lambda = 0, zeros=NULL,min.err=1e-06,max.iter=100) {
    dta=t(dta)
    m = ncol(dta)
    n = nrow(dta)
    mdta = t(rep(1,n))%*%dta/n
    cdta = dta-rep(1,n)%*%mdta
    vdta = (t(rep(1,n))%*%dta^2/n)-mdta^2
    sddta = sqrt(n/(n-1))*sqrt(vdta)
    cdta = cdta/(rep(1,n)%*%sddta)
    S = t(cdta)%*%cdta/(n-1)
    
    if (nbf == 0) {
        B = NULL
        Psi = rep(1, m)
    }
    
    if (nbf > 0) {
        print(paste("Fitting Factor Analysis Model with", nbf,"factors"))
        
        svddta = svd(cdta/sqrt(nrow(dta)-1),nv=n)
        evalues = (svddta$d[1:nbf])^2
        evectors = svddta$v[,1:nbf]
        
        if (nbf > 1)
        B = evectors[, 1:nbf] %*% diag(sqrt(evalues[1:nbf]))
        if (nbf == 1)
        B = matrix(evectors, nrow = m, ncol = 1) * sqrt(evalues[1])
        
        Psi = 1 - (B^2%*%rep(1,nbf))[,1]
        iter = 0
        stop = FALSE
        
        while (!stop) {
            
            iter = iter + 1
            iS = ifa(Psi,B)
            xiSB = cdta%*%iS$iSB
            Cyz = t(cdta)%*%xiSB/(n-1)
            Czz = t(iS$iSB)%*%Cyz+diag(nbf)-t(B)%*%iS$iSB
            mzeros<-matrix(0,nrow=m,ncol=nbf)
            result<-.C("cycliccd",as.double(S), as.double(Cyz), as.double(Czz),
            as.integer(nrow(Cyz)), as.integer(ncol(Czz)), as.double(Psi), as.integer(mzeros),
            as.double(lambda),as.double(min.err), as.integer(max.iter), retourPsi=rep(as.double(0),m),
            retourB=matrix(as.double(0),nrow=m, ncol=ncol(Cyz)))
            faccd=list(Psi=result$retourPsi,B=result$retourB)
            Bnew = faccd$B
            Psinew = faccd$Psi
            crit = mean(abs(sql(B) - sql(Bnew))) ; print(crit)
            stop = (crit<min.err)|(iter>=max.iter)
            B = Bnew
            Psi = Psinew
            
        }
        
    }
    S = cov2cor(diag(Psi) + B %*% t(B))
    res = list(S = S, B = B, Psi = Psi)
    return(res)
    
}



sql = function(M) {
   M2 = M%*%t(M)
   M2[lower.tri(M2)]
}