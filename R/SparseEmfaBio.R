SparseEmfaBio = function (dta, zeros, min.err=1e-03) {
    if (sum(rownames(zeros)==rownames(dta))!=nrow(dta)) {
        stop("Row names of constrains matrix and of the expression dta don't match")
    }
    nbf=ncol(zeros)
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
        crit = 1
        
        while (crit > min.err) {
            iter = iter + 1
            iS = ifa(Psi,B)
            xiSB = cdta%*%iS$iSB
            Cyz = t(cdta)%*%xiSB/(n-1)
            Czz = t(iS$iSB)%*%Cyz+diag(nbf)-t(B)%*%iS$iSB
            Bnew = SparseMstep(Czz,Cyz,zeros)
            Psinew = 1 - (Bnew^2%*%rep(1,nbf))[,1]
            crit = mean((Psi - Psinew)^2) ; print(crit)
            B = Bnew
            Psi = Psinew
        }
        
        sB = scale(t(B), center = FALSE, scale = sqrt(Psi))
        G = solve(diag(nbf) + sB %*% t(sB))
        sB = scale(t(B), center = FALSE, scale = Psi)
        Factors = cdta%*%t(sB)%*%t(G)
    }
    B[zeros==0] = 0
    S = cov2cor(diag(Psi)+B%*%t(B))
    return(S)
}
