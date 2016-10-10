bicfa = function(dta,B){
    dta = t(dta)
    n = nrow(dta)
    k=log(n)
    nbf = ncol(B)
    S = cor(dta)
    Psi = 1-apply(B^2,1,sum)
    Sigma = diag(Psi)+B%*%t(B)
    nbpar = sum(abs(B)>1e-06)
    ldet = determinant(Sigma,logarithm=TRUE)$modulus
    iSigma = ifa(Psi,B)
    dr = n*nbf*log(2*pi)+n*ldet+n*sum(diag(S%*%iSigma$iS))
    bic = dr+k*nbpar
    list(bic=bic)
}