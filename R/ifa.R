ifa = function(Psi, B) {
    if (class(B) == "numeric")
    B = matrix(B, ncol = 1)
    q = ncol(B)
    Phi = rep(0, length(Psi))
    Phi[abs(Psi) > 1e-05] = 1/Psi[abs(Psi) > 1e-05]
    PhiB = tcrossprod(Phi, rep(1, q))
    PhiB = PhiB * B
    G = diag(q) + t(B) %*% PhiB
    GinvtPhiB = tcrossprod(solve(G), PhiB)
    Phib2 = tcrossprod(PhiB, t(GinvtPhiB))
    iS = diag(Phi) - Phib2
    PhiB2 = crossprod(PhiB, B)
    GinvtPhiB2 = crossprod(solve(G), PhiB2)
    Phib2 = tcrossprod(PhiB, t(GinvtPhiB2))
    iSB = PhiB - Phib2
    return(list(iS = iS, iSB = iSB))
}
