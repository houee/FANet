SparseMstep = function(Czz,Cyz,zeros) {
   m = nrow(zeros)
   nbf = ncol(zeros)
   id0 = (1:m)[apply(zeros==0,1,sum)!=0]
   hatB = Cyz%*%solve(Czz)
   if (length(id0)>0) {
      ld = lapply(id0,function(i,z,iCzz,Cyz) {
         res = matrix(0,nrow=nrow(iCzz),ncol=ncol(iCzz))
         res[z[i,]==0,z[i,]==0] = solve(iCzz[z[i,]==0,z[i,]==0])
         res = (iCzz-iCzz%*%res%*%iCzz)%*%Cyz[i,]
         return(res) 
      },z=zeros,iCzz=solve(Czz),Cyz=Cyz) 
      hatB[id0,] = matrix(unlist(ld),nrow=length(id0),ncol=nbf,byrow=TRUE)
   }
   return(hatB)
}
