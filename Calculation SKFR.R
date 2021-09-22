skfr.cal<-function(X,k,s,gt)
{
   v=c()
   for(i in 1:30){
      print(i)
      
      #v[i]=sparse.km.fr.1(X,k,s,t(initial(t(X),k)),gt)$NMI
      
      #v[i]=sparse.km.fr.2(X,k,s,gt)$NMI
      
      v[i]=NMI(KMeansSparseCluster(X,k,
                         KMeansSparseCluster.permute(X,k,20)$bestw)[[1]]$Cs,gt)
      
      #v[i]=km.pp(X,k,gt)$NMI
   }
   
   print("HI")
   return(mean(v))
}
