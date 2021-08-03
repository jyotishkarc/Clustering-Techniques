sim.gaussian <- function(prob,q){
  
  i <- 1
  a = rep(0,length(prob))
  X = c()
  prob2 = prob
  for(k in 2:length(prob)){
    prob2[k] = prob2[k] + prob2[k-1]
  }
  
  prob2=c(0,prob2)
  
  while(i<=500){
    
    b=runif(1)
    for(j in 1:length(prob)){
      if(prob2[j] <= b && b <= prob2[j+1]) { a[i]=j }
    }
    
    X[i] = rnorm(1,q[a[i]])
    i = i+1
  }
  return(list("gtPi"=a,"data"=X))
}
