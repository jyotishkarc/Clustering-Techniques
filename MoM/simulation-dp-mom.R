sim.dpm <- function(k)
{
  mu=matrix(0,k,5,byrow=T)
  gen=sample(1:100000,k)
  for(i in 1:k)
  {
    for(j in 1:5){
      mu[i,j]=gen[i]%%10
      gen[i]=floor(gen[i]/10)}
  }
  data <- matrix(0, k*125, 5)
  for(i in 1:k)
  {
    for(j in 1:100)
    {
      data[(100*(i-1)+j),]=rmvnorm(1,mu[i,],diag(0.01,5))
    }
  }
  min.data <- min(data)
  max.data <- max(data)
  for(i in 1:k)
  {
    for(j in 1:25)
    {
      data[(k*100+25*(i-1)+j),]=c(runif(1,min.data-1,max.data+1),runif(1,min.data-1,max.data+1),runif(1,min.data-1,max.data+1),runif(1,min.data-1,max.data+1),runif(1,min.data-1,max.data+1))
      for(p in 1:5)
      {
        if(data[(k*100+25*(i-1)+j),p]>=min.data & data[(k*100+25*(i-1)+j),p]<=max.data)
        {
          j=j-1
          break
        }
      }
    }
  }
  v1=rep(1:k,each=100)
  v2=rep(k+1,(k*25))
  gt=c(v1,v2)
  return(list(data,gt))
}