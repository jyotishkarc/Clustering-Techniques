f=function(k){
#k=9
p=2
x=seq(0,1,length.out = k)
y=x
z=outer(x,y)
M=matrix(0,k^2,2)
for(i in 1:k){
  for(j in 1:k){
    M[k*(i-1)+j,]=c(x[i],y[j])
  }
}
X=data_generate(600,M,rep(1/(k^2),k^2),0.05,0.05)
ground_truth=X$label
X=X$data
X=rbind(X,rand(500,2)*50)
return(list(X,ground_truth))
}
