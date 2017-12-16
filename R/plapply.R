plapply<-function(X, FUN, ...){
if (exists("mclapply"))
  {o<-mclapply(X, FUN,...)}
else
  {o<-lapply(X, FUN,...)}
o
}
