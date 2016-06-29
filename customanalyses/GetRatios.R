```{r}
FnInflate     <-  function(balloon, target.dim) {
  #Check number of input dimensions
  dim.balloon     <-  dim(balloon)
  if(length(target.dim)!=length(dim.balloon)) {
    stop("Input array does not have the same dimensionality as indicated in target.dim.")
  }
  #Identify the dimensions along which to inflate the array
  same.dims               <-  dim.balloon == target.dim             #A boolean array indicating whether the arrays correspond in that dimension
  corresponding.dims      <-  (1:length(target.dim))[same.dims]     #Dimensions that correspond between the target dim and the input array
  non.corresponding.dims  <-  (1:length(target.dim))[!same.dims]     #Dimensions that do not correspond between the target dim and the input array
  #Check that the non-corresponding dimensions are all of size 1 in the input array - we don't know what to do if this is not the case
  if(prod(dim(balloon)[non.corresponding.dims] == 1)==FALSE) {
    stop("Input array not suitable for inflation as non-corresponding dimensions are not of size 1.")
  }
  #Create temp array of 1s
  inflated.balloon  <-  array(1,dim=target.dim)
  #And inflate using sweep
  inflated.balloon  <-  sweep(inflated.balloon, corresponding.dims, balloon, "*")
  #Return results
  return(inflated.balloon)
}
library(plyr)
library(reshape2)
rats<-sapply(arrays, "[[", "ratios")

total<-data.frame()
for(i in 1:length(arrays)){
ratios<-rats[[i]]
genes<-arrays[[i]]$genes
input<-arrays[[i]]$input
rat<-ratios


m<-melt(rats[[1]],id.vars=c("day","mouse","gene"))


colnames(m)<-c("day","mouse","gene","ratio")
m$day<-m$day+3
m$input=input[m$gene]
m$gene=genes[m$gene]

data_wide <- dcast(m, gene~  mouse + day)
total<-rbind.fill(total,data_wide)

}
nrow(total)
colnames(total)
```
