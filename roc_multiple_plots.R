library("pROC")
data<-read.table("train",sep="\t",header=T)


aa=list(0)
aaa=list(0)
di=(dim(data)[2])/2
dii=di-1
roc_rose <- plot(roc(data[,1], data[,2]), col = "blue",legacy.axes=TRUE,yaxs="i",xaxs="i")
cols=c("green","red","khaki4","cyan")
pp=3
ss=1
for (x in 1:dii){
  yy=pp+1
  print(pp)
  print(yy)
  print(ss)
  roc_rose <- plot(roc(data[,pp], data[,yy]), col =cols[ss] ,legacy.axes=TRUE, yaxs="i",xaxs="i",add = TRUE)
  pp=pp+2
  ss=ss+1
}
s=0
p=1
for (i in 1:di){
  j=p+1
  print(p)
  print(j)
  print(s)
  aa[i]=auc(data[,p], data[,j])
  p=p+2
  s=s+1
}


name=c("SVM", "Random Forest", "Naive Bayes", "SMO", "j48")
ii=1
ss="AUC="
for (ii in 1:length(name))
{
  aaa[ii]=paste(c(name[ii],ss,round(as.numeric(aa[ii]),2)), collapse = " ")
}

legend("bottomright", legend = aaa, col = c("blue","green","red","khaki4","cyan"), lty = 1,lwd=2)





