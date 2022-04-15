cor.diff.test = function(x1, x2, y1, y2) {

  r1 = x1
  r2 = x2
  n1 = y1
  n2 = y2
  fisher = ((0.5*log((1+r1)/(1-r1)))-(0.5*log((1+r2)/(1-r2))))/((1/(n1-3))+(1/(n2-3)))^0.5

  p.value = (2*(1-pnorm(abs(fisher))))

  result= list(
    "p.value.twosided" = as.numeric(p.value),
    "p.value.onesided" = as.numeric(p.value) / 2
  )
  cat(paste(sep="",
          "diffence: p(one-sided)=", format(result$p.value.onesided, digits=3), ", p(two-sided)=", format(result$p.value.twosided, digits=3), "\n"
  ))
  return(result);
}
a=read.table("aa",header=F)
for(i in 1:nrow(a)) {
    s =a[i,]
    cor.diff.test(s$V1,s$V2,s$V3,s$V4)
	
    # do stuff with row
}




