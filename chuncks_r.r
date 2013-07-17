library(rockchalk)

## Squared semi-partial correlation coefficient

unlist(lapply(split(res,res$scaling),function(x) getDeltaRsquare(lm(bprod~S1+S2+S3+log(compratio,10),x))))