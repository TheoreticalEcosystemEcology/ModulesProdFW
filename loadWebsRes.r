load(ResWebs.Rdata)

res = aggregate(.~w+scaling,res,mean,na.rm=TRUE)

## Scaled
Un = subset(res,res$scaling==0,drop=TRUE)
Mod = lm(bm~S5,Un)
bm.pr = predict.lm(Mod)
Un$bm.pr = bm.pr