Outputs = list.files(path='output', pattern='WEB')
Results = NULL
for (o in Outputs){
	sim = read.table(paste('output/',o,sep=''))
	Results = rbind(Results,sim)
}

res = Results

MotId = c(paste('S',c(1:5),sep=''),paste('D',c(1:8),sep=''))
colnames(res) = 
c('w','repl','scaling','L',MotId,'S','b','bm','bprod','bpred')

for(i in c(1:nrow(res))){
	Mo = res[i,MotId]
	res[i,MotId] = Mo / sum(Mo)
}

res$compratio = res$S4/res$S5

res$co = res$L / res$S^2

lapply(split(res,res$scaling),function(x) anova(lm(bm~S1+S2+S3+S4+S5,x)))
lapply(split(res,res$scaling),function(x) summary(lm(bm~S1+S2+S3+S4+S5,x)))

lapply(split(res,res$scaling),function(x) anova(lm(bprod~S1+S2+S3+S4+S5,x)))
lapply(split(res,res$scaling),function(x) summary(lm(bprod~S1+S2+S3+S4+S5,x)))

lapply(split(res,res$scaling),function(x) anova(lm(bpred~S1+S2+S3+S4+S5,x)))
lapply(split(res,res$scaling),function(x) summary(lm(bpred~S1+S2+S3+S4+S5,x)))

## Using compratio

lapply(split(res,res$scaling),function(x) anova(lm(bm~S1+S2+S3+log(compratio,10),x)))
lapply(split(res,res$scaling),function(x) summary(lm(bm~S1+S2+S3+log(compratio,10),x)))

lapply(split(res,res$scaling),function(x) anova(lm(bprod~S1+S2+S3+log(compratio,10),x)))
lapply(split(res,res$scaling),function(x) summary(lm(bprod~S1+S2+S3+log(compratio,10),x)))

lapply(split(res,res$scaling),function(x) anova(lm(bpred~S1+S2+S3+log(compratio,10),x)))
lapply(split(res,res$scaling),function(x) summary(lm(bpred~S1+S2+S3+log(compratio,10),x)))

sim = aggregate(.~w+scaling,res,mean,na.rm=TRUE)
sim = sim[order(sim$bm),]

Scaled = subset(sim,sim$scaling==2)
UnScaled = subset(sim,sim$scaling==0)

Scaled$prank = rank(Scaled$bm)
write.table(Scaled,file='w-scaled.dat', row.names=FALSE,quote=FALSE)

UnScaled$prank = rank(UnScaled$bm)
write.table(UnScaled,file='w-unscaled.dat', row.names=FALSE,quote=FALSE)
