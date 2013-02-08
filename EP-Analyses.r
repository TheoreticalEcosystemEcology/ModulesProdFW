Outputs = list.files(path='output', pattern='EPweb')
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

res = res[-which.max(res$bm),]

res$compratio = res$S4/res$S5

res$co = res$L / res$S^2

res$bm = log(res$bm,10)
res$bprod = log(res$bprod,10)
res$bpred = log(res$bpred,10)

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

sim = res
sim = sim[order(sim$bm),]
sim$prank = rank(sim$bm)
write.table(sim,file='ep.dat', row.names=FALSE,quote=FALSE)