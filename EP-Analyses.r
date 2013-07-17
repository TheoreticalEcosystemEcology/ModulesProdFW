library(car)
Outputs = list.files(path='output', pattern='EPweb')
Simu = list.files(path='output',pattern="^EP[1234567890]+", ignore.case=TRUE)
Results = NULL
for (o in Outputs){
	sim = read.table(paste('output/',o,sep=''))
	Results = rbind(Results,sim)
}

Simul = NULL
for (s in Simu)
{
   sim = read.table(paste('output/',s,sep=''))
   Simul = rbind(Simul,sim)
}

res = Results

MotId = c(paste('S',c(1:5),sep=''),paste('D',c(1:8),sep=''))
colnames(res) = c('w','repl','scaling','L',MotId,'S','b','bm','bprod','bpred','tlav','tlvar','tlmax','tlmed','dsk','odsk','idsk','pp','tp')

colnames(Simul) = c('w','repl','scaling','L',MotId,'S','b','bm','bprod','bpred')
Simul = aggregate(Simul,list(w=Simul$w),mean)
colnames(Simul) = paste(colnames(Simul),'_s',sep='')
colnames(Simul)[1] = 'w'

res = merge(res,Simul)

for(i in c(1:nrow(res))){
	Mo = res[i,MotId]
	res[i,MotId] = Mo / sum(Mo)
}

res = res[-which.max(res$bm),]

res$compratio = log(res$S4/res$S5+1,10)

res$co = res$L / res$S^2

res$bm = log(res$bm,10)
res$bprod = log(res$bprod,10)
res$bpred = log(res$bpred,10)

anova(lm(bm~S1+S2+S3+compratio+S+co+tlav+tlvar+tlmax+dsk+odsk+idsk+pp+tp,res))
summary(lm(bm~S1+S2+S3+compratio+S+co+tlav+tlvar+tlmax+dsk+odsk+idsk+pp+tp,res))
vif(lm(bm~S1+S2+S3+compratio+S+co+tlav+tlvar+tlmax+dsk+odsk+idsk+pp+tp,res))

MOD = lm(bm~S1+S2+S3+compratio+S+co+tlav+tlvar+tlmax+dsk+odsk+idsk+pp+tp+D1+D2+D3+D4+D5+D6+D7+D8,res)
summary(MOD)
anova(MOD)
vif(MOD)

anova(lm(bprod~S1+S2+S3+compratio+S+co+tlav+tlvar+tlmax+dsk+odsk+idsk+pp+tp,res))
summary(lm(bprod~S1+S2+S3+compratio+S+co+tlav+tlvar+tlmax+dsk+odsk+idsk+pp+tp,res))

anova(lm(bpred~S1+S2+S3+compratio+S+co+tlav+tlvar+tlmax+dsk+odsk+idsk+pp+tp,res))
summary(lm(bpred~S1+S2+S3+compratio+S+co+tlav+tlvar+tlmax+dsk+odsk+idsk+pp+tp,res))

sim = res
sim = sim[order(sim$bm),]
sim$prank = rank(sim$bm)
write.table(sim,file='ep.dat', row.names=FALSE,quote=FALSE)
