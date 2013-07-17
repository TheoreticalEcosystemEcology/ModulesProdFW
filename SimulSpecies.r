library(rockchalk)

load('webSp.Rdata')
Files = paste('output/web',webSp,'-SPECIES.dat',sep='')
#Outputs = list.files(path='output', pattern='web88-SPECIES')

Results = NULL
for (o in Files){
	sim = read.table(o)
	Results = rbind(Results,sim)
}

res = Results

MotId = paste('S',c(1:5),sep='')
colnames(res) = 
	c('w','repl','scaling','id',MotId,'tl','bm','gen','vul')

for(i in c(1:nrow(res))){
	Mo = res[i,MotId]
	res[i,MotId] = Mo / sum(Mo)
}

res$compratio = res$S4/res$S5

res = subset(res,res$bm >= 0.001)
res$isprod = res$gen == 0

#lapply(split(res,list(res$scaling)),function(x) anova(lm(bm~S1+S2+S3+S4+S5,x)))
#lapply(split(res,list(res$scaling)),function(x) summary(lm(bm~S1+S2+S3+S4+S5,x)))

#lapply(split(res,list(res$scaling)),function(x) anova(lm(bm~S1+S2+S3+compratio,x)))
#lapply(split(res,list(res$scaling)),function(x) summary(lm(bm~S1+S2+S3+compratio,x)))

Fval = function(an){
	Fv = an$`F value`
	Fv = Fv[-length(Fv)]
	print(Fv)
	return(Fv / sum(Fv))
}

lapply(split(res,res$scaling),function(x) anova(lm(bm~S1+S2+S3+S4+S5+gen+vul,x)))
lapply(split(res,res$scaling),function(x) summary(lm(bm~S1+S2+S3+S4+S5+gen+vul,x)))

lapply(split(res,res$scaling),function(x) anova(lm(bm~S1+S2+S3+compratio,x)))
lapply(split(res,res$scaling),function(x) summary(lm(bm~S1+S2+S3+compratio,x)))

plot(res$S4,res$bm,pch=20)