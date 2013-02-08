#####
# Loading the diet matrices
#####
load("./DIET.Rdata")
# Predators are columns
for(i in c(1:length(DIET))){
	mat = DIET[[i]]
	mat[mat>0] = 1
	EdgeList = NULL
	for(pred in c(1:ncol(mat))){
		for(prey in c(1:nrow(mat))){
			if(mat[prey,pred] == 1){
				EdgeList = rbind(EdgeList,c(pred,prey))
			}
		}
	}
	EdgeList = as.data.frame(EdgeList)
	write.table(EdgeList,file=paste('../webs/EP',i,'.txt',sep=''),row.names=F,col.names=F)
}

#####
# Loading the biomass vectors
#####
load('./Sp_B.Rdata')
for(i in c(1:length(Sp_B))){
	write.table(Sp_B[[i]],file=paste('../webs/EP',i,'.biomass',sep=''),row.names=F,col.names=F)
}