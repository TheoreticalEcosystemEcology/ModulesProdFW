library(igraph)

##igraph tests

angle.mean = function(angles){
	X = NULL
	Y = NULL
	for(a in angles){
		X = c(X,cos(a))
		Y = c(Y,sin(a))
	}
	n = length(angles)
	return(atan2(sum(Y)/n,sum(X)/n))
}

layout.foodweb = function(w, tlfunc = min, ...)
{
	# Who is a primary producer ?
	PrimProd = (degree(w,mode='out') == 0)
	# Get trophic level
	fullSP = shortest.paths(w)
	SP = shortest.paths(w,V(w)[ names(PrimProd[!PrimProd]) ],V(w)[ names(PrimProd[PrimProd]) ])
	V(w)[ names(PrimProd[PrimProd]) ]$y = 0
	LE = apply(SP,1, tlfunc, ...)
	for(i in c(1:length(LE))){
		V(w)[ names(LE)[i] ]$y = LE[i]
	}
	# Put species at random x positions
	V(w)$x = sample(seq(from=-1,to=1,length=length(V(w))))
	#
	# Get x positions
	for(rep in c(1:5))
	{
		for(v1 in sample(c(1:length(V(w))))){
			V(w)[v1]$x = mean(V(w)[fullSP[v1,]==1]$x)
		}
		for(uy in unique(V(w)$y)){
			OkNodes = V(w)[V(w)$y == uy]$x
			NOk = length(OkNodes)
			V(w)[V(w)$y == uy]$x = seq(from=-NOk,to=NOk,length=NOk)[rank(OkNodes,ties.method='random')]
		}
	}
	return(cbind(V(w)$x,V(w)$y))
}

layout.foodweb.radial = function(w, tlfunc = min, ...)
{
	cscale = function(v,m=0,M=1)
	{
		v <- v-min(v)
		v <- v/max(v)
		v <- v*(M-m)
		v <- v+m
		return(v)
	}
	## First we get the linear positions
	LAY = layout.foodweb(w, tlfunc, ...)
	## Second we re-order the TL
	LAY[,2] = rev(sort(unique(LAY[,2])))[match(LAY[,2],sort(unique(LAY[,2])))] + 1
	## The, for each TL, we change the X position into an angle
	for(TL in unique(LAY[,2])){
		OK = LAY[,1][LAY[,2] == TL]
		if(length(OK) == 1){
			LAY[,1][LAY[,2] == TL] = 0
		} else {
			NotUse = (2*pi / (length(OK)+1))
			LAY[,1][LAY[,2] == TL] = cscale(LAY[,1][LAY[,2] == TL],-pi,pi-NotUse)
		}
	}
	## And we convert to linear
	X = LAY[,2] * cos(LAY[,1])
	Y = LAY[,2] * sin(LAY[,1])
	return(cbind(X,Y))
}

w = read.table('webs/w27.txt')
w[,1] = w[,1]-1
w[,2] = w[,2]-1
w = graph.data.frame(w)

LAY = layout.foodweb(w,min)
plot(as.undirected(w), layout=LAY, edge.curved=TRUE, vertex.label = NA, vertex.size = 2.5, edge.color = rgb(0.5, 0.5, 0.5, 0.1), vertex.frame.color = NA)

LAYR = layout.foodweb.radial(w,median)
plot(as.undirected(w), layout=LAYR, edge.curved=TRUE, vertex.label = NA, vertex.size = 2.5, edge.color = rgb(0.5, 0.5, 0.5, 0.1), vertex.frame.color = NA, vertex.color = 'black')

pdf(file='FoodWebs.pdf',width=20)
par(mfcol=c(1,3))
for(f in list.files(path='webs/')){
	w = read.table(paste('webs/',f,sep=''))
	w[,1] = w[,1]-1
	w[,2] = w[,2]-1
	w = graph.data.frame(w)
	LAYR = layout.foodweb.radial(w,min)
	plot(as.undirected(w), layout=LAYR, edge.curved=FALSE, vertex.label = NA, vertex.size = 2.5, edge.color = rgb(0.5, 0.5, 0.5, 0.2), vertex.frame.color = NA, vertex.color = 'black')
	LAY = layout.foodweb(w,min)
	plot(as.undirected(w), layout=LAY, edge.curved=FALSE, vertex.label = NA, vertex.size = 2.5, edge.color = rgb(0.5, 0.5, 0.5, 0.2), vertex.frame.color = NA, vertex.color = 'black')
	plot(as.undirected(w), layout=layout.fruchterman.reingold, edge.curved=FALSE, vertex.label = NA, vertex.size = 2.5, edge.color = rgb(0.5, 0.5, 0.5, 0.2), vertex.frame.color = NA, vertex.color = 'black')
}
dev.off()