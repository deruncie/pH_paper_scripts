#Bubble plots of GO categories require a measure of distance between the terms. I use the SimRel distance from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1559652/, which takes account of the term tree as well as the term sizes.

calc_SimREL = function(Gene_to_cat, ontology){
#This function calculates a matrix of simRel distances for the ontology in your species. It takes into account both the tree topology and the number of genes in the term 
#	Gene_to_cat: dataframe listing the genes in each category (Category ID in column 3, Category Name in column 2, GeneName in column 1)
#	ontology: imported ontology tree from the ontoCAT package.
	tryCatch(library(ontoCAT),finally = return())

	#calculate the number of genes in each category
	cat_stats = tapply(1:dim(Gene_to_cat)[1],Gene_to_cat[,3],function(x) data.frame(Length = length(x),term = Gene_to_cat[x[1],2]))
	r = c()
	for(go in names(cat_stats)) r=rbind(r,data.frame(GO=go,cat_stats[[go]]))
	cat_stats = r
	cat_stats$GO = as.character(cat_stats$GO)
	cat_stats$term = as.character(cat_stats$term)
	
	#calculate frequency of each term (p), assuming all genes are part of the root	
	cat_stats$p = cat_stats$Length/length(unique(Gene_to_cat[,1]))	
	#calculate information content (IC)
	cat_stats$IC = -log10(cat_stats$p)
	
	
	#calculate SimREL distance for each pair of categories	
	n=dim(cat_stats)[1]
	simrel = matrix(0,n,n)
	colnames(simrel)=rownames(simrel) = cat_stats$term
	print(n)
	for(i in 1:n){
		print(i)
		id1 = cat_stats$GO[i]
		ancestors1 = c(id1,sapply(getAllTermParentsById(ontology,id1),function(x) sub("_",":",getAccession(x),fixed=T)))
		id1_IC = cat_stats$IC[i]
		id1_genes = Gene_to_cat[Gene_to_cat[,3]==id1,1]
		for(j in i:n){
			id2 = cat_stats$GO[j]
			id2_genes = Gene_to_cat[Gene_to_cat[,3]==id2,1]
			
			if((cat_stats$IC[i]+cat_stats$IC[j])==0) next
			ancestors2 = c(id2,sapply(getAllTermParentsById(ontology,id2),function(x) sub("_",":",getAccession(x),fixed=T)))
			common_ancestors = intersect(ancestors1,ancestors2)
			if(length(common_ancestors)==0)	next
			common_ancestors = data.frame(A= common_ancestors,p=cat_stats$p[match(common_ancestors,cat_stats$GO)],IC=cat_stats$IC[match(common_ancestors,cat_stats$GO)])
			common_ancestors = na.omit(common_ancestors)
			if(dim(common_ancestors)[1]==0)	next
			A = common_ancestors	
			simrel[i,j] = simrel[j,i] = max(2*A$IC/(cat_stats$IC[i]+cat_stats$IC[j])*(1-A$p))
		}
	}
	return(simrel)	
}

categories_MDS = function(d, method='cmdscale', n_rand_starts=1,k=2,niter=1e3,tol=1e-10,trace=F,p=2,ncores=1){
#This function use multidimensional scaling to display the category distance matrix in 2D. It finds an optimal projection using the method provided. Since most MDS methods only find a local optimum, this function can try multiple random starting values.
#	d: distance matrix
#	method: MDS method. 'sammon': generally my preference. Uses sammon function, tries it n_rand_starts times and chooses the version with the lowest stress. 'cmdscale': uses cmdscale function. 'isoMDS': uses isoMDS function. Also tries it n_rand_starts times and chooses the version with the minimum stress. 'PCA': uses principle components analysis. Returns the largest 2 eigenvectors of d.
#	n_rand_starts: If method finds a local optimum, how many random starting conditions should be tested?
#	k: number of dimensions to project d onto.
#	niter: maximum iterations per starting condition
#	tol: convergence criteria
#	trace: should optimization function provide intermediate output?
#	p: see isoMDS
	library(parallel)
	library(MASS)
	if(method=='sammon'){
		n = dim(d)[1]
		d=as.dist(d)
		mdss = mclapply(1:n_rand_starts,function(i) {
			if(i == 1){
				y <-cmdscale(d,k)
			} else{
				y= matrix(rnorm(n*2),nc=2)
			}
			mds = try(sammon(d,y=y,k=k,niter=niter,tol=tol,trace=trace),silent=F)
			if(class(mds)=='try-error') mds = list(stress=Inf)			
			return(mds)
		},mc.cores=ncores)
		stress = sapply(mdss,function(x) x$stress)
		mds = mdss[[match(min(stress),stress)]][[1]]
	}
	if(method=='cmdscale'){
		mds = cmdscale(d,k)
	}
	if(method=='isoMDS'){
		n = dim(d)[1]
		d=as.dist(d)
		mdss = mclapply(1:n_rand_starts,function(i) {
			if(i == 1){
				y <-cmdscale(d,k)
			} else{
				y= matrix(rnorm(n*2),nc=2)
			}
			mds = try(isoMDS(d,y=y,k=k,maxit=niter,tol=tol,trace=trace,p=p),silent=F)
			if(class(mds)=='try-error') mds = list(stress=Inf)			
			return(mds)
		},mc.cores=ncores)
		stress = sapply(mdss,function(x) x$stress)
		mds = mdss[[match(min(stress),stress)]][[1]]
	}
	if(method=='PCA'){
		mds = svd(d)$u[,1:2]
		rownames(mds) = rownames(d)
	}
	return(mds)
}

bubbleplot = function(categories, size, coordinates, circle_scale = 0.1, max_size=NULL, color_scale, 
	circle_colors, border_colors = 'grey50', borderwd = 1.3, line_colors = 1, lwd = .3,
	text_colors = 1,text_cex = .6,text_font=1,xlim= 4*range(coordinates[,1]),ylim = xlim ,
	interactive_labels=F,wrap_text=T,add=F,label_method=1,drawlegend=T,
	main = NULL
	){
# This function draws the bubble plots. The main challenge is drawing the labels in a way that they don't overlap. Two options are given: 1) interactive: it draws all the bubbles and then lets you place labels with your mouse. 2) It draws the labels in a circle around the figure.
#	Note: most parameters can be given as single values, in which case the same value will be applied to all circles..
#	categories: vector of text names of the categories
#	size: number of genes in each category
#	coordinates: matrix of x,y coordinates for each category
#	circle_scale: coefficient multiplies the diameters of the circles for drawing
#	color_scale: vector of colors for the legend
#	circle_colors: vector of colors for filling the circles
#	border_colors: border colors for the circles
#	borderwd: border width factor for the circles
#	line_colors: vector of line colors for connecting text and circles
#	lwd: line width factor for connecting text and circles
#	text_colors: text colors of the labels
#	text_cex: text size factor for the labels
#	text_font: see text
#	xlim,ylim: dimensions of plot. Set ylim=xlim to make plot square the the circle size to have a consistent meaning on the two axes
#	interactive_labels: If True, the script will draw bubble outlines, fill circles one at a time, and wait for a mouse-click to position the next label. If false, uses label_method to automatically apply labels to bubbles.
#	wrap_text: Should long category labels be wrapped across lines? Wrapping happens for names wider than (xlim[2]-xlim[1])/3, but this can be modified below
#	add: Should these circles be added to the existing plot?
#	label_method: 1: labels are drawn in a circle around bubbles. 2: label positions are calculated using the maptools package. Only applies if interactive_labels=F
	library(plotrix)
	if(!add || dev.cur() == 1) plot(NA,NA,xlim=xlim,ylim=ylim,pty='s',main = main)
	
	if(length(circle_colors)==1) circle_colors = rep(circle_colors,length(categories))
	if(length(size)==1) size = rep(size,length(categories))
	if(length(border_colors)==1) border_colors = rep(border_colors,length(categories))
	if(length(borderwd)==1) borderwd = rep(borderwd,length(categories))
	if(length(line_colors)==1) line_colors = rep(line_colors,length(categories))
	if(length(lwd)==1) lwd = rep(lwd,length(categories))
	if(length(text_colors)==1) text_colors = rep(text_colors,length(categories))
	if(length(text_cex)==1) text_cex = rep(text_cex,length(categories))
	if(length(text_font)==1) text_font = rep(text_font,length(categories))
	
	scale = circle_scale*max(c(xlim[2]-xlim[1],ylim[2]-ylim[1]))
	radius = (1/(1*pi)* size)^(1/2)
		
	legend = c(100,1000,10000)
	names(legend) = legend
	legend = sqrt(legend/pi)
	
	if(is.null(max_size)) max_size = max(c(radius,legend))
	
	legend = legend/max_size*scale
	radius=radius/max_size*scale
	
	plot_order = order(-radius)
	interactive_label_order = order(sign(coordinates[,1]),atan(coordinates[,2]/coordinates[,1]))
	
	if(wrap_text){
		max_width = (xlim[2]-xlim[1])/3
		for(i in 1:length(categories)){
			text_width = strwidth(categories[i],cex=text_cex[i])
			if(text_width > max_width){
				elements = strsplit(categories[i],' ')[[1]]
				if(length(elements)==1) next
				categories[i] = elements[1]
				for(j in 2:length(elements)){
					text_lines = strsplit(categories[i],'\n')[[1]]
					if(strwidth(text_lines[length(text_lines)],cex=text_cex[i],font=text_font[i])>max_width){
						categories[i] = paste(categories[i],'\n',elements[j],sep='')
					}else {
						categories[i] = paste(categories[i],' ',elements[j],sep='')
					}
				}
			}			
		}
	}

	if(interactive_labels){
		if(length(dev.list())<2) quartz()
		if(add){
			dev.set(2)
			dev.copy(which=3)
			dev.set(3)
		} else{
			dev.set(3)
			plot(NA,NA,xlim=xlim,ylim=ylim)
		}
		for(i in 1:dim(coordinates)[1]){
			symbols(coordinates[i,1], coordinates[i,2],circle=radius[i],add=T,inches=F,bg= NULL,fg=border_colors[i],lwd=borderwd[i])
		}
		text_coords = matrix(0,length(categories),3)
		for(i in interactive_label_order){
			draw.circle(coordinates[i,1],coordinates[i,2],radius = radius[i],border = border_colors[i],col=circle_colors[i],lwd = borderwd[i])
			print(categories[i])
			text_width = strwidth(categories[i],cex=text_cex[i],font=text_font[i])
			new_pt = locator(n=1)
			adj=c(.5,.5)
			line_x=0
			if(new_pt$x>coordinates[i,1]){
				adj[1]=0
				line_x = new_pt$x
			} else {
				adj[1]=1			
				line_x = new_pt$x# + text_width/2					
			}
			text_coords[i,] = c(x=new_pt$x,y=new_pt$y,adj[1])
			text(text_coords[i,1],text_coords[i,2],categories[i],cex= text_cex[i],col=text_colors[i],adj=adj)
			lines(c(line_x, coordinates[i,1]),c(text_coords[i,2],coordinates[i,2]),lwd=lwd,col=line_colors[i])			
		}
		dev.set(2)
		if(!add) plot(NA,NA,xlim=xlim,ylim=ylim)
		for(i in plot_order){
			draw.circle(coordinates[i,1],coordinates[i,2],radius = radius[i],border = border_colors[i],col=circle_colors[i],lwd = borderwd[i])
			line_x = text_coords[i,1]
			adj=c(text_coords[i,3],.5)
			text(text_coords[i,1],text_coords[i,2],categories[i],cex= text_cex[i],col=text_colors[i],adj=adj,font=text_font[i])
			lines(c(line_x, coordinates[i,1]),c(text_coords[i,2],coordinates[i,2]),lwd=lwd,col=line_colors[i])			
		}
	}else {
		# plot(NA,NA,xlim=xlim,ylim=ylim)
		if(label_method==1){
			library(TeachingDemos)
			distances = max((coordinates[,1]^2 + coordinates[,2]^2)^.5)*1.3
			angles = atan(coordinates[,2]/coordinates[,1])
			text_coords = list()
			text_coords$x = sign(coordinates[,1])*distances*cos(angles)
			text_coords$y = sign(coordinates[,1])*distances*sin(angles)	
			text_height = 1.5*strheight('AAAA',cex=text_cex,font=text_font)
			text_coords$y[text_coords$x>0] = spread.labs(text_coords$y[text_coords$x>0],mindiff = text_height)
			text_coords$y[text_coords$x<0] = spread.labs(text_coords$y[text_coords$x<0],mindiff = text_height)
			text_coords$adj = c(1,0)[(sign(coordinates[,1])==1)+1]		
			text_coords$line_x = text_coords$x
		} else if(label_method==2) {
			library(maptools)
			text_coords=pointLabel(coordinates[,1], coordinates[,2],rownames(coordinates),cex=text_cex,doPlot=F)
			for(i in 1:length(text_coords$x)){
				text_width = strwidth(categories[i],cex=text_cex[i],font=text_font[i])
				line_x=0
				if(text_coords$x[i]>coordinates[i,1]) {
					text_coords$line_x[i] = text_coords$x[i] - text_width/2
				} else {
					text_coords$line_x[i] = text_coords$x[i] + text_width/2
				}
				text_coords$adj[i] = 0.5
			}

		}
		for(i in plot_order){
			draw.circle(coordinates[i,1],coordinates[i,2],radius = radius[i],border = border_colors[i],col=circle_colors[i],lwd = borderwd[i])
			text(text_coords$x[i],text_coords$y[i],categories[i],cex= text_cex[i],col=text_colors[i],adj = text_coords$adj[i],font=text_font[i])
			lines(c(text_coords$line_x[i], coordinates[i,1]),c(text_coords$y[i],coordinates[i,2]),lwd=lwd,col=line_colors[i])
		}
	}
	
	#add colorbar
	height = 0.05*(ylim[2]-ylim[1])
	width = 0.05*(xlim[2]-xlim[1])
	n_colors = length(color_scale)
	y = ylim[1]
	x = (xlim[2]+xlim[1])/2 - height*n_colors/2
	a=sapply(1:n_colors,function(i) rect(x+(i-1)*width,y,x+i*width,y+height,col=color_scale[i],border=0))	
	a=sapply(1:n_colors,function(i) text(x+(i-.5)*width,y+1.5*height,names(color_scale)[i],adj=c(.5,.5),cex=.5))
	
	
	#add legend		
	if(drawlegend){
		x=xlim[1]#+(0:2)*3*max(legend)
		y=ylim[1]
		x=rep(xlim[1],3)
		y=rep(ylim[1],3)
		for(i in 1:length(x)){
			draw.circle(x[i],y[i],radius = legend[i])
		}
		text(x+legend,y,names(legend),srt=90,pos=4)
	}
}

# # #load ontology
# library(ontoCAT)
# cat='BP'
# Gene_to_cat = read.delim(paste("/Volumes/JunkDNA/druncie/Urchin_genomes/annotations/SPUs/SPUs_build7/Gene_to_",cat,".txt",sep=''),stringsAsFactors=F,h=F)
# Gene_to_cat = Gene_to_cat[Gene_to_cat[,1] %in% genes,]
# Gene_to_cat = Gene_to_cat
# if(cat == 'PC') {
# 	ontology = getOntology('/Volumes/JunkDNA/druncie/Urchin_genomes/annotations/PANTHER/Protein_class_relationship.obo')
# 	cat_stats$GO = paste("GO:",cat_stats$GO,sep="")
# } else {
# 	ontology = getOntology("/Volumes/JunkDNA/druncie/Urchin_genomes/annotations/PANTHER/PANTHERGOslim_modRel.obo")
# }

# #calculate simrel distances
# simrel = calc_SimREL(Gene_to_cat, ontology)

# #choose categories
# categories = rownames(simrel)[30:40]
# size = seq(30,1000,length=length(categories))

# #calculate mds distances
# d=1-simrel[categories,categories]
# niter=1e3
# tol=1e-10
# trace=F
# p=2
# n_rand_starts=100
# method='sammon'
# coordinates = categories_MDS(d,method='sammon',n_rand_starts=10,k=2,niter=1e3,tol=1e-10,trace=F,ncores=4)

# # #
# # categories = rownames(coordinates)[1:10]
# # size=seq(100,1000,length=length(categories))
# # coordinates = coordinates[1:10,]

# library(RColorBrewer)
# colMap = brewer.pal(4,'Greens')
# circle_colors = sample(colMap,length(categories),replace=T)
# color_scale = colMap
# names(color_scale)=1:length(colMap)

# bubbleplot(categories, size, coordinates, color_scale, circle_colors,interactive_labels=T,add=F)
# # dev.off()	


