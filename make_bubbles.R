
#load datatable linking genes to categories
Categories = read.delim('gene_to_cat_list.txt',stringsAsFactors=F)
cats = unique(Categories$Category)
catList = lapply(cats,function(cat) Categories$V1[Categories$Category==cat])
names(catList) = cats




#Bubble plots of factors

#bubble plots
make_bubbles = function(
	categories,
	scores,
	cat_type='BP',
	main = NULL,
	color_range = NULL,
	nbins = 10,	
	colMap = brewer.pal(n_colors,'PuOr'),
	border_colors = 'grey70',
	border_wd = 0.5,
	text_font = 1,
	text_cex = 0.2,
	text_colors = 1,
	lab_color = 'grey70',
	coordinates = NULL,
	circle_scale = 0.2,
	simrel_path = 'Semantic_Similarity_scores/'
	){

	simrel = as.matrix(read.delim(paste(simrel_path,cat_type,'_simrel_scores.dat',sep=""),h=T,check.names=F))
	rownames(simrel)=colnames(simrel)

	#find how many genes were in each category. 
	categories_fixedNames = sub('/','__',categories,fixed=T)		#sometimes, R replaces a / with a __. Sometimes it doesn't.
	sizes = sapply(categories_fixedNames,function(x) length(catList[[paste(cat_type,x,sep=':')]]))
	categories = sub('__','/',categories,fixed=T)

	if(is.null(coordinates)){
		#calculate mds distances
		d=1-simrel[categories,categories]
		if(length(categories)>2){
			coordinates = categories_MDS(d,method='sammon',n_rand_starts=40,k=2,niter=1e3,tol=1e-10,trace=F,ncores=1)
		} else{
			coordinates = rbind(c(-d[2,1]/2,0),c(d[2,1]/2,0))
		}
	}
	# recover()

	names(categories) = NULL
	categories = sub('\n',' ',categories,fixed=T)
	ylim=c(-1.3,1.3)
	distances = max((coordinates[,1]^2 + coordinates[,2]^2)^.5)*1.3
	angles = atan(coordinates[,2]/coordinates[,1])
	text_coords = list()
	text_coords$x = sign(coordinates[,1])*distances*cos(angles)
	text_coords$y = sign(coordinates[,1])*distances*sin(angles)	
	text_height = diff(ylim)/sum(text_coords$x>0)#0*.15*strheight('AA\nAA',cex=.3,font=text_font)
	text_coords$y[text_coords$x>0] = spread.labs(text_coords$y[text_coords$x>0],mindiff = text_height,min=ylim[1],max=ylim[2])
	text_coords$y[text_coords$x<0] = spread.labs(text_coords$y[text_coords$x<0],mindiff = text_height,min=ylim[1],max=ylim[2])
	text_coords$adj = c(1,0)[(sign(coordinates[,1])==1)+1]		
	# text_coords$line_x = text_coords$x
	text_coords = data.frame(Cat = categories,text_coords)
	rownames(text_coords) = NULL
	text_size = .75*min(c(diff(sort(text_coords$y[text_coords$x>0])),diff(sort(text_coords$y[text_coords$x<0]))))
#
	data = data.frame(x = coordinates[,1],y = coordinates[,2],size = (sizes)^(1/1),score = scores)
	p = ggplot(data,aes(x=x,y,y)) + xlim(-3,3) + ylim(ylim)
	# p = p + geom_point(aes(x=x+.01,y=y-0.01,size = size),color='grey70')
	p = p + geom_point(aes(size = size,fill=score),pch=21,color='grey70') + scale_fill_gradient2(low='purple',mid='white',high='orange') + scale_size_area(max_size=30,breaks=c(100,1000))
	p = p + geom_text(data=text_coords,aes(x=x,y=y,label = Cat,hjust = adj),size=text_size*60)
	# p = p + annotate('text',x=text_coords$x,y=text_coords$y,label = text_coords$Cat,size=1.5,adj = text_coords$adj)
	p = p + annotate('segment',x=text_coords$x,y=text_coords$y,xend = coordinates[,1],yend = coordinates[,2],size=.1)
	p = p + theme_bw() + xlab('') + ylab('')
	print(p)

	return(coordinates)
}

#create two vectors: one of category names, matching the names in the PANTHER, the other of scores for the corresponding factors
categories = NULL
scores = NULL

#the color_range should span the scores. If scores can be positive or negative, I think it's best to ensure that color_range is symmetric around zero
color_range = c(-1,1)
n_colors = 10
try({
factor_cat_coordinates = make_bubbles(		#the coordinates of the bubbles plotted are returned. These can be passed back to another call to make the bubbles plotted in the same place
	categories = categories,
	scores = scores,
	cat_type=cat_type,
	main = NULL,
	color_range = color_range,
	colMap = brewer.pal(n_colors,'PuOr')[n_colors:1],
	nbins = n_colors,
	text_cex = .4,
	circle_scale = 0.15,
	coordinates = NULL
)})
