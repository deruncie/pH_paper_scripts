# script to calculate SimREL scores.

require(ontoCAT)
require(parallel)
# require(maptools)
outdir = 'Semantic_Similarity_scores_20150531'
type='SimREL'
try(dir.create(outdir))


Calc.Similarity.Scores = function(gene_to_cat, ontology,type = 'SimREL'){
	require(ontoCAT)
	require(parallel)
	# gene_to_cat is a data.frame with 3 columns: GeneID, Category_description, Category_ID
		# all genes in gene_to_cat will be mapped

	# calculate category frequencies
	categories = unique(gene_to_cat$Category_ID)
	category_statistics = lapply(categories,function(cat) data.frame(Category_ID = cat, 
																	 Length = sum(gene_to_cat$Category_ID == cat), 
																	 Category_description = gene_to_cat$Category_description[gene_to_cat$Category_ID == cat][1]))
	category_statistics = do.call(rbind,category_statistics)
	
	#calculate frequency based on total number of genes linked to terms
	category_statistics$p = category_statistics$Length / length(unique(gene_to_cat$GeneID)) 
	category_statistics$IC = -log10(category_statistics$p)

	category_statistics$Category_ID = as.character(category_statistics$Category_ID)
	category_statistics$Category_description = as.character(category_statistics$Category_description)
	if(cat == 'PC') category_statistics$Category_ID  = paste('GO',category_statistics$Category_ID,sep='_')

	n=dim(category_statistics)[1]
	similatiry_matrix = matrix(0,n,n)
	colnames(similatiry_matrix)=rownames(similatiry_matrix) = category_statistics$Category_description

	category_parents = lapply(1:nrow(category_statistics),function(i) {
		id1 = category_statistics$Category_ID[i]
		return(c(id1,sapply(getAllTermParentsById(ontology,id1),function(x) sub("_",":",getAccession(x),fixed=T))))
	})
	names(category_parents) = category_statistics$Category_ID

	for(i in 1:n){
		print(i)
		id1 = category_statistics$Category_ID[i]
		ancestors1 = category_parents[[id1]]
		id1_IC = category_statistics$IC[i]
		id1_genes = gene_to_cat$GeneID[gene_to_cat$Category_ID==id1]
		scores = rep(0,n)
		for(j in i:n){
			id2 = category_statistics$Category_ID[j]
			id2_genes = gene_to_cat$GeneID[gene_to_cat$Category_ID==id2]
			num_intersecting_genes = length(intersect(id1_genes,id2_genes))

			# if((category_statistics$IC[i]+category_statistics$IC[j])==0) next shouldn't be needed, right?

			ancestors2 = category_parents[[id2]]
			common_ancestors = intersect(ancestors1,ancestors2)
			if(length(common_ancestors)==0)	next
			# check if any common ancestors are not represented by any genes. If they are not, then the score is zero
			common_ancestors_index = match(common_ancestors,category_statistics$Category_ID)
			common_ancestors = data.frame(A= common_ancestors,
										  p=category_statistics$p[common_ancestors_index],
										  IC=category_statistics$IC[common_ancestors_index]
										  )
			common_ancestors = na.omit(common_ancestors)
			if(dim(common_ancestors)[1]==0)	next

			A = common_ancestors
			if(type == 'SimOF'){
				# Olivier Fedrigo's similarity Score
				scores[j] = max(num_intersecting_genes/c(category_statistics$L[i],category_statistics$L[j]))
			}
			if(type == 'SimRES'){
				scores[j] = max(-log10(A$p))
			}
			if(type == 'SimLIN'){
				scores[j] = max(2*A$IC/(category_statistics$IC[i]+category_statistics$IC[j]))
			}
			if(type == 'SimREL'){
				scores[j] = max(2*A$IC*(1-A$p)/(category_statistics$IC[i]+category_statistics$IC[j]))
			}
		}
		similatiry_matrix[i,] = scores
	}
	similatiry_matrix[lower.tri(similatiry_matrix)] = t(similatiry_matrix)[lower.tri(similatiry_matrix)]
	return(similatiry_matrix)
}

load('SPUs_annotation_info.Robj')
for(cat in  c('BP','MF','CC','PC')){
	print(cat)
	spus = names(geneInfo)
	spu_to_cat = read.delim(paste("Gene_to_",cat,".txt",sep=''),stringsAsFactors=F,h=F)
	spu_to_cat = spu_to_cat[spu_to_cat[,1] %in% spus,]
	names(spu_to_cat) = c('GeneID', 'Category_description', 'Category_ID')
	gene_to_cat = spu_to_cat


	if(cat != 'PC') ontology = getOntology("/Volumes/JunkDNA/druncie/Urchin_genomes/annotations/PANTHER/PANTHERGOslim_modRel.obo")
	if(cat == 'PC') ontology = getOntology("/Volumes/JunkDNA/druncie/Urchin_genomes/annotations/PANTHER/Protein_class_relationship.obo")

	similatiry_matrix = Calc.Similarity.Scores(gene_to_cat, ontology, type=type)
	write.table(similatiry_matrix,file=sprintf('%s/%s_%s_scores.dat',outdir,cat,type),sep="\t",quote=F)	
}



