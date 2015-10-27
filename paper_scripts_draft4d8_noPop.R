library(GSVA)
library(edgeR)
library(glmnet)
library(reshape)
library(qvalue)
source('modified_Limma_fns.R')

gene_counts = as.matrix(read.delim("Data/Gene_counts.txt",check.names=F))
ercc_counts = as.matrix(read.delim("Data/ERCC_counts.txt",check.names=F))
no_feature_counts = as.matrix(read.delim("Data/No_feature_counts.txt",check.names=F))

sample_info = read.delim("Data/Sample_info_updated.txt")
sample_info$Culture = paste("2_",sample_info$Culture,sep="")
sample_info = sample_info[match(colnames(gene_counts),sample_info$Culture),]
sample_info[,11:19] = scale(sample_info[,11:19])

design = model.matrix(~Treatment + Male + Female + Male:Treatment + Female:Treatment + Male:Female,data = sample_info)
test_coefs = list(
				Treatment   = which(grepl('Treatment',colnames(design))),
				Male        = which(grepl('Male',colnames(design))),
				Female      = which(grepl('Female',colnames(design))),
				Male_Trt    = grep("TreatmentLow:Male",colnames(design)),
				Female_Trt  = grep("TreatmentLow:Female",colnames(design)),
				Male_Female = which(grepl("Male",colnames(design)) & grepl("Female",colnames(design)))
			  )

trait_tests = function(Y,design,test_coefs,do_voom=F,lib.size,norm.factors,sample.weights = NULL){
	# do lmFit. Run voom if specified
	if(do_voom){ 
		# incoming data are RNAseq counts.
		Y = voomWithQualityWeights(Y,design=design,lib.size = lib.size*norm.factors,normalize.method='none',plot=F)
		f = lmFit_mod(Y,design,test_coefs = test_coefs)
		prior_trend = T
	} else{
		# already normalized and variance stabilized
		if(!is(Y,'list')) {
			# Y is a matrix.
			if(length(dim(Y)) != 2) {
				# Y is a vector
				Y = list(M = matrix(Y,nr=1))
			} else{
				Y = list(M = Y)
			}
		}
		if(is.null(Y$weights)) Y$weights = array(1,dim = dim(Y$M))
		if(!is.null(sample.weights)){
			Y$weights = t(apply(Y$weights,1,'*',sample.weights))
		}		
		f = lmFit_mod(Y,design,test_coefs = test_coefs)
		prior_trend = F
	}
	trait_names = rownames(Y$M)

	f_eBayes = eBayes_mod(f,trend=prior_trend)
	Fs = f_eBayes$Fs
	p_vals = f_eBayes$Fs.p.values
	rownames(p_vals) = trait_names


	q_vals = apply(p_vals,2,function(p){
		q = tryCatch(qvalue(p,pi0.method='bootstrap'),error=function(e) p.adjust(p,method='BH'))
		if(class(q) == 'numeric') {
			q = list(qvalue = q,
					method = 'BH')
		} else{
			q$method = 'qvalue'
		}
		return(q)
	})

	# calculate percent variance explained by factors
	total_var = apply(Y$M,1,var,na.rm=T)
	fitted_var = apply(f_eBayes$coefficients %*% t(design),1,var)
	coef_SS = t(apply(Fs*f_eBayes$s2.post,1,'*',sapply(test_coefs,length))) #
	coef_var = coef_SS/rowSums(Y$weights)
	perc_total = coef_var / total_var
	perc_fitted = coef_var/fitted_var
	rownames(perc_total) = rownames(perc_fitted) = trait_names

	pi0s = sapply(q_vals,function(x) x$pi0)
	qval_methods = sapply(q_vals,function(x) x$method)
	q_vals = do.call(cbind,lapply(q_vals,function(x) x$qvalue))
	colnames(q_vals) = paste(names(test_coefs),qval_methods,sep='_')
	rownames(q_vals) = trait_names

	results = list(	eB_results = f_eBayes,
					p_vals = p_vals,
					q_vals = q_vals,
					pi0s = pi0s,
					Fs=Fs,
					qval_methods = qval_methods,
					perc_total = perc_total,
					perc_fitted = perc_fitted
				  )
	if(do_voom) {
		results$sample.weights = Y$sample.weights
		results$logCPM = Y$E
	}

	return(results)
}

# first, do gene expression tests
lib.size = apply(rbind(gene_counts,ercc_counts,no_feature_counts),2,sum)
norm.factors = calcNormFactors(gene_counts,method='TMM',lib.size = lib.size)

low_genes = rowSums(gene_counts) < 10
gene_results = trait_tests(Y=gene_counts[!low_genes,],design=design,test_coefs = test_coefs,do_voom=T,lib.size,norm.factors)
 

# now, do GSVA:
load('Data/urchinCats.Robj')
urchinCats = read.delim('Data/urchin_cat_list.txt',stringsAsFactors=F)
urchinCats_sub = urchinCats #[urchinCats$Cat!='Hand',]
cats = unique(urchinCats_sub$Category)
catList = lapply(cats,function(cat) urchinCats_sub$V1[urchinCats_sub$Category==cat])
names(catList) = cats

# do GSVA on logCPM values with min cat size=10 and max cat size=1000
cat_enrichPANTHER = gsva(expr = gene_results$logCPM, gset.idx.list=catList, 
							method='gsva',rnaseq=F,abs.ranking=F,min.sz=10,max.sz=1000,
							no.bootstraps=0,parallel.sz=1,mx.diff=T,tau=1,kernel=T,verbose=T)
catTypes = t(sapply(rownames(cat_enrichPANTHER$es),function(x) strsplit(x,':')[[1]]))
M = data.frame(cat_enrichPANTHER$es,check.names=F)
# rescale to unit variance
M = apply(M,2,'/',apply(M,1,sd))
Ys = lapply(unique(catTypes[,1]),function(cat){
	index = catTypes[,1] == cat
	return(list(M = M[index,],probes = catTypes[index,2]))
	})
names(Ys) = unique(catTypes[,1])

gsva_results = lapply(names(Ys),function(class){
	trait_tests(Ys[[class]],design=design,test_coefs = test_coefs,do_voom=F,lib.size,norm.factors,sample.weights = gene_results$sample.weights)
})
names(gsva_results) = names(Ys)


# now, trait results
#Growth rate
GRs = read.delim("Data/Sam_GR.txt",h=F,stringsAsFactors=F)
GRs[,1] = paste('2_',GRs[,1],sep='')
GRs = GRs[match(sample_info$Culture,GRs[,1]),]
GR_results = trait_tests(GRs[,2],design=design,test_coefs = test_coefs,do_voom=F) # need to have at least 2 rows of data. So repeat it 2x

# Survival
cum_survival = as.matrix(read.delim('Data/Daily_counts_summary_chain5.txt'))
cum_survival = cum_survival[,match(sample_info$Culture,paste('2_',colnames(cum_survival),sep=''))]
cum_survival = t(cum_survival)
days = which(colSums(is.na(cum_survival)) == 0)
days = days[days >= 4]
cum_survival = cum_survival[,days]/cum_survival[,1]
colnames(cum_survival) = paste('cum_survival_',days,sep='')

pheno_traits = data.frame(
						Growth_Rate = GRs[,2], 
						ASIN_cum_survival_6_9 = asin(sqrt(cum_survival[,'cum_survival_9']/cum_survival[,'cum_survival_5'])), 
						ASIN_cum_survival_9_21 = asin(sqrt(cum_survival[,'cum_survival_21']/cum_survival[,'cum_survival_9']))
					)

survival_results = trait_tests(t(pheno_traits[,-1]),design=design,test_coefs = test_coefs,do_voom=F)


# Lasso multiple regressions
pH = as.numeric(factor(sample_info$Treatment))-1

lasso_fits = list()
# first GSVA traits
for(cat_type in unique(catTypes[,1])){
	# collect and standardize predictors
	X = t(cat_enrichPANTHER$es[grep(paste(cat_type,':',sep=''),rownames(cat_enrichPANTHER$es)),])

	# add pH as a predictor
	X = cbind(pH=as.numeric(factor(sample_info$Treatment))-1,X)

	# don't penalize the coefficient of pH
	penalty=c(0,rep(1,dim(X)[2]-1))

	lasso_fits[[cat_type]] = apply(pheno_traits,2,function(y) {
		a=cv.glmnet(X,y,nfolds=length(y),family='gaussian',alpha=1,nlambda=100,lambda.min.ratio = 0.001,penalty.factor=penalty,grouped=FALSE)
		})
}
# Then gene expression
X = cbind(pH,t(gene_counts[!low_genes,]))
penalty=c(0,rep(1,dim(X)[2]-1))
lasso_fits[['GE']] = apply(pheno_traits,2,function(y) {
	a=cv.glmnet(X,y,nfolds=length(y),family='gaussian',alpha=1,nlambda=100,lambda.min.ratio = 0.001,penalty.factor=penalty,grouped=FALSE)
	})


save(sample_info,GRs,days,catList,design,test_coefs,
	gene_results,gsva_results,GR_results,survival_results,
	lasso_fits,
	file='draft_4d8_noPop_trait_testing.RData')
