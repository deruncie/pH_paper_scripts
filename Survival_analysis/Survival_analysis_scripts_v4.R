options(stringsAsFactors=F)
larvae_counts = read.csv('cleaned_count_data.csv')
culture_names = unique(larvae_counts$Culture)

base_name="_chain5"
folder=paste("survival_plots",base_name,sep="")
source("Survival_function_v4.R")
dir.create(folder)
setwd(folder)
library(parallel)
ncores=detectCores()

priors = list(  
				prior_l0_shape       = 10,
				prior_start          = 50000,
				prior_M_shape        = 1,
				prior_M_rate         = 1,
				truncation_point_low = 00
				)
# culture_survival('',priors,plotPrior=T)

daily_survival = mclapply(culture_names,function(culture) culture_survival(subset(larvae_counts,Culture == culture),priors,c(),nIter=1000,thin=10,burnin=1000),mc.cores=ncores)
names(daily_survival) = culture_names

total_days = max(larvae_counts$Day)
daily_counts_summary = sapply(culture_names,function(culture) {
	X=rep(NA,total_days); 
	X[larvae_counts$Day[larvae_counts$Culture==culture]] = colMeans(t(apply(daily_survival[[culture]],1,cumprod))); 
	X})

daily_survival_summary = sapply(culture_names,function(culture) {
	X=rep(NA,total_days); 
	X[larvae_counts$Day[larvae_counts$Culture==culture]] = colMeans(daily_survival[[culture]]); 
	X})

setwd('../../Data')
save(daily_survival,file=paste("Daily_survival_samples",base_name,".Robj",sep=""))
write.table(daily_counts_summary,file=paste("Daily_counts_summary",base_name,".txt",sep=""),sep="\t",row.names=F,quote=F,,col.names=T)
write.table(daily_survival_summary,file=paste("Daily_survival_summary",base_name,".txt",sep=""),sep="\t",row.names=F,quote=F,col.names=T)
# write.table(larvae_counts,file=paste("Augmented_larvae_counts",base_name,".txt",sep=""),sep="\t",row.names=T,quote=F)
setwd("..")
