
larvae_counts = read.csv('cleaned_count_data.csv')
culture_names = levels(larvae_counts$Culture)
base_name="_chain4"
folder=paste("survival_plots",base_name,sep="")

load(paste("Daily_survival_samples",base_name,".Robj",sep=""))


plot_posterior_mean = function(culture_data, X){
	culture = culture_data$Culture[1]
	culture_data = culture_data[order(culture_data$Day),]
	#exclude days that you don't want to analyze
	
	data  = culture_data[,c("Total_1","Total_2")]
	vols  = culture_data$Vol_Used
	fact  = culture_data$Density_multiplication_factor
	scale = vols/(5000*fact)
	days_between_samples = diff(culture_data$Day)
	X = t(apply(X,1,'*',c(1,days_between_samples)))
	nDays = dim(data)[1]
#plot results
	fitted = colMeans(t(apply(X,1,cumprod)))
	a=as.data.frame(data)
	a=a/scale
	plot(rep(1:nDays,2),unlist(a),ylim=c(0,max(c(unlist(a),fitted))),main=culture)
	lines(1:nDays,rowMeans(a))
	abline(v=(2:nDays)[fact[2:nDays] != fact[2:nDays-1]])
	points(fitted,col=4)	
}

pdf(sprintf('%s/Fitted_survival_curves.pdf',folder))
for(culture in culture_names){
	# pdf(paste("Fitted_survival_",culture,".pdf"))
	plot_posterior_mean(subset(larvae_counts,Culture == culture), daily_survival[[culture]])
	# dev.off()
}
dev.off()

library(coda)
#test autocorelation
acfs=sapply(daily_survival,function(x) apply(x,2,function(y) autocorr(mcmc(y))[3]))



culture_survival_curves = sapply(daily_survival,function(x) t(apply(x,1,cumprod)))
names(culture_survival_curves) = culture_names

days=lapply(culture_names,function(cult) larvae_counts$Day[larvae_counts$Culture==cult])
names(days) = culture_names

fitted_pops=matrix(NA,nrow=length(culture_names),max(unlist(days)))
rownames(fitted_pops) = culture_names
for(culture in as.character(culture_names)) fitted_pops[culture,days[[culture]]] = colMeans(culture_survival_curves[[culture]])

pdf(sprintf("%s/Survival_plots.pdf",folder))
	#fitted_pops=fitted_pops[substr(rownames(fitted_pops),1,1)!= "Z",]
	fitted_pops_std =fitted_pops/apply(fitted_pops,1,max,na.rm=T)
	fitted_pops_std =fitted_pops/fitted_pops[,1]
	
	plot(NA,NA,xlim=c(1,dim(fitted_pops_std)[2]),ylim=c(0,1),ylab="Percent Surviving",xlab = "Days post-fertalization")
		treatments = rep(1,dim(fitted_pops_std)[1])
		treatments[grep("Low",rownames(fitted_pops_std))]=2
		treatments=factor(substr(rownames(fitted_pops_std),4,6))
		Female=factor(substr(rownames(fitted_pops_std),1,1))
		Male=factor(substr(rownames(fitted_pops_std),2,2))
		names(treatments)=culture_names
		for(culture in culture_names){
			data=na.omit(cbind(1:max(unlist(days)), fitted_pops_std[culture,]))
			lines(data,col= c(1,2)[as.numeric(treatments[culture]=="Low")+1])
		}
		legend('topright',legend=c('Ctl','Low'),col=1:2,lty=1)

	D=c("A","B","C","D")
	N=c("W","X","Y","Z")

	Pop = c(rep("D",length(D)),rep("N",length(N)))
	names(Pop)=c(D,N)
	Pop[substr(rownames(fitted_pops_std),2,2)]

	plot(apply(fitted_pops_std[treatments=="Ctl",],2,mean,na.rm=T),ylim=c(0,1),ylab="Percent Surviving",xlab = "Days post-fertalization")
		points(1:23,apply(fitted_pops_std[treatments=="Low",],2,mean,na.rm=T),col=2)
		legend('topright',legend=c('Ctl','Low'),col=1:2,pch=1)

	plot(apply(fitted_pops_std[Pop[substr(rownames(fitted_pops_std),2,2)]=="N",],2,mean,na.rm=T),ylim=c(0,1),ylab="Percent Surviving",xlab = "Days post-fertalization")
		points(1:23,apply(fitted_pops_std[Pop[substr(rownames(fitted_pops_std),2,2)]=="D",],2,mean,na.rm=T),col=2)
		legend('topright',legend=c('Norway','Denmark'),col=1:2,pch=1)


dev.off()
