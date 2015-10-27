rTruncGamma = function(n,trunc_point_low,trunc_point_high,shape,rate){
	trunc_p_low = pgamma(trunc_point_low,shape,rate)
	trunc_p_high = pgamma(trunc_point_high,shape,rate)
	inv_sample = runif(n,trunc_p_low,trunc_p_high)
	sample = qgamma(inv_sample,shape,rate)
	return(sample)
}

culture_survival=function(culture_data,priors,days_no_mortality=c(),plotPrior=F,nIter=10000,thin=100,burnin=1000){
	# culture = culture to analyze
	# days_no_mortality - mortality is assumed to be zero from this day to the next
	# model:
			# 1) Y_ij \sim Pois(X_i*vol_i/(original_volume * dilution_factor)). 
				# X_i is the fitted number of larvae (had larvae not been removed for gene expression sampling, or other reasons).
				# The current volume of the culture is not relevant. The dilution factor accounts for when water was added.
				# Ultimately, X_i is based on the larva concentration. 
			# 2) X_i = X_0 * prod_{j=1}^{i}(M_j). M_i are survival rates from day j-1 to day j.
			# 3) M_i \sim Ga(a,b)|M_i <= 1. Survival must be between 0 and 1.
	
	prior_start=priors$prior_start
	prior_l0_shape = priors$prior_l0_shape
	prior_l0_rate = prior_l0_shape/priors$prior_start
	prior_M_shape = priors$prior_M_shape
	prior_M_rate = priors$prior_M_rate
	truncation_point_low = priors$truncation_point_low
	if(plotPrior){
		hist(rgamma(1e6,prior_l0_shape,prior_l0_rate))
		hist(rTruncGamma(1e6,truncation_point_low,1,prior_M_shape,prior_M_rate))
		return()
	}

	culture = culture_data$Culture[1]
	culture_data = culture_data[order(culture_data$Day),]
	#exclude days that you don't want to analyze

	data  = culture_data[,c("Total_1","Total_2")]
	vols  = culture_data$Vol_Used
	fact  = culture_data$Density_multiplication_factor
	scale = vols/(5000*fact)
	days_between_samples = diff(culture_data$Day)
	
	
	nDays = dim(data)[1]
	Ms_to_sample = 2:nDays
	Ms_to_sample = Ms_to_sample[Ms_to_sample %in% (days_no_mortality +1) == F]
	
	

	l0 = rgamma(1,10, 10/prior_start)
	# Ms = c(l0,rTruncGamma(nDays-1,truncation_point_low,1,prior_M_shape,prior_M_rate))
	Ms = c(l0,rep(1,nDays-1))
	X=c()

	
	cumsum_larvae = cumsum(rowSums(data)[nDays:1])
	cumsum_larvae = cumsum_larvae[length(cumsum_larvae):1]
		
		
	for(j in 1:(nIter*thin+burnin)){
		if(j %% 1000 == 0) {
			# print(j)
			# points(cumprod(Ms))
		}
		scaled_product = scale[nDays]
		for(i in nDays:2){
			if(i %in% Ms_to_sample){
				shape = cumsum_larvae[i]+prior_M_shape
				rate = 2*prod(Ms[1:(i-1)])*scaled_product+prior_M_rate		
				Ms[i] = rTruncGamma(1,truncation_point_low,1,shape,rate)	
				if(Ms[i] < truncation_point_low) recover()
				# if(Ms[i] == 0) recover()
			}
			scaled_product = scale[i-1] + Ms[i]*scaled_product
		}
		Ms[1] = rgamma(1,shape = cumsum_larvae[1] + prior_l0_shape , rate = 2*scaled_product + prior_l0_rate)	
		if(j>burnin && j %% thin == 0) X=rbind(X,Ms)
	}
	
	X = t(apply(X,1,function(x) x^1/c(1,days_between_samples)))

	return(X)
}
