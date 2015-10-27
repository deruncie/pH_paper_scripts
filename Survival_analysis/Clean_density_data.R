library(gdata)
larvae_counts = read.xls("Larvae_counts_Experiment_2_cleaned.xls",stringsAsFactors=F)
larvae_counts = larvae_counts[larvae_counts$Date != "15-Apr",]
larvae_counts = larvae_counts[-grep("E_",larvae_counts$Culture),]
larvae_counts = larvae_counts[substr(larvae_counts$Culture,1,1)!="Z",]
larvae_counts$Date=as.character(larvae_counts$Date)
larvae_counts$Culture_volume=NA
larvae_counts = larvae_counts[,c('Day','Date','Time','Culture','good_1','bad_1','Total_1','good_2','bad_2','Total_2','Culture_volume','Vol_Used')]
# larvae_counts$Culture_dilution=NA
# larvae_counts$Percent_larvae_removed=NA

culture_volume=read.xls("Experiment_2_Culture_Volume_by_date",stringsAsFactors=F)
date_water_added=read.xls("Expt_2_Date_water_added.xls",stringsAsFactors=F)
date_water_added[,4]=sub(" ","",date_water_added[,4])

culture_names = unique(larvae_counts$Culture)
exclude_days = c()
total_days = max(larvae_counts$Day)

for(culture in culture_names){

	culture_volume=5000 	# starting culture volume
	culture_dilution=1 		# starting culture dilution
	days = sort(unique(larvae_counts$Day[larvae_counts$Culture==culture]))	# days density measured
	for(day in days){
		data_index=larvae_counts$Culture==culture & larvae_counts$Day==day
		if(sum(data_index) != 1){
			#make sure not multiple measures per day
			# if so, pick the first, but return a comment
			print(larvae_counts[data_index,])
			data_index=(1:dim(larvae_counts)[1])[larvae_counts$Culture==culture & larvae_counts$Day==day][1]
		}

		percent_larvae_removed=0
		date=larvae_counts$Date[data_index]
		if(date==date_water_added[date_water_added[,1]==culture,4]) {
			#sampling for gene expression: 1000ml of larvae removed. water added to 5000. So, culture is diluted by this amount, and then brought back to 5000ml
			culture_dilution = culture_dilution * 5000/(culture_volume - 1000)
			culture_volume = 5000
		}

		#Culture_volume is the daily volume of water during sampling. culture_dilution accounts for changing larval concentration by dilution. Multiplying by this number predicts the concentration had the dilution step (either sampling for gene expression and re-filling, or just refilling) hadn't happened.
		larvae_counts$Culture_volume[data_index] = culture_volume
		larvae_counts$Density_multiplication_factor[data_index] = culture_dilution
		# print(paste(culture_dilution,culture_volume))
		
		#after sampling, subtract sample volume from running total
		if(day != max(days)) culture_volume = culture_volume - 2*larvae_counts$Vol_Used[data_index]
		if(date=="28-Mar") {
			#on this date, water filled back to 5000ml after sampling
			culture_dilution = culture_dilution*(5000/culture_volume)
			culture_volume=5000
		}
	}
}
#to predict original larvae numbers (with no mortality): orig = Vol_used/(5000* Density_multiplication_factor)
larvae_counts$Estimated_count_1 = larvae_counts$Total_1*larvae_counts$Density_multiplication_factor *5000/larvae_counts$Vol_Used
larvae_counts$Estimated_count_2 = larvae_counts$Total_2*larvae_counts$Density_multiplication_factor *5000/larvae_counts$Vol_Used

if(length(exclude_days)>0) larvae_counts = larvae_counts[larvae_counts$Day %in% exclude_days==F,]

write.table(larvae_counts,file = 'cleaned_count_data.csv',sep=',',row.names=F)
