rm(list = ls(.GlobalEnv), envir = .GlobalEnv)

library("ggplot2")
library("reshape")
library("scales")
source("functions.r")

############################################################
#mainfile
data_dir <- '../_ipy/Data';
file_names_rep <- c("sp101_controls_negative_clust.csv");
data_rep <- read.csv(file=paste(data_dir,"/",file_names_rep[1],sep=""),sep = ",",header=TRUE);

############################################################
data_temp <- data_rep;
temp_frame <- as.data.frame(matrix(ncol=2, nrow=nrow(data_temp)));
for (i in 1:nrow(data_temp))
{
	if (data_temp[i,]$readct_neg>0)
	{
		R <- data_temp[i,]$repreadct;
		q <- data_temp[i,]$readct_neg/data_temp[i,]$repreadct_neg;
		n <- data_temp[i,]$readct;
		mle_stats <- mle_estimation(R,q,n)
		temp_frame[i,1] = mle_stats[1]
		temp_frame[i,2] = mle_stats[2]		
	} else
	{
		temp_frame[i,1] = data_temp[i,]$readct/data_temp[i,]$repreadct;
		temp_frame[i,2] = 0	
	}

}

colnames(temp_frame, do.NULL = TRUE, prefix = "col")

############################################################
#change for appended column id
colnames(temp_frame) <- c('mle','std')
############################################################

data_all <- data.frame(data_temp,temp_frame)

write.csv(data_all,paste0(data_dir,"/",gsub(".csv","",file_names_rep),'_mle.csv'));

