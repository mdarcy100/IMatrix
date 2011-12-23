##############################################################################
#  M.  D'Arcy - first attempt at making imatrix code somewhat generalizable
#  11/2011
#############################################################################
library(Biobase)
library(yaml)
library(base)
library(impute)
library(rJava)  #for hashmaps
setwd("/Users/mdarcy100/Desktop/MTroester/Fahti/")
source("fsam-log.R")
#########################################################################
# a function - will put in source directory
basicfilter<-function(fulldata,basicrow,firstcol,firstrow){
	fulldata_basicfilter<-fulldata[apply(is.na(fulldata),1,mean)<basicrow,,];
	print(paste("number of rows after very basic filter: >", basicrow*100," missing ", nrow(fulldata_basicfilter)))

	fulldata_colfilter<-fulldata_basicfilter[,apply(is.na(fulldata_basicfilter),2,mean)<firstcol,drop=F];
	print(paste("number of columns after main sample filter: (",firstcol*100,"% missing)", ncol(fulldata_colfilter)))
	
	fulldata_filtered<-fulldata_colfilter[apply(is.na(fulldata_colfilter),1,mean)<firstrow,,];
	print(paste("number of rows after ",firstrow*100,"% missing row filter: ", nrow(fulldata_filtered)))
	
	return(fulldata_filtered)
}
#########################################################################
#  STEP 1: process the raw data 
#############################################################################
#read in initialization file (need yaml to do this)
i_config = yaml.load_file("imatrix.yml") #contains all configurations information for imatrix

#read in file (may be able to put this in the initialization file)
#exp_raw<-read.delim("11.8.11Rawdata_2_10842.txt", header=TRUE,sep="\t",check.names=FALSE, stringsAsFactors=FALSE)
exp_raw<-read.delim("11.8.11Rawdata_2_10842.txt", header=TRUE,sep="\t")

dim(exp_raw)
#[1] 51712   162

#first remove rows columns with too much missing data
# read these values in from the configuration file
#i_config$strictprobe_mis = throw out probes with > this percentage missing (for experiments done on multiple platforms)
#i_config$array_mis = throw out arrays with this percentage missing (generally 70%-80%)
#i_config$probe_mis = throw out probes with this percentage missing (generally 70%-80%)

### there was a LOT of missing data in this dataset!
filtered_exp_raw<-basicfilter(exp_raw,as.numeric(i_config$strictprobe_mis),as.numeric(i_config$array_mis),as.numeric(i_config$probe_mis))
dim(filtered_exp_raw)

####impute missing data before moving on (will average samples later in the process)
filtered_exp_raw.imputed <- impute.knn(filtered_exp_raw,k=10, rng.seed=98765)$data

numcols=i_config$imatrix_col  #how many columns are expcted
medialist<-i_config$media
medialist<-unlist(strsplit(i_config$media, "\\,"))

##########################################################################
#initialize java in R for hashmaps
.jinit()
rmf_mono_hm<-.jnew("java/util/HashMap") ### create a list of rmf monocultures in media
cline_media_hm<-.jnew("java/util/HashMap")####create a list of the media that various experiments were done in 
rmf_position_hm<-.jnew("java/util/HashMap") #keeps the position in rmf master matrix for quick look up

#need to use jrcall instead of jcall since jrcall uses reflection to get types? (someone on internet wrote this)

##########################################################################
#iterate through media list, create a new vector of monocultures of RMF in media
#extract monoculture RMFs and average them for use in imatrix method
mList <- list()

for(i in 1:length(medialist)){ #names of all media used
	print(medialist[i])		
	mList[[medialist[i]]]<-i_config$RMFmono[[medialist[i]]]	#print(mList[[medialist[i]]])
	.jrcall(rmf_mono_hm,"put",medialist[i],i_config$RMFmono[[medialist[i]]]) #all the column names should be listed - there were a few repeats.  that is ok
}

##########################################################################
####create a list of the media that various experiments were done in 
full_cell_medialist<-unlist(strsplit(i_config$CoCulturesMedia,"\\;"))

for(i in 1:length(full_cell_medialist)){
	#print(full_cell_medialist[i])
	cell<-unlist(strsplit(full_cell_medialist[i],"\\:"))
	.jrcall(cline_media_hm,"put",cell[1],cell[2])
	print(cell)	
}
## make sure that this worked
keySet <- .jrcall(cline_media_hm,"keySet")
an_iter<-.jrcall(keySet,"iterator")
aList<-list()
while(.jrcall(an_iter,"hasNext")){
	key<-.jrcall(an_iter,"next")
	print(key)
	value<-.jrcall(cline_media_hm,"get",key)
	print(value)
	}
##########################################################################



#convert to an R list? make sure get same result for both
keySet <- .jrcall(rmf_mono_hm,"keySet")
an_iter<-.jrcall(keySet,"iterator")
aList<-list()
while(.jrcall(an_iter,"hasNext")){
	key<-.jrcall(an_iter,"next")
	print(key)
	aList[[key]]<-.jrcall(rmf_mono_hm,"get",key)
	print(aList[[key]])
	}

#extract out columns of matrix for monoculture rmfs
#average and create new dataset with those columns
#put the the column of the rmf in 
##########################################################################
# average monoculture RMFs
count = 1
rmf_matrix<-c()
print(colnames(filtered_exp_raw))  #check the column names
keySet <- .jrcall(rmf_mono_hm,"keySet")
an_iter<-.jrcall(keySet,"iterator")
while(.jrcall(an_iter,"hasNext")){
	key<-.jrcall(an_iter,"next")
	print(key)
	value<-.jrcall(rmf_mono_hm,"get",key)
	print(value)
	v_value<-as.vector(unlist(strsplit(value,"\\;")))	#print(v_value)
	rmf_cols= which(colnames(filtered_exp_raw) %in% v_value)
	#average the values and create a new column vector with the average
	print(rmf_cols)	
	cor(filtered_exp_raw[,rmf_cols])
	#put the position in the column in a hashmap
	print(key)
	print(count)
	.jrcall(rmf_position_hm,"put",key,as.character(count)) #all the column names should be listed)
	count = count+1
}
#######################################################################################################
#pull out monoculture data from cell lines.  there should be as many as the length of: cline_media_hm
# average all monocultures of the same cell line.  create a list like above.
#Name like: C_RMF:MCF10DCIS(0:1,48h) => C_string:string(0:Y,48h) .  X = 0 & Y = 1 => monoculture of the cell line.  find all the monocultures of the cell line
#######################################################################################################
rmf_columns<-list()
for(i in 1:length(mList)){
	print(mList[[i]])
	#average values
}
#######################################################################################################
# pull out co-cultures from the list.  i'm sure we lost some in the missingness filtering process
# Name like: C_RMF:MCF10DCIS(1:2,48h) => C_string:string(X:Y,48h) where X and Y are the proportions of the co-culture
#######################################################################################################


#######################################################################################################
# create a master list of RMF monoculture, cell line, co-culture with appropriate names.  can use this file
# to start the buess method.
# write out in case something goes wrong
#######################################################################################################


#########################################################################
#  STEP 2: remove low variability probes
#############################################################################


#########################################################################
#  STEP 3: create the IMatrix from the master list
#############################################################################