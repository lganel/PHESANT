##
## this file sorts the results by -
## 1) combining results if they were run in parts
## 2) combining the flow counts files into one flow, if they were run in parts
## 3) adding the descriptions of the variables from the variable info file

#install.packages("optparse");
library("optparse")

option_list = list(
	make_option(c("-t", "--test"), action="store_true", default=FALSE, help="run test pheWAS on test data (see test subfolder) [default= %default]"),
  	make_option(c("-r", "--resDir"), type="character", default=NULL, help="resDir option should specify directory where results files should be stored", metavar="character"),
	make_option(c("-b", "--numParts"), type="integer", default=NULL, help="number of phenotype parts (used to parellise)"),
	make_option(c("-v", "--variablelistfile"), type="character", default=NULL, help="variablelistfile file name (should be tab separated)", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
  
if (opt$test==TRUE) {
	
	opt$resDir = '../test/results/';
	opt$variablelistfile = '../test/variable-lists/outcome-info3.txt';
        	
} else {
	
	if (is.null(opt$variablelistfile)){
          print_help(opt_parser)
          stop("variablelistfile argument must be supplied", call.=FALSE)
        }
	if (is.null(opt$resDir)){
          print_help(opt_parser)
          stop("resDir argument must be supplied", call.=FALSE)
        }	
}

source("combineFlowCounts.r")
combineFlowCounts();

source("combineAndSortResults.r")
combineAndSortResults();

# add the name of the variable as listed in the phenotype info file
source("addVariableDescriptions.r");
addVariableDescriptions();

write.table(resultsAll, file=paste(opt$resDir,"results-combined.txt",sep=""), row.names=FALSE, quote=FALSE, sep="\t");


