# Quick and dirty code to clean up the PHESANT data when restricting to men and women

library(data.table)

#root_file <- "ukb11214_final_july_reference_QC_more_phenos_and_corrected."
root_file <- "/gscmnt/gc2802/halllab/lganel/ukbb/mt_cn/phesant/040220/output."

# Need the phenotype summary file as well.
# Use run_summarise_phenotypes_cloud.r if these have not yet been created.

# These are just used to pull out the males and females.
#root_file_males <- "new_and_corrected_PHESANT_phenos/ukb11214_rw_restricted_to_QC_males_more_phenos_and_corrected."
root_file_males <- "/gscmnt/gc2802/halllab/lganel/ukbb/mt_cn/phesant/040220/output.males"
root_file_females <- "/gscmnt/gc2802/halllab/lganel/ukbb/mt_cn/phesant/040220/output.females"
# This is for checking I got the right answer when in 'summarise_phenotypes_output_july copy'

phenotype_info_file <- "/gscmnt/gc2802/halllab/lganel/ukbb/mt_cn/outcome_info_final_round3_w_custom.tsv"
coding_info_file <- "/gscmnt/gc2719/halllab/users/lganel/src/PHESANT/variable-info/data-coding-ordinal-info-nov2019-update.txt"

minimum_bin <- function(df_column, minimum=100) {
	
	if(any(df_column=="", na.rm=TRUE)){
		df_column <- df_column[-which(df_column=="")]
		if(length(df_column) == 0) return(FALSE)
	}
	
	if(sum(is.na(df_column)) == length(df_column)) return(FALSE)
	result <- table(df_column)

	if(length(result) == 1) {
		# This is the case when everyone is in a single bin.
		return(FALSE)
	} else {
		# Otherwise, ask the size of the smallest bin.
		return(min(result) > minimum)
	}
}

keep_cat_ordered <- function(df_column, minimum=5000) {
	
	if(any(df_column=="", na.rm=TRUE)){
		df_column <- df_column[-which(df_column=="")]
		if(length(df_column) == 0) return(FALSE)
	}

	if(sum(!is.na(df_column)) >= 5000) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}

get_raw_cts <- function(cts_variable, file_header, file="/gscmnt/gc2802/halllab/lganel/ukbb/mt_cn/ukb41102.nealesamps.nowithdraw.phesant_head.noempty.tsv", single=TRUE)
{	
	# Get matches in the header file.
	match <- grep(paste0('x', cts_variable, '_'), header)

	if(single) {
		if(length(match)==1) {
			return(match)
		} else {
			return(c())
		}
	} else {
		if(length(match)==1) {
			return(c())
		} else {
			return(match)
		}
	}
}

reassign_all_values <- function(cts_variables, phenofile) {
	where <- phenofile$FieldID %in% cts_variables
	changes <- cbind(phenofile$FieldID[where], phenofile$DATA_CODING[where])
	changes <- changes[!is.na(changes[,2]),]
	# create a list of variables for each coding
	codings <- list()
	for(i in names(table(changes[,2]))) {
		codings[[i]] <- changes[which(changes[,2]==strtoi(i)),1]
	}
	return(codings)
}

get_reassignment <- function(reassignment) {
	matrix(as.integer(unlist(strsplit(strsplit(reassignment, split='\\|')[[1]], split='='))), ncol=2, byrow=TRUE)
}

make_the_changes <- function(reassignment_matrix, cts_variables, data_frame)
{
	for(i in 1:nrow(reassignment_matrix)) {
		matches <- data_frame[cts_variables] == reassignment_matrix[i,1]
		if(length(matches) > 0)
			data_frame[cts_variables][matches] <- reassignment_matrix[i,2]
	}
	return(data_frame)
}

change_values <- function(codings, data_frame, coding_info_file) {
	# Need to read in the encoding and determine the changes to make
	reassignments <- fread(coding_info_file, sep=',', header=TRUE, data.table=FALSE)
	reassignments <- reassignments[reassignments$dataCode %in% names(codings),]
	reassignment_matrices <- sapply(reassignments$reassignments, get_reassignment)
	names(reassignment_matrices) <- reassignments$dataCode
	codings <- codings[order(names(codings))]
	reassignment_matrices <- reassignment_matrices[order(names(reassignment_matrices))]

	for(i in 1:length(reassignment_matrices)) {
		data_frame <- make_the_changes(reassignment_matrices[[i]], as.character(codings[[i]]), data_frame)
	}

	return(data_frame)
}

average_over_cts_multi <- function(cts_variable, data_frame)
{
	where <- grep(cts_variable, names(data_frame))
	column <- rowMeans(data_frame[,where], na.rm=TRUE)
	column[is.nan(column)] <- NA
	return(column)
}

irnt <- function(cts_variable) {
    set.seed(1234) # This is the same as was used by PHESANT - for checking.
    n_cts <- length(which(!is.na(cts_variable)))
    quantile_cts <- (rank(cts_variable, na.last = "keep", ties.method = "random") - 0.5) / n_cts
    # use the above to check, but also use frank for the real thing
    cts_IRNT <- qnorm(quantile_cts)	
    return(cts_IRNT)
}

look_for_logical <- function(column) {
	return("TRUE" %in% column | "FALSE" %in% column)
}

# Let's read in the phenotype information file
phenofile <- fread(phenotype_info_file, sep='\t', header=TRUE)

# Want to pull out variables that end up as cts variables in the both_sex PHESANT file.
# First, let's get the header.
file <- "/gscmnt/gc2802/halllab/lganel/ukbb/mt_cn/ukb41102.nealesamps.nowithdraw.phesant_head.noempty.tsv"
header <- strsplit(system(paste0("head -n 1 ", file), intern=TRUE), split='\t')[[1]]

single_cts_columns <- c()
multi_cts_columns <- c()

for(i in 1:4) {
	pheno_summary <- paste0(root_file, i, "_phenosummary.tsv")
	cts_variables <- system(paste("grep IRNT", pheno_summary, "| cut -f1 -d'\t'"), intern=TRUE)

	# These columns can be written to a new file as is (no need to perform any averaging)
	single_cts_columns <- c(single_cts_columns, unlist(lapply(cts_variables, get_raw_cts, header, single=TRUE)))
	multi_cts_columns <- c(multi_cts_columns, unlist(lapply(cts_variables, get_raw_cts, header, single=FALSE)))
}

# awk out the single cts columns and write them to a file.
# Include the userId
single_cts_columns <- c(1, single_cts_columns)
multi_cts_columns <- c(1, multi_cts_columns)
outfile_single <- "/gscmnt/gc2802/halllab/lganel/ukbb/mt_cn/ukb41102.nealesamps.nowithdraw.phesant_head.noempty_cts_single.tsv"
print(paste0("awk -F $'\t' -v OFS=$'\t' '{print ", paste0("$", single_cts_columns, collapse=","), "}' ", file, " > ", outfile_single))
outfile_multi <- "/gscmnt/gc2802/halllab/lganel/ukbb/mt_cn/ukb41102.nealesamps.nowithdraw.phesant_head.noempty_cts_multi.tsv"
print(paste0("awk -F $'\t' -v OFS=$'\t' '{print ", paste0("$", multi_cts_columns, collapse=","), "}' ", file, " > ", outfile_multi))
