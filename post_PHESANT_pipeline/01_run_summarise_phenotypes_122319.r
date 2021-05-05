library(data.table)
source("/gscmnt/gc2719/halllab/users/lganel/src/PHESANT/summarise_phenotypes.r")

filename_root <- "/gscmnt/gc2802/halllab/lganel/ukbb/mt_cn/phesant/output."

phenotype_info_file <- "/gscmnt/gc2802/halllab/lganel/ukbb/mt_cn/outcome_info_final_round2_w_custom.tsv"

coding_info_file <- "/gscmnt/gc2719/halllab/users/lganel/src/PHESANT/variable-info/data-coding-ordinal-info-nov2019-update.txt"

get_codings <- function(folder) {
  to_read <- paste0(folder, "/", dir(folder))
  codings_tables <- list()
  for(i in to_read) {
    name <- gsub(paste0(folder, "/coding(.*).tsv"), "\\1", i)
    codings_tables[[name]] <- read.table(i, header=TRUE, sep='\t', quote="", stringsAsFactors=FALSE)
  }
  return(codings_tables)
}
codings_tables <- get_codings('/gscmnt/gc2719/halllab/users/lganel/src/PHESANT/WAS/codings')

for (i in 1:4)
{
	# Both sexes
	hist_filename <- paste0(filename_root, i, "_hist")
	pheno_summary <- paste0(filename_root, i, "_phenosummary.tsv")

	filename <- paste0(filename_root, i)
	tsv_filename <- paste(filename, ".tsv", sep="")
	log_filename <- paste0(filename_root, i)
	log_file <- paste(log_filename, ".log", sep="")
	tsv_data <- read.table(tsv_filename, header=TRUE, sep='\t')
	names(tsv_data)[1] <- "userId"

	outcome_info <- read.table(phenotype_info_file, sep='\t', quote="", comment.char="", header=TRUE)
	print(paste0("both sexes ", i))
	summary_11214 <-  get_hists_and_notes(hist_filename, tsv_data, log_file, outcome_info, codings_tables, samples_for_inclusion=TRUE, check=FALSE, start_column=4)
	write.table(summary_11214, file=pheno_summary, col.names=TRUE, row.names=TRUE, sep='\t', quote=FALSE)
}
