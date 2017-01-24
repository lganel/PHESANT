
# returns boolean var - whether this field denotes the trait of interest
getIsExposure <- function(varName) {

	idx=which(vl$phenoInfo$FieldID==varName)
        isExposure = vl$phenoInfo$EXPOSURE_PHENOTYPE[idx]
        if (!is.na(isExposure) & isExposure!="") {
		return(TRUE)
    	}
	return(FALSE)

}