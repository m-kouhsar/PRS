
# Tests whether a set of columns exist in a data.frame
#' @param dta A data.frame
#' @param cols A vector of columns to look for
#' @return TRUE if all dta has cols, FALSE otherwise
hasAllColumns <- function (dta, cols) {
  sum(cols %in% colnames(dta)) == length(cols)
}

# Verify that an existing file to append scores to meet the requirements
#' @param fileName A file to verify
#' @param scoreName The name of the score to place in the file
#' @param mergeIDs Variables to use when merging score.
verifyScoreOutputFile <- function (fileName, scoreName, mergeIDs) {
  require(bigreadr, quietly = T)
  if (!file.exists(fileName)) stop('File', fileName, 'does not exist.\n')
  # Inspect that the necessary columns to merge are available in this file (IID and FID)
  tmp <- bigreadr::fread2(fileName, nrows=2)
  cnames <- colnames(tmp)
  if (!hasAllColumns(tmp, mergeIDs)) {
    stop('Necessary ID columns (', paste0(mergeIDs, collapse=', '),') for merging not found in ', fileName, '. Found: ', paste0(cnames, collapse=', '))
  }
  if (scoreName %in% cnames) warning('A column with the designated score name (',scoreName,') already exists in output file ', fileName)
}

# Calculate counts of missingness in the genotype matrix, i.e., in the columns of the matrix.
#' @param genoMat A FBM.code256 bigSNP genotype matrix (eg, <bigSNP>$genotypes)
#' @return A list of counts of missingness per genotype across all individuals
countMissingGenotypes <-  function(genoMat, cores=nb_cores()) {
  require(bigsnpr)
  big_apply(genoMat, a.FUN=function (x, ind) colSums(is.na(x[,ind])), a.combine='c', ncores = cores)
}

# Replace missing genotypes with zero
#' @param genoMat A FBM.code256 bigSNP genotype matrix (eg, <bigSNP>$genotypes)
#' @return A FBM.code256 bigSNP genotype matrix
zeroMissingGenotypes <- function(genoMat) {
  if (class(genoMat) != "FBM.code256") stop('Argument must be a FBM.code256 object. Received ', class(genoMat))
  genoMat$copy(code=c(0,1,2, rep(0, 253)))
}

# Test if x is NA or NULL
isVarNAorNULL <- function (x) {
  all(is.na(x)) | all(is.null(x))
}

# Test if x is numeric
# Haven't found a better way to deal with warnings without doing text based filtering.
#' @param x A variable to test.
isNumeric <- function(x) {
  suppressWarnings(!is.na(as.numeric(x)))
}

# Test if there is only numbers in a variable
isOnlyNumeric <-  function (x) {
  sum(isNumeric(x)) == length(x)
}

# Get the indices of a vector for elements that are numeric
#' @param x A vector
getNumericIndices <- function(x) {
  which(isNumeric(x))
}

# Complement sumstats with missing information.
#' @param sumstats A data.frame with sumstats
#' @param reference A data.frame with reference data to complement sumstats
#' @param colRsidSumstats Name of column with RSID/SNP ID in sumstats
#' @param colRsidRef Name of column with RSID/SNP ID in reference data
#' @param colsKeepReference Vector of columns to merge with sumstats
#' @return A handsome data.frame with merged data of sumstats and reference
complementSumstats <- function(sumstats, reference, colRsidSumstats='SNP', colRsidRef='ID', colsKeepReference=c('CHR','POS')) {
  require(data.table, quietly=T) # merge in this package is supposedly faster than base::merge and works the same
  if (!colRsidSumstats %in% colnames(sumstats)) stop('SNP ID column in sumstats (', colRsidSumstats, ') not found')
  if (!colRsidRef %in% colnames(reference)) stop('SNP ID column in reference data (', colRsidRef, ') not found')
  nrRows <- nrow(sumstats)
  colsRef <- c(colRsidRef, colsKeepReference)
  colsRefIntersect <- intersect(colsRef, colnames(reference)) 
  colsRefMissing <- setdiff(colsRef, colsRefIntersect)
  if (length(colsRefMissing) > 0) stop(paste0('Column(s) ', paste0(colsRefMissing, collapse=','), ' was not found in reference data'))
  res <- data.table::merge.data.table(sumstats, reference[,colsRef], by.x=colRsidSumstats, by.y=colRsidRef, all.x=T)
  if (nrow(res) != nrRows) warning('The merge resulted in ', nrow(res), ' rows but input contained ', nrRows, ' rows')
  res
}

# Get effective sample size from various sources in user input
# Note that arguments effectiveSampleSize or cases and controls will override
# effective sample size if provided as a column in sumstats.
#' @param sumstats A data.frame with sumstats
#' @param effectiveSampleSize Precalculated effective sample size
#' @param cases No of cases for a binary trait
#' @param controls No of controls for a binary trait
#' @param colN Column containing effective sample size in sumstats
#' @return Either integer or a vector of integers
getEffectiveSampleSize <- function (sumstats, traitType,effectiveSampleSize=NULL, cases=NULL, controls=NULL,colN=NULL, colNCases=NULL , colNControls = NULL) {
  if(isVarNAorNULL(traitType)) stop("--trait-type must be provided (b for binary or q for quantitative")
  traitType <- tolower(traitType)
  argsCcNA <- isVarNAorNULL(cases) + isVarNAorNULL(controls)
  colCCNA <- isVarNAorNULL(colNCases) + isVarNAorNULL(colNControls)
  cases <- as.numeric(cases)
  controls <- as.numeric(controls)
  esInSumstats <- ifelse(isVarNAorNULL(colN), F, colN %in% colnames(sumstats))
  
  if(!(traitType %in% c("b","q"))){
    stop("Invalid trait type!")
  }
  if(traitType=="b"){
    if(argsCcNA != 0 & colCCNA != 0) stop("Provide both number of cases and controls (as columns in sumstats or as arguments with --n-cases and --n-controls)")
    
    if(argsCcNA == 0) esOut <- 1/((1/cases) + (1/controls))
    if(colCCNA == 0){
      if(!(colNCases %in% colnames(sumstats))) stop(colNCases," is not in sumsats column names")
      if(!(colNControls %in% colnames(sumstats))) stop(colNControls," is not in sumsats column names")
      esOut <- 1/((1/sumstats[, colNCases]) + (1/sumstats[, colNControls]))
      }
  }
  if(traitType=="q"){
    if (isVarNAorNULL(effectiveSampleSize) && !esInSumstats) 
      stop("Effective sample size has not been provided as an argument and no such column was found in the sumstats (column ", colN, ")")
    if (esInSumstats) esOut <- sumstats[, colN]
    if (!isVarNAorNULL(effectiveSampleSize)) {
      if (argsCcNA < 2) stop('Do not provide both --effective sample size and --n-cases/--n-controls')
      esOut <- as.numeric(effectiveSampleSize)
      if (!is.numeric(esOut)) stop('Effective sample size needs to be numeric, received: ', effectiveSampleSize)
    }
  }
  
  esOut
}

# Rename columns in data.frame
#' @param df A data.frame
#' @param old_names A vector of old column names
#' @param new_names A vector of new column names
#' @return A data.frame with renamed columns
rename_columns <- function(df, old_names, new_names) {
  if (length(old_names) != length(new_names)) {
    stop("The length of the old_names and new_names lists must be the same.")
  }
  colnames(df)[match(old_names, colnames(df))] <- new_names
  return(df)
}
