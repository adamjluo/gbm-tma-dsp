# Normalization Methods
HOUSEKEEPERS <- c(
  "C1orf43", "GPI", "OAZ1", "POLR2A", "PSMB2", "RAB7A",
  "SDHA", "SNRPD3", "TBC1D10B", "TPM4", "TUBB", "UBB"
)
setMethod(
  "normalize", "NanoStringGeoMxSet",
  function(object, norm_method = c("q3", "neg", "hk"),
           fromElt = "exprs", toElt,
           housekeepers = HOUSEKEEPERS, ...) {
    norm_method <- match.arg(norm_method)
    switch(norm_method,
           "q3" = {
             q3Norm(object,
                          toElt = toElt, fromElt = fromElt, ...
             )
           },
           "neg" = {
             negNorm(object,
                     toElt = toElt, fromElt = fromElt, ...
             )
           },
           "hk" = {
             hkNorm(object,
                    toElt = toElt, fromElt = fromElt,
                    housekeepers = housekeepers, ...
             )
           },
    )
  }
)

q3Norm <- function(object, toElt, fromElt) {
  qs <- apply(exprs(object), 2, function(x) stats::quantile(x, probs = 0.75, type = 6))
  pData(object)[["q3normFactors"]] <- qs / ngeoMean(qs)
  assayDataElement(object, toElt, validate = TRUE) <- sweep(assayDataElement(object, fromElt), 2L, qs / ngeoMean(qs), FUN = "/")
  return(object)
}

negNorm <- function(object, averageType = "mean", toElt, fromElt) {
  if (!featureType(object) == "Target") {
    stop("Error: Negative Background normalization is for collapsed data set.
        Run function aggregateCounts() to collapse the probes to targets.\n")
  }
  if (is.null(fData(object)[["Module"]])) {
    stop("Error: Module is not specified in the object. Check your GeoMxSet object.\n")
  }
  
  pools <- as.list(unique(fData(object)[["Module"]]))
  pool_neg_norm <- lapply(
    pools,
    function(pool) {
      pool_neg <- fData(object)[which(fData(object)$CodeClass == "Negative" &
                                        fData(object)$Module == pool), "TargetName"]
      if (length(pool_neg) < 1) {
        stop(paste0(
          "Error: No negative could be located for probe pool ",
          pool, ".\n",
        ))
        }
      if (length(pool_neg) > 1) {
        stop(paste0(
          "Error: More than one negative was located for probe pool ",
          pool, ".\n"
        ))
        }
      pool_targets <- fData(object)[which(fData(object)$Module == pool), "TargetName"]
      if (averageType == "mean") {
        AT <- mean(exprs(object[pool_neg,]))
      } else if (averageType == "geomean") {
        AT <- exp(mean(log(exprs(object[pool_neg,]))))
      } else if (averageType == "median") {
        break
      }
      pool_neg_factors <-
        exprs(object[pool_neg,])/AT
      pool_counts <- as.matrix(exprs(object[pool_targets,])) %*%
        diag(1 / pool_neg_factors[1:ncol(pool_neg_factors)])
      norm_list <- list(normFactors = pool_neg_factors, norm_exprs = pool_counts)
      return(norm_list)
      }
    )
  pData(object)[["negnormFactors"]] <- t(do.call(rbind, lapply(pool_neg_norm, "[[", 1)))
  neg_norm_df <- data.frame(do.call(rbind, lapply(pool_neg_norm, "[[", 2)))
  colnames(neg_norm_df) <- colnames(exprs(object))
  neg_norm_df <- neg_norm_df[rownames(exprs(object)), ]
  assayDataElement(object, toElt) <- as.matrix(neg_norm_df)
  return(object)
}

hkNorm <- function(object, averageType, toElt, fromElt, housekeepers) {
  if (!featureType(object) == "Target") {
    stop("Housekeeping normalization is for collapsed data set.
            Run function aggregateCounts() to collapse the probes to targets.\n")
  } else {
    if(analyte(object) == "Protein" & all(housekeepers == HOUSEKEEPERS)){
      housekeepers <- hkNames(object)
    }
    hksubset <- subset(object, subset = TargetName %in% housekeepers)
    hks <- apply(exprs(hksubset), 2, function(x) ngeoMean(x))
    ## Save the normfactors in desired pData element
    if (toElt != "exprs_norm") {
      pData(object)[[paste(toElt, "hkFactors", sep = "_")]] <- hks / ngeoMean(hks)
    } else {
      pData(object)[["hknormFactors"]] <- hks / ngeoMean(hks)
    }
    assayDataElement(object, toElt) <- sweep(assayDataElement(object, fromElt), 2L, hks / ngeoMean(hks), FUN = "/")
    return(object)
  }
}


