# Load relevant packages
packages <- c("GSVA", "GSEABase")
lapply(packages, library, character.only = TRUE)

classify_subtype <- function(exprData, geneSets) {
  par <- zscoreParam(exprData, geneSets)
  es <- gsva(par)
  max_score <- apply(X = es, MARGIN = 2, FUN = max) |>
    as.list()
  subtype <- character()
  for (i in 1:length(max_score)) {
    index <- which(es == max_score[i], arr.ind=TRUE)
    subtype <- append(subtype, rownames(es)[index[1]])
  }
  result <- list(max_score, subtype)
  return(result)
}

# Run GBM subtyping algorithm (GSVA) with gene sets from Verhaak et al. (2010)
verhaak_gene_lists <- c("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/Verhaak/MSigDB-CL.xml",
                        "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/Verhaak/MSigDB-MES.xml",
                        "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/Verhaak/MSigDB-NE.xml",
                        "/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/Verhaak/MSigDB-PN.xml")
verhaak_gs <- getBroadSets(verhaak_gene_lists)
verhaak_result <- classify_subtype(exprData = target_dsp_data@assayData[["Q3_Norm"]], geneSets = verhaak_gs)

# Run GBM subtyping algorithm (GSVA) with gene sets from Wang et al. (2017)
wang_MES <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/GBMsubtypeprofiles.3.xlsx", sheet = "Subtype Signatures", range = "A5:A54", col_names = FALSE) |>
  pull()
wang_PN <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/GBMsubtypeprofiles.3.xlsx", sheet = "Subtype Signatures", range = "I5:I54", col_names = FALSE) |>
  pull()
wang_CL <- read_excel("/Users/adamluo/Desktop/College/Lab [Shenderov]/Data Analysis/gbm-tma-dsp/import-files/GBMsubtypeprofiles.3.xlsx", sheet = "Subtype Signatures", range = "Q5:Q54", col_names = FALSE) |>
  pull()
wang_gs <- list(wang_MES, wang_PN, wang_CL)
names(wang_gs) <- c("MESENCHYMAL", "PRONEURAL", "CLASSICAL")
wang_result <- classify_subtype(exprData = target_dsp_data@assayData[["Q3_Norm"]], geneSets = wang_gs)

# Rename subtype data
rename_subtype <- function(result) {
  for (i in 1:length(result[[2]])) {
    if (grepl("CLASSICAL", result[[2]][i], fixed = TRUE)) {
      result[[2]][i] <- "CL"
    } else if (grepl("PRONEURAL", result[[2]][i], fixed = TRUE)) {
      result[[2]][i] <- "PN"
    } else if (grepl("NEURAL", result[[2]][i], fixed = TRUE)) {
      result[[2]][i] <- "NE"
    } else if (grepl("MESENCHYMAL", result[[2]][i], fixed = TRUE)) {
      result[[2]][i] <- "ME"
    }
  }
  return(result)
}
verhaak_result <- rename_subtype(verhaak_result)
wang_result <- rename_subtype(wang_result)

# Transfer GSVA results to S4 object
target_dsp_data[["Verhaak et al. Score"]] <- unlist(verhaak_result[1])
target_dsp_data[["Verhaak et al. Subtype"]] <- unlist(verhaak_result[2])
target_dsp_data[["Wang et al. Score"]] <- unlist(wang_result[1])
target_dsp_data[["Wang et al. Subtype"]] <- unlist(wang_result[2])

