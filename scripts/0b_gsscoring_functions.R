compDiapauseScore <- function(data.exprs.norm.z, diasigFile, method, datasetID){

  diaGenes <- toupper(diasigFile$gene)

  int <- intersect(rownames(data.exprs.norm.z), diaGenes)
  message("Number of diapause genes present in this dataset: ", length(int), " out of ", length(diaGenes), " genes.")

  # keepRows <- which(rownames(data.exprs.norm.z) %in% diaGenes)
  # keepRows <- diaGenes
  # keepRows <- intersect(rownames(data.exprs.norm.z), diaGenes)

  keepRows <- which(rownames(data.exprs.norm.z) %in% int)
  message("Are the indices of common genes extracted correctly: ", identical(rownames(data.exprs.norm.z)[which(rownames(data.exprs.norm.z) %in% int)], int))

  if (method=="weighting"){
    # make sure the right indices are extracted and multiplied
    # aa = data.exprs.norm.z[keepRows, ] * diasigFile$weight
    aa = data.exprs.norm.z[keepRows, ]
    bb = tibble::rownames_to_column(as.data.frame(aa), var = "gene")

    cc = diasigFile %>% dplyr::inner_join(bb, by = "gene")

    dd = cc$weight * cc[,3:ncol(cc)]
    rownames(dd) = cc$gene

    dd <- na.omit(dd)
    dia.score <- colMeans(dd)

    #Saving results in a new table
    dia.df <- data.frame(sample = colnames(data.exprs.norm.z), dia.score=dia.score)
  }
  if (method == "gsva"){  # Scoring using GSVA or ssGSEA
    res <- as.data.frame(GSVA::gsva(as.matrix(data.exprs.norm.z),
                                    list(diaGenes),
                                    method = "gsva"))
    dia.df <- data.frame(sample=colnames(res),
                         dia.score = as.numeric(as.vector(res[1,1:ncol(res)])))
  }
  if(method == "ssgsea"){
    res <- as.data.frame(GSVA::gsva(as.matrix(data.exprs.norm.z),
                                    list(diaGenes),
                                    method = "ssgsea"))
    dia.df <- data.frame(sample=colnames(res),
                         dia.score = as.numeric(as.vector(res[1,1:ncol(res)])))
  }

  # Accounting for differences in sample names in datasets
  if (datasetID=="GSE14333"){
    dia.df$sample <-sub(pattern = "\\..*", "", dia.df$sample) # GSE14333
  }
  else if(datasetID=="GSE38832"){
    dia.df$sample <- substr(dia.df$sample, start = 1, stop = 12)  # GSE38832
  }
  else{
    dia.df <- dia.df
  }

  return(dia.df)
}



compAutophagyScore <- function(data.exprs.norm.z, autoGenes, scoreType, datasetID){

  autoGenes <- toupper(autoGenes)

  int <- intersect(rownames(data.exprs.norm.z), autoGenes)
  message("Number of autophagy genes present in this dataset: ", length(int), " out of ", length(autoGenes), " genes.")

  # keepRows.ap <- which(rownames(data.exprs.norm.z) %in% autoGenes)
  # keepRows.ap <- autoGenes
  # keepRows.ap <- intersect(rownames(data.exprs.norm.z), autoGenes)
  keepRows.ap <- which(rownames(data.exprs.norm.z) %in% int)
  message("Are the indices of common genes extracted correctly: ", identical(rownames(data.exprs.norm.z)[which(rownames(data.exprs.norm.z) %in% int)], int))

  # make sure the right indices are extracted and multiplied
  bb = data.exprs.norm.z[keepRows.ap, ]
  bb <- na.omit(bb)
  ap.score <- vector(mode = "numeric", length = ncol(bb))
  if (scoreType == "wang"){
    for (kk in 1:ncol(bb)){
      ap.score[kk] = 0.1991*bb[which(rownames(bb) %in% "DLC1"),kk] + 0.2351*bb[which(rownames(bb) %in% "FKBP1B"),kk] + 0.0973*bb[which(rownames(bb) %in% "PEA15"),kk] + 0.5989*bb[which(rownames(bb) %in% "DNAJB1"),kk] + 0.1977*bb[which(rownames(bb) %in% "VAMP7"),kk] + 0.1584*bb[which(rownames(bb) %in% "PEX14"),kk]
    }
    #Saving results in a new table
    ap.df <- data.frame(sample = colnames(data.exprs.norm.z), ap.score=ap.score)
  }
  else{
    ap.score <- colMeans(bb)
    #Saving results in a new table
    ap.df <- data.frame(sample = colnames(data.exprs.norm.z), ap.score=ap.score)
  }
  if (scoreType == "gsva"){  # Scoring using GSVA or ssGSEA
    res <- as.data.frame(GSVA::gsva(as.matrix(data.exprs.norm.z),
                                    list(autoGenes),
                                    method = "gsva"))
    ap.df <- data.frame(sample=colnames(res),
                        ap.score = as.numeric(as.vector(res[1,1:ncol(res)])))
  }
  if (scoreType == "ssgsea"){  # Scoring using GSVA or ssGSEA
    res <- as.data.frame(GSVA::gsva(as.matrix(data.exprs.norm.z),
                                    list(autoGenes),
                                    method = "ssgsea"))
    ap.df <- data.frame(sample=colnames(res),
                        ap.score = as.numeric(as.vector(res[1,1:ncol(res)])))
  }

  # Accounting for differences in sample names in datasets
  if (datasetID=="GSE14333"){
    ap.df$sample <- sub("\\..*", "", ap.df$sample) # GSE14333
    }else if(datasetID=="GSE38832"){
      ap.df$sample <- substr(ap.df$sample, start = 1, stop = 12) # GSE38832
    }
  else{
    ap.df <- ap.df
  }


  return(ap.df)
}
