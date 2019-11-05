
filename = "data/counts_subread.txt"
source = "subread"
symbol = FALSE
skiplines = 1
mfl_num = c(415)

counts2fpkm <- function(filename, source, symbol = FALSE, skiplines = 1, mfl_num)
{
  
  # Load data  ------------------------------------------------------------
  
  dat_subread <- read.table(filename,skip = skiplines, header = T)
  dat_subread <- dat_subread[,c("Geneid", names(dat_subread)[ncol(dat_subread)])]
  names(dat_subread) <- c("GeneID","Counts")
  
  # Update with biomart  ------------------------------------------------------------
  
  library("biomaRt")
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset("hsapiens_gene_ensembl", mart)

  ens <- as.character(dat_subread$GeneID)
  ens2 <- gsub('\\..+$', '', ens)
  
  annotLookup2 <- getBM(
    mart=mart,
    attributes=c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id", "chromosome_name","start_position","end_position"),
    filter="ensembl_gene_id",
    values=ens2,
    uniqueRows=TRUE)
  
  annotLookup2$length <- annotLookup2$end_position - annotLookup2$start_position
  
  # temp <- annotLookup2
  
  # Merging   ------------------------------------------------------------
  
  dat_subread$GeneID <- gsub('\\..+$', '', dat_subread$GeneID)
  
  # Remove duplicate agregating as sum of copies
  masde1 <- names(table(dat_subread$GeneID)[table(dat_subread$GeneID)>1])
  temp3 <- dat_subread[dat_subread$GeneID%in%masde1,]
  temp3b <- aggregate(temp3$Counts, by = list(temp3$GeneID), sum)
  unique_name <- unique(temp3$GeneID)
  dat_subread2 <- dat_subread[-which(dat_subread$GeneID%in%unique_name),]
  names(temp3b) <- names(dat_subread)
  dat_subread <- rbind(dat_subread2, temp3b)
  print(temp3)
  
  rownames(dat_subread) <- dat_subread[,1]

  gene_names2 <- merge(dat_subread, annotLookup2, by.x = "GeneID", by.y = "ensembl_gene_id")

  
  # Prepare for FPKM ------------------------------------------------------------
  
  library(countToFPKM)
  
  # repetir si solo hay una columna de conteos
  mfl <- c(rep(mfl_num,2))

  # annot3$ensembl_gene_id[annot3$ensembl_gene_id%in%rownames(total4b)]
  nam_delete <- names(table(gene_names2$ensembl_gene_id)[table(gene_names2$ensembl_gene_id)>1])
 
  print(nam_delete)
  
  gene_names2 <- unique(gene_names2)
  
  
  gene_names3 <- gene_names2[!gene_names2$ensembl_gene_id%in%nam_delete,]
  
  temp6 <- gene_names2[gene_names2$ensembl_gene_id%in%nam_delete,]
  temp6$GeneID <- paste0(temp6$GeneID,c("a","b"))
  gene_names4 <- rbind(gene_names3,temp6)
  
  df_calc <- data.frame(gene_names4$Counts,gene_names4$Counts)
  rownames(df_calc) <- gene_names4$GeneID
  colnames(df_calc) <- c("Sample_1","Sample_2")
  
  ## fpkm ------------------------------------------------------------
  
  fpkm_matrix <- countToFPKM::fpkm (df_calc, gene_names4$length, mfl)
  dim(fpkm_matrix)
  
  table(apply(fpkm_matrix, 1, sum) == 0)
  
  # save(fpkm_matrix, file="Cuartiles.RData")
  
  fpkm_matrix2 <- data.frame(fpkm_matrix)
  fpkm_matrix2$Gene_id <- rownames(fpkm_matrix)
  x <- fpkm_matrix2
  
  for(g in 1:2)
  {
    temp1 <- x[,g]
    qt <- quantile(temp1[temp1!=0], probs=0:4/4)
    q1 <- as.numeric(qt[1])
    q2 <- as.numeric(qt[2])
    q3 <- as.numeric(qt[3])
    q4 <- as.numeric(qt[4])
    
    qtal <- ifelse(temp1 <= q2, 1,
                   ifelse(temp1 <= q3 & temp1 > q2, 2,
                          ifelse(temp1 <= q4 & temp1 > q3, 3,
                                 ifelse(temp1 > q4, 4, NA))))
    table(qtal)
    x[,paste0("Cuartiles_",names(x)[g])] <- qtal
    qtal <- NA
  }
  
  x3 <- merge(x, gene_names4, by.x = "Gene_id", by.y = "GeneID")
  return(x3)
}

counts2fpkm (filename = "data/counts_subread.txt", source = "subread", symbol = FALSE, skiplines = 1, mfl_num = c(415))

# openxlsx::write.xlsx(x3,  file = "20191031_subread.xlsx")

filename = "data/tumor_rna_counts_union"
source = "htseq"
symbol = FALSE
skiplines = 0
mfl_num = c(415)

counts2fpkm <- function(filename, source, symbol = FALSE, skiplines = 1, mfl_num)
{
  
  # Load data  ------------------------------------------------------------
  
  dat_subread <- read.table(filename,skip = skiplines, header = F)
  # dat_subread <- dat_subread[,c("Geneid", names(dat_subread)[ncol(dat_subread)])]
  names(dat_subread) <- c("GeneID","Counts")
  
  # Update with biomart  ------------------------------------------------------------
  
  library("biomaRt")
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  
  ens <- as.character(dat_subread$GeneID)
  ens2 <- gsub('\\..+$', '', ens)
  
  annotLookup2 <- getBM(
    mart=mart,
    attributes=c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id", "chromosome_name","start_position","end_position"),
    filter="ensembl_gene_id",
    values=ens2,
    uniqueRows=TRUE)
  
  annotLookup2$length <- annotLookup2$end_position - annotLookup2$start_position
  
  # temp <- annotLookup2
  
  # Merging   ------------------------------------------------------------
  
  dat_subread$GeneID <- gsub('\\..+$', '', dat_subread$GeneID)
  
  # Remove duplicate agregating as sum of copies
  masde1 <- names(table(dat_subread$GeneID)[table(dat_subread$GeneID)>1])
  temp3 <- dat_subread[dat_subread$GeneID%in%masde1,]
  temp3b <- aggregate(temp3$Counts, by = list(temp3$GeneID), sum)
  unique_name <- unique(temp3$GeneID)
  dat_subread2 <- dat_subread[-which(dat_subread$GeneID%in%unique_name),]
  names(temp3b) <- names(dat_subread)
  dat_subread <- rbind(dat_subread2, temp3b)
  print(temp3)
  
  rownames(dat_subread) <- dat_subread[,1]
  
  gene_names2 <- merge(dat_subread, annotLookup2, by.x = "GeneID", by.y = "ensembl_gene_id")
  
  
  # Prepare for FPKM ------------------------------------------------------------
  
  library(countToFPKM)
  
  # repetir si solo hay una columna de conteos
  mfl <- c(rep(mfl_num,2))
  
  # annot3$ensembl_gene_id[annot3$ensembl_gene_id%in%rownames(total4b)]
  nam_delete <- names(table(gene_names2$ensembl_gene_id)[table(gene_names2$ensembl_gene_id)>1])
  
  print(nam_delete)
  
  gene_names2 <- unique(gene_names2)
  
  
  gene_names3 <- gene_names2[!gene_names2$ensembl_gene_id%in%nam_delete,]
  
  temp6 <- gene_names2[gene_names2$ensembl_gene_id%in%nam_delete,]
  temp6$GeneID <- paste0(temp6$GeneID,c("a","b"))
  gene_names4 <- rbind(gene_names3,temp6)
  
  df_calc <- data.frame(gene_names4$Counts,gene_names4$Counts)
  rownames(df_calc) <- gene_names4$GeneID
  colnames(df_calc) <- c("Sample_1","Sample_2")
  
  ## fpkm ------------------------------------------------------------
  
  fpkm_matrix <- countToFPKM::fpkm (df_calc, gene_names4$length, mfl)
  dim(fpkm_matrix)
  
  table(apply(fpkm_matrix, 1, sum) == 0)
  
  # save(fpkm_matrix, file="Cuartiles.RData")
  
  fpkm_matrix2 <- data.frame(fpkm_matrix)
  fpkm_matrix2$Gene_id <- rownames(fpkm_matrix)
  x <- fpkm_matrix2
  
  for(g in 1:2)
  {
    temp1 <- x[,g]
    qt <- quantile(temp1[temp1!=0], probs=0:4/4)
    q1 <- as.numeric(qt[1])
    q2 <- as.numeric(qt[2])
    q3 <- as.numeric(qt[3])
    q4 <- as.numeric(qt[4])
    
    qtal <- ifelse(temp1 <= q2, 1,
                   ifelse(temp1 <= q3 & temp1 > q2, 2,
                          ifelse(temp1 <= q4 & temp1 > q3, 3,
                                 ifelse(temp1 > q4, 4, NA))))
    table(qtal)
    x[,paste0("Cuartiles_",names(x)[g])] <- qtal
    qtal <- NA
  }
  
  x3 <- merge(x, gene_names4, by.x = "Gene_id", by.y = "GeneID")
  return(x3)
}

# openxlsx::write.xlsx(x3, file = "20191031_htseq.xlsx")

filename = "data/Hepa_pac6_quant3p.cnt"
source = "quant3p"
symbol = TRUE
skiplines = 0
mfl_num = c(415)

counts2fpkm <- function(filename, source, symbol = FALSE, skiplines = 1, mfl_num)
{
  
  # Load data  ------------------------------------------------------------
  
  dat_subread <- read.table(filename,skip = skiplines, header = T)
  # dat_subread <- dat_subread[,c("Geneid", names(dat_subread)[ncol(dat_subread)])]
  dat_subread <- data.frame(rownames(dat_subread),dat_subread[,1])
  names(dat_subread) <- c("GeneID","Counts")
  
  # Update with biomart  ------------------------------------------------------------
  
  if(!symbol)
  {
    library("biomaRt")
    mart <- useMart("ENSEMBL_MART_ENSEMBL")
    mart <- useDataset("hsapiens_gene_ensembl", mart)
    
    ens <- as.character(dat_subread$GeneID)
    ens2 <- gsub('\\..+$', '', ens)
    
    annotLookup2 <- getBM(
      mart=mart,
      attributes=c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id", "chromosome_name","start_position","end_position"),
      filter="ensembl_gene_id",
      values=ens2,
      uniqueRows=TRUE)
    
    annotLookup2$length <- annotLookup2$end_position - annotLookup2$start_position
    
    # temp <- annotLookup2
    
    # Merging   ------------------------------------------------------------
    
    dat_subread$GeneID <- gsub('\\..+$', '', dat_subread$GeneID)
    
    # Remove duplicate agregating as sum of copies
    masde1 <- names(table(dat_subread$GeneID)[table(dat_subread$GeneID)>1])
    temp3 <- dat_subread[dat_subread$GeneID%in%masde1,]
    temp3b <- aggregate(temp3$Counts, by = list(temp3$GeneID), sum)
    unique_name <- unique(temp3$GeneID)
    dat_subread2 <- dat_subread[-which(dat_subread$GeneID%in%unique_name),]
    names(temp3b) <- names(dat_subread)
    dat_subread <- rbind(dat_subread2, temp3b)
    print(temp3)
    
    rownames(dat_subread) <- dat_subread[,1]
    
    gene_names2 <- merge(dat_subread, annotLookup2, by.x = "GeneID", by.y = "ensembl_gene_id")
  }
  
  
  # Prepare for FPKM ------------------------------------------------------------
  
  library(countToFPKM)
  
  # repetir si solo hay una columna de conteos
  mfl <- c(rep(mfl_num,2))
  
  # obtener longitudes con biomart
  x4 <- x3[grep("_",x3$chromosome_name,invert = TRUE),]
  x4 <- x4[x4$hgnc_symbol != "",]
  
  x5 <- x4[which(x4$hgnc_symbol%in%names(table(x4$hgnc_symbol))[table(x4$hgnc_symbol)>1]),]
  x6 <- x4[-which(x4$hgnc_symbol%in%names(table(x4$hgnc_symbol))[table(x4$hgnc_symbol)>1]),]
  
  x5b <- x5[order(x5$hgnc_symbol),]
  ken5 <- c(1,3,6,7,9,11,13,15,17,19,21,23,25)
  x7 <- x5b[-ken5,]
  
  x8 <- rbind(x6,x7)


  mm <- merge(dat_subread, x8[,c("hgnc_symbol","length")],by.x = "GeneID",by.y = "hgnc_symbol")
  
  if(!symbol)
  {
    # annot3$ensembl_gene_id[annot3$ensembl_gene_id%in%rownames(total4b)]
    nam_delete <- names(table(gene_names2$ensembl_gene_id)[table(gene_names2$ensembl_gene_id)>1])
    
    print(nam_delete)
    
    gene_names2 <- unique(gene_names2)
    
    
    gene_names3 <- gene_names2[!gene_names2$ensembl_gene_id%in%nam_delete,]
    
    temp6 <- gene_names2[gene_names2$ensembl_gene_id%in%nam_delete,]
    temp6$GeneID <- paste0(temp6$GeneID,c("a","b"))
    gene_names4 <- rbind(gene_names3,temp6)
    
    df_calc <- data.frame(gene_names4$Counts,gene_names4$Counts)
    rownames(df_calc) <- gene_names4$GeneID
    colnames(df_calc) <- c("Sample_1","Sample_2")
  } else {
    df_calc <- data.frame(mm$Counts,mm$Counts)
    rownames(df_calc) <- mm$GeneID
    colnames(df_calc) <- c("Sample_1","Sample_2")
  }
  
 
  
  ## fpkm ------------------------------------------------------------
  
  fpkm_matrix <- countToFPKM::fpkm (df_calc, mm$length, mfl)
  dim(fpkm_matrix)
  
  table(apply(fpkm_matrix, 1, sum) == 0)
  
  # save(fpkm_matrix, file="Cuartiles.RData")
  
  fpkm_matrix2 <- data.frame(fpkm_matrix)
  fpkm_matrix2$Gene_id <- rownames(fpkm_matrix)
  x <- fpkm_matrix2
  
  for(g in 1:2)
  {
    temp1 <- x[,g]
    qt <- quantile(temp1[temp1!=0], probs=0:4/4)
    q1 <- as.numeric(qt[1])
    q2 <- as.numeric(qt[2])
    q3 <- as.numeric(qt[3])
    q4 <- as.numeric(qt[4])
    
    qtal <- ifelse(temp1 <= q2, 1,
                   ifelse(temp1 <= q3 & temp1 > q2, 2,
                          ifelse(temp1 <= q4 & temp1 > q3, 3,
                                 ifelse(temp1 > q4, 4, NA))))
    table(qtal)
    x[,paste0("Cuartiles_",names(x)[g])] <- qtal
    qtal <- NA
  }
  
  df_calc$GeneID <- rownames(df_calc)
  x3 <- merge(x, df_calc, by.x = "Gene_id", by.y = "GeneID")
  return(x3)
}

# openxlsx::write.xlsx(x3, file = "20191031_quant3.xlsx")









# 
# 
# 
# load("G:/Mi unidad/Hepatil/neoantigen_detection/20191009_prueba_template.RData")
# 
# # df_all_fin
# # STON1\x3bSTON1-GTF2A1L
# # TRIM73\x3bTRIM74
# # USP17L11\x3bUSP17L18\x3bUSP17L20
# 
# df_all_fin$Gene_name[grep("STON1",df_all_fin$Gene_name)] <- "STON1-GTF2A1L"
# df_all_fin$Gene_name[grep("TRIM74",df_all_fin$Gene_name)] <- "TRIM73"
# df_all_fin$Gene_name[grep("USP17L11",df_all_fin$Gene_name)] <- "USP17L11"
# 
# df_all_fin$Gene_name[grep("FAM208B",df_all_fin$Gene_name)] <- "TASOR2"
# 
# which(df_all_fin$Gene_name%in%c("FAM208B","TRIM74"))
# 
# grep("GKN1",x3$hgnc_symbol,value = T)
# grep("GKN1",df_all_fin$Gene_name)
# 
# for(t in 1:nrow(df_all_fin))
# {
#   print(paste0("t: ",t))
#   
#   x_temp <- x3[x3$hgnc_symbol == df_all_fin$Gene_name[t],]
#   str_temp <- strsplit(df_all_fin$variant_key[t]," ")
#   
#   print(str_temp)
#   
#   if(nrow(x_temp)>0)
#   {
#     if(nrow(x_temp) == 1)
#     {
#       temp5 <- strsplit(str_temp[[1]], ":")
#       chr <- temp5[[1]][1]
#       pos <- as.numeric(temp5[[1]][2])
#       
#       fpkm_val <- x_temp$Sample_1[x_temp$chromosome_name == chr & x_temp$start_position < pos & x_temp$end_position > pos]
#       cuart_val <- x_temp$Cuartiles_Sample_1[x_temp$chromosome_name == chr & x_temp$start_position < pos & x_temp$end_position > pos]
#       
#       if(length(frpk_val)>0 & !identical(fpkm_val, numeric(0)))
#       {
#         df_all_fin$Cuartil[t] <- cuart_val
#         df_all_fin$FPKM[t] <- fpkm_val
#       }
#       
#     } else {
#       for(g in 1:nrow(x_temp))
#       {
#         temp5 <- strsplit(str_temp[[1]], ":")
#         chr <- temp5[[1]][1]
#         pos <- as.numeric(temp5[[1]][2])
#         
#         fpkm_val <- x_temp$Sample_1[x_temp$chromosome_name == chr & x_temp$start_position < pos & x_temp$end_position > pos]
#         cuart_val <- x_temp$Cuartiles_Sample_1[x_temp$chromosome_name == chr & x_temp$start_position < pos & x_temp$end_position > pos]
#         
#         if(length(frpk_val)>0 & !identical(fpkm_val, numeric(0)))
#         {
#           df_all_fin$Cuartil[t] <- cuart_val
#           df_all_fin$FPKM[t] <- fpkm_val
#         }
#       }
#     }
#   }
# }
