
# vector_names <- as.character(ht_sub$Geneid)

ens2genenames <- function(vector_names)
{
  library("biomaRt")
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  
  ens <- as.character(vector_names)
  ens2 <- gsub('\\..+$', '', ens)
  
  annotLookup2 <- getBM(
    mart=mart,
    attributes=c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id", "chromosome_name","start_position","end_position"),
    filter="ensembl_gene_id",
    values=ens2,
    uniqueRows=TRUE)
  
  annotLookup2$length <- annotLookup2$end_position - annotLookup2$start_position
  return(annotLookup2)
}


# LOAD DATA -------------------------------------------------

setwd("D:/HEPATIL/Pac22/")

dir()

subread_raw_file <- "counts_subread.txt"
cufflink_raw_file <- "genes.fpkm_tracking"
quant_raw_file <- "Hepa_pac6_quant3p.cnt"
htseq_raw_file <- "tumor_rna_counts_union"

subread_raw <- read.table(subread_raw_file, header = T)
head(subread_raw)
cufflink_raw <- read.table(cufflink_raw_file, header = T, sep = "\t")
head(cufflink_raw)

f <- cufflink_raw[cufflink_raw$Gene.Name!="-",]
table(cufflink_raw$FPKM>0)

quant_raw <- read.table(quant_raw_file, header = T)
quant_raw$gen_names <- rownames(quant_raw)
head(quant_raw)
htseq_raw <- read.table(htseq_raw_file)
head(htseq_raw)

filename = subread_raw
source = "subread"
symbol = FALSE
skiplines = 1
mfl_num = c(415)

counts2fpkm_subread <- function(filename, symbol = FALSE, skiplines = 1, mfl_num)
{
  
  # Load data  ------------------------------------------------------------
  
  dat_subread <- filename
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

semi_subread <- counts2fpkm_subread (filename = subread_raw, symbol = FALSE, skiplines = 1, mfl_num = c(415))



# from FPMK to cuartiles --------------------------------------------------


fpkm2cuartiles_cuff <- function(cufflink_raw)
{
  cf <- cufflink_raw
  g <- 10
  temp1 <- cf[,g]
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
  cf[,paste0("Cuartiles")] <- qtal
  return(cf)
  # openxlsx::write.xlsx(cf, file = "20191105_conteos_cufflinks2.xlsx")
} 
 
semi_cuff <- fpkm2cuartiles_cuff(cufflink_raw)
semi_cuff2 <- fpkm2cuartiles_cuff(cufflink_raw[cufflink_raw$gene_short_name!="-",])


# openxlsx::write.xlsx(x3,  file = "20191031_subread.xlsx")

filename = htseq_raw
source = "htseq"
symbol = FALSE
skiplines = 0
mfl_num = c(415)

counts2fpkm_htseq <- function(filename, symbol = FALSE, mfl_num = c(415))
{
  
  # Load data  ------------------------------------------------------------
  
  dat_subread <- filename
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

semi_htseq <- counts2fpkm_htseq (htseq_raw)

# Quant desde symbolo -----------------------------------------------------


filename = quant_raw
source = "quant3p"
symbol = TRUE
skiplines = 0
mfl_num = c(415)

counts2fpkm_quant <- function(filename, symbol = TRUE, mfl_num)
{
  
  # Load data  ------------------------------------------------------------
  
  dat_subread <- quant_raw
  # dat_subread <- dat_subread[,c("Geneid", names(dat_subread)[ncol(dat_subread)])]
  dat_subread <- data.frame(rownames(dat_subread),dat_subread[,1])
  names(dat_subread) <- c("GeneID","Counts")
  
  # Update with biomart  ------------------------------------------------------------
  

    library("biomaRt")
    mart <- useMart("ENSEMBL_MART_ENSEMBL")
    mart <- useDataset("hsapiens_gene_ensembl", mart)
    
    ens2 <- as.character(dat_subread$GeneID)
    # ens2 <- gsub('\\..+$', '', ens)
    
    annotLookup2 <- getBM(
      mart=mart,
      attributes=c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id", "chromosome_name","start_position","end_position"),
      filter="hgnc_symbol",
      values=ens2,
      uniqueRows=TRUE)
    
    annotLookup2$length <- annotLookup2$end_position - annotLookup2$start_position
    
    # temp <- annotLookup2
    
    # Merging   ------------------------------------------------------------
    
    # dat_subread$GeneID <- gsub('\\..+$', '', dat_subread$GeneID)
    
    # Remove duplicate agregating as sum of copies
    # masde1 <- names(table(dat_subread$GeneID)[table(dat_subread$GeneID)>1])
    # temp3 <- dat_subread[dat_subread$GeneID%in%masde1,]
    # temp3b <- aggregate(temp3$Counts, by = list(temp3$GeneID), sum)
    # unique_name <- unique(temp3$GeneID)
    # dat_subread2 <- dat_subread[-which(dat_subread$GeneID%in%unique_name),]
    # names(temp3b) <- names(dat_subread)
    # dat_subread <- rbind(dat_subread2, temp3b)
    # print(temp3)
    
    rownames(dat_subread) <- dat_subread[,1]
    
    gene_names2 <- merge(dat_subread, annotLookup2, by.x = "GeneID", by.y = "hgnc_symbol")
  
    gene_names2 <- gene_names2[nchar(gene_names2$chromosome_name)<3,]
    # gene_names4 <- gene_names2[nchar(gene_names2$chromosome_name)>3,]
  
  # Prepare for FPKM ------------------------------------------------------------
  
  library(countToFPKM)
  
  # repetir si solo hay una columna de conteos
  mfl <- c(rep(mfl_num,2))
  
  # obtener longitudes con biomart
  # x4 <- x3[grep("_",x3$chromosome_name,invert = TRUE),]
  # x4 <- x4[x4$hgnc_symbol != "",]
  # 
  # x5 <- x4[which(x4$hgnc_symbol%in%names(table(x4$hgnc_symbol))[table(x4$hgnc_symbol)>1]),]
  # x6 <- x4[-which(x4$hgnc_symbol%in%names(table(x4$hgnc_symbol))[table(x4$hgnc_symbol)>1]),]
  # 
  # x5b <- x5[order(x5$hgnc_symbol),]
  # ken5 <- c(1,3,6,7,9,11,13,15,17,19,21,23,25)
  # x7 <- x5b[-ken5,]
  # 
  # x8 <- rbind(x6,x7)
  # 
  # 
  # mm <- merge(dat_subread, x8[,c("hgnc_symbol","length")],by.x = "GeneID",by.y = "hgnc_symbol")
  
  # if(!symbol)
  # {
  #   # annot3$ensembl_gene_id[annot3$ensembl_gene_id%in%rownames(total4b)]
    nam_delete <- names(table(gene_names2$GeneID)[table(gene_names2$GeneID)>1])

    print(nam_delete)

    gene_names2 <- unique(gene_names2)


    gene_names3 <- gene_names2[!gene_names2$GeneID%in%nam_delete,]

    temp6 <- gene_names2[gene_names2$GeneID%in%nam_delete,]
    temp6$GeneID <- paste0(temp6$GeneID,c("a","b"))
    gene_names4 <- rbind(gene_names3,temp6)

    df_calc <- data.frame(gene_names4$Counts,gene_names4$Counts)
    rownames(df_calc) <- gene_names4$GeneID
    colnames(df_calc) <- c("Sample_1","Sample_2")
  # } else {
  #   df_calc <- data.frame(mm$Counts,mm$Counts)
  #   rownames(df_calc) <- mm$GeneID
  #   colnames(df_calc) <- c("Sample_1","Sample_2")
  # }
  
  # df_calc <- data.frame(gene_names2$Counts, gene_names2$Counts)
  # colnames(df_calc) <- c("Sample_1","Sample_2")
  # rownames(df_calc) <- gene_names2$GeneID
  
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
  
  df_calc$GeneID <- rownames(df_calc)
  x3 <- merge(x, df_calc, by.x = "Gene_id", by.y = "GeneID")
  x4 <- x3[,c(1,2,4,6)]
  names(x4) <- c("Gene_id","FPKM","Cuartiles","Counts")
  return(x4)
}

# openxlsx::write.xlsx(x3, file = "20191031_quant3.xlsx")

semi_quant <- counts2fpkm_quant(quant_raw,mfl_num = c(415))


# Juntar todos los conteos ------------------------------------------------

from_multinames_to_rows <- function(dataset, colu_name)
{
  colu_position <- which(names(dataset) == colu_name)
  
  which_delete <- grep(",",dataset[,colu_name])
  dataset$delete <- 0  
  dataset$delete[which_delete] <- 1
  
  add_1 <- NULL
  
  for(v in which_delete)
  {
    print(v)
    temp <- as.character(dataset[v,colu_name])
    st <- strsplit(temp,",")  
    len1 <- length(st[[1]])
    if(len1 > 1)
    {
      for(b in 1:len1)
      {
        # cufflink2$borrar[v] <- 1
        add_1 <- c(add_1, as.matrix(dataset[v,-which(names(dataset)==colu_name)]), st[[1]][b])
      }
    }
  }
  
  gehitu <- data.frame(matrix(add_1, ncol = ncol(dataset), byrow = T))
  # gehitu$delete <- 0
  dataset2 <- dataset[dataset$delete == 0,]
  
  izenak <- names(dataset)[-colu_position]
  names(gehitu) <- c(izenak, colu_name)
  gehitu$delete <- 0
  
  dataset3 <- rbind(dataset2, gehitu)
  return(dataset3)
}

# dataset <- semi_cuffb
# colu_name <- "Symbol"



head(semi_cuff)
semi_cuffb <- semi_cuff[,c(1,5,7,10,14)]
names(semi_cuffb) <- c("Gene_id","Symbol","cuff_locus","cuff_FPKM","cuff_Cuartiles")
head(semi_cuffb)
semi_cuffb10 <- from_multinames_to_rows (dataset = semi_cuffb, colu_name = "Symbol")
semi_cuffb10 <- semi_cuffb10[,-which(names(semi_cuffb10)=="delete")]


head(semi_cuff2)
semi_cuff2b <- semi_cuff2[,c(1,5,7,10,14)]
names(semi_cuff2b) <- c("Gene_id","Symbol","cuff2_locus","cuff2_FPKM","cuff2_Cuartiles")
head(semi_cuff2b)
semi_cuff2b10 <- from_multinames_to_rows (dataset = semi_cuff2b, colu_name = "Symbol")
semi_cuff2b10 <- semi_cuff2b10[,-which(names(semi_cuff2b10)=="delete")]


head(semi_htseq)
semi_htseq2 <- semi_htseq[,c(1,7,2,4,6)]
semi_htseq2$ht_locus <- paste0(semi_htseq$chromosome_name,":",semi_htseq$start_position,"-",semi_htseq$end_position)
names(semi_htseq2) <- c("Gene_id","Symbol","ht_FPKM","ht_Cuartiles","ht_Counts","ht_locus")
head(semi_htseq2)
semi_htseq2 <- semi_htseq2[which(semi_htseq2$ht_locus%in%grep("^C|^G",semi_htseq2$ht_locus, value=TRUE, invert = TRUE)),]

head(semi_quant)
names(semi_quant) <- c("Symbol","Quan_FPKM","Quan_Cuartiles","Quan_Counts")
head(semi_quant)

head(semi_subread)
semi_subread2 <- semi_subread[,c(1,7,2,4,6)]
semi_subread2$ht_locus <- paste0(semi_subread$chromosome_name,":",semi_subread$start_position,"-",semi_subread$end_position)
names(semi_subread2) <- c("Gene_id","Symbol","SR_FPKM","SR_Cuartiles","SR_Counts","SR_locus")
head(semi_subread2)
semi_subread2 <- semi_subread2[which(semi_subread2$SR_locus%in%grep("^C|^G",semi_subread2$SR_locus, value=TRUE, invert = TRUE)),]

# merging

m1_total <- merge(semi_subread2, semi_quant, by = "Symbol", all = TRUE)
# quitando los vacioes en symbol
m1_logico <- merge(semi_subread2[semi_subread2$Symbol!="",], semi_quant, by = "Symbol")

semi_cuff_c <- data.frame(semi_cuff2b10,semi_cuffb10$cuff_Cuartiles[semi_cuffb10$Symbol!="-"])
head(semi_cuff_c)
names(semi_cuff_c)[6] <- "cuff_Cuartiles"
semi_cuff_c$cuff2_locus <- as.character(semi_cuff_c$cuff2_locus)
semi_cuff_c <- semi_cuff_c[which(semi_cuff_c$cuff2_locus%in%grep("^H|^G|^C",semi_cuff_c$cuff2_locus, value=TRUE, invert = TRUE)),]

# m2_total <- merge(semi_cuff_c, semi_htseq2, by = "Symbol", all = TRUE)
# quitando los vacioes en symbol
# table(semi_cuffb$Symbol!="-")
m2_logico <- merge(semi_cuff_c, semi_htseq2, by = "Symbol", all.y = TRUE) 

# m3_total <- merge(m1_total,m2_total, by = "Symbol", all = TRUE)
m3_logico <- merge(m1_logico,m2_logico, by = "Symbol", all = TRUE)

# openxlsx::write.xlsx(m2_total, file="subread_quant_merged.xlsx")
# openxlsx::write.xlsx(m1_total, file="cuff_htseq_merged.xlsx")

# openxlsx::write.xlsx(m3_logico, file="four_together.xlsx")

# filtrar con nuestros datos

nuestros <- openxlsx::read.xlsx("20190525_listadeGenes-2.xlsx")
ngenes <- unique(nuestros$All)

mfinal <- m3_logico[m3_logico$Symbol%in%ngenes,]
ngenes[!ngenes%in%m3_logico$Symbol]

head(mfinal)
mfinal2 <- mfinal[,c(1,2,15,6,11,19,5,9,18,3,7,12,16,4,8,13,14,17)]
head(mfinal2)


