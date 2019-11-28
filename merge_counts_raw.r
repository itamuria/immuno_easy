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


# FPKM + Cuartiles por separados ------------------------------------------



# htseq + subread
ht_sub <- merge(subread_raw[,c(1,6,7)], htseq_raw, by.x = "Geneid", by.y = "V1")
names(ht_sub) <- c("Geneid","SubRead_Length","SubRead_Count","htseq_Count")

# cuff + quant
table(cufflink_raw$gene_short_name!="-")
cufflink2 <- cufflink_raw[cufflink_raw$gene_short_name!="-",c(5,7,10)]
cufflink2_no <- cufflink_raw[cufflink_raw$gene_short_name=="-",]

cufflink2$gene_short_name <- as.character(cufflink2$gene_short_name)
cufflink2$borrar <- 0
add_1 <- NULL

for(v in 1:nrow(cufflink2))
{
  print(v)
  temp <- cufflink2$gene_short_name[v]
  st <- strsplit(temp,",")  
  len1 <- length(st[[1]])
  if(len1 > 1)
  {
    for(b in 1:len1)
    {
      cufflink2$borrar[v] <- 1
      add_1 <- c(add_1,
                 st[[1]][b], 
                cufflink2$locus[v], 
                cufflink2$FPKM[v])
     }
  }
}

gehitu <- data.frame(matrix(add_1, ncol = 3, byrow = T))
gehitu$borrar <- 0
names(gehitu) <- names(cufflink2)
cufflink3 <- rbind(cufflink2, gehitu)
cufflink4 <- cufflink3[cufflink3$borrar!=1,]

cuf_qua <- merge(cufflink4, quant_raw, by.x = "gene_short_name", by.y = "gen_names", all = TRUE)
names(cuf_qua)[c(3,5)] <- c("Cuff_FPKM","QUant_Counts")
head(cuf_qua)

# detectar gene names

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

# tarda
ens_gens <- ens2genenames (as.character(ht_sub$Geneid))
head(ens_gens)
head(ht_sub)
ht_sub$Geneid <- as.character(ht_sub$Geneid)
ht_sub$Geneid <- gsub('\\..+$', '', ht_sub$Geneid)

ht_sub2 <- merge(ht_sub, ens_gens, by.x = "Geneid", by.y = "ensembl_gene_id", all.x = TRUE)
head(ht_sub2)


head(cuf_qua)
kenduizena <- grep("HSCHR",cuf_qua$locus, value = TRUE)
kenduizena2 <- grep("^HG",cuf_qua$locus, value = TRUE)
kenduizena6 <- grep("^G",cuf_qua$locus, value = TRUE)
kenduizena3 <- grep("PATCH",cuf_qua$locus, value = TRUE)
kenduizena5 <- grep("^__",cuf_qua$gene_short_name, value = TRUE)
kenduizena4 <- unique(c(kenduizena, kenduizena2, kenduizena3,kenduizena6))

cuf_qua2 <- cuf_qua[-which(cuf_qua$locus%in%kenduizena4),]
cuf_qua2 <- cuf_qua2[-which(cuf_qua2$gene_short_name%in%kenduizena5),]


kenduizena7 <- grep("HSCHR",ht_sub2$chromosome_name, value = TRUE)
kenduizena8 <- grep("CHR_HG",ht_sub2$chromosome_name, value = TRUE)
kenduizena9 <- grep("^G",ht_sub2$chromosome_name, value = TRUE)
kenduizena11 <- grep("KI",ht_sub2$chromosome_name, value = TRUE)
kenduizena12 <- unique(c(kenduizena7, kenduizena8, kenduizena9,kenduizena11))

ht_sub3 <- ht_sub2[-which(ht_sub2$chromosome_name%in%kenduizena12),]


four_together <- merge(cuf_qua2, ht_sub3, by.x = "gene_short_name", by.y = "hgnc_symbol", all= TRUE)



# FPKM --------------------------------------------------------------------


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




















