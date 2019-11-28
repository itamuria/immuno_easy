library(openxlsx)

ht <- read.xlsx("data/20191031_conteos_htseq.xlsx")
quan <- read.xlsx("data/20191031_conteos_quant3.xlsx")
sub <- read.xlsx("data/20191031_conteos_subread.xlsx")
cuf <- read.xlsx("data/20191105_conteos_cufflinks2.xlsx")

genes <- read.xlsx("data/LISTA_GENES_HEPA_6_SECUENCIAS_PEPTIDOS.xlsx")
genes_names2 <- genes$Gene_name
# genes_names2 <- gsub("\\x3bSTON1-GTF2A1L","",genes_names2)

ht[ht$hgnc_symbol%in%genes_names2,]

table(ht$hgnc_symbol%in%genes_names2)
table(quan$Gene_id%in%genes_names2)
table(sub$hgnc_symbol%in%genes_names2)
table(cuf$gene_short_name%in%genes_names2)

head(ht)
head(quan)
head(sub)
head(cuf)

ht2 <- ht[ht$hgnc_symbol%in%genes_names2,c("hgnc_symbol","Counts","Sample_1","Cuartiles_Sample_1")]
names(ht2) <- c("hgnc_symbol","Counts_ht","FPKM_ht","Cuartiles_ht")
quan2 <- quan[quan$Gene_id%in%genes_names2,]
names(quan2) <- c("Gene_id","Counts_quan","FPKM_quan","Cuartiles_quan")
sub2 <- sub[sub$hgnc_symbol%in%genes_names2,c("hgnc_symbol","Counts","Sample_1","Cuartiles_Sample_1")]
names(sub2) <- c("hgnc_symbol","Counts_sub","FPKM_sub","Cuartiles_sub")
cuf2 <- cuf[cuf$gene_short_name%in%genes_names2,c("gene_short_name","FPKM","Cuartiles")]
names(cuf2) <- c("gene_short_name","FPKM_cuf","Cuartiles_cuf")

m0 <- merge(genes[,1:2],ht2, by.x = "Gene_name", by.y = "hgnc_symbol", all.x = TRUE)

m1 <- merge(m0, quan2, by.x = "Gene_name", by.y = "Gene_id", all.x = TRUE)
m2 <- merge(m1, sub2, by.x = "Gene_name", by.y = "hgnc_symbol", all.x = TRUE)
m3 <- merge(m2, cuf2, by.x = "Gene_name", by.y = "gene_short_name", all.x = TRUE)

m3$Suma <- m3$Cuartiles_ht + m3$Cuartiles_quan + m3$Cuartiles_sub + m3$Cuartiles_cuf

write.xlsx(m3, file="20191105_merged_cuartiles.xlsx")
