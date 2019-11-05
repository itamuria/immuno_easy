mut38_2_excel <- function(vcffile)
{
  # vcffile <- mut38
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  fox <- grep("FOXOG",vcffile$V9)
  pgt <- grep("PGT",vcffile$V9)
  
  table(vcffile$V9)
  
  vcffile$fox <- ifelse(rownames(vcffile)%in%fox,"Fox",NA)
  vcffile$pgt <- ifelse(rownames(vcffile)%in%pgt,"Pgt",NA)
  
  vcffile$group <- 0
  vcffile$group[vcffile$fox == "Fox" &  is.na(vcffile$pgt)] <- 1
  vcffile$group[vcffile$fox == "Fox" &  vcffile$pgt == "Pgt"] <- 2
  vcffile$group[is.na(vcffile$fox) &  is.na(vcffile$pgt)] <- 3
  vcffile$group[vcffile$pgt == "Pgt" &  is.na(vcffile$fox)] <- 4
  vcffile$ID <- 1:(nrow(vcffile))
  
  table(vcffile$group)
  name_group <- names(table(vcffile$group))
  
  # Extract names
  u1 <- unique(vcffile$V9[vcffile$fox == "Fox" &  vcffile$pgt == "Pgt"])
  # if(!is.na(u1))
  # {
  u1 <- u1[!is.na(u1)] 
  names9 <- strsplit(u1,":")[[1]]
  len9 <- length(names9)
  
  # GROUP 1 
  g1_10 <- as.character(vcffile$V10[vcffile$group==1])
  g1_11 <- as.character(vcffile$V11[vcffile$group==1])
  id_1 <- vcffile$ID[vcffile$group==1]
  
  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = 9,byrow = TRUE))
  col10c <- cbind(col10b[,1:6],rep(NA,nrow(col10b)),rep(NA,nrow(col10b)),col10b[,7:9])
  names(col10c) <- names9
  
  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = 9,byrow = TRUE))
  col11c <- cbind(col11b[,1:6],rep(NA,nrow(col11b)),rep(NA,nrow(col11b)),col11b[,7:9])
  names(col11c) <- names9
  
  # GROUP 2
  g2_10 <- as.character(vcffile$V10[vcffile$group==2])
  g2_11 <- as.character(vcffile$V11[vcffile$group==2])
  id_2 <- vcffile$ID[vcffile$group==2]
  
  col210 <- unlist(strsplit(g2_10,":"))
  col210b <- data.frame(matrix(col210,ncol = 11,byrow = TRUE))
  col210c <- col210b
  names(col210c) <- names9
  
  col211 <- unlist(strsplit(g2_11,":"))
  col211b <- data.frame(matrix(col211,ncol = 11,byrow = TRUE))
  col211c <- col211b
  names(col211c) <- names9
  # }
  
  # 
  # if(sum(name_group==c(1,2))==0)
  # {
  #   # GROUP 3
  #   g3_10 <- as.character(vcffile$V10[vcffile$group==3])
  #   g3_11 <- as.character(vcffile$V11[vcffile$group==3])
  #   id_3 <- vcffile$ID[vcffile$group==3]
  #   
  #   col310 <- unlist(strsplit(g3_10,":"))
  #   col310b <- data.frame(matrix(col310,ncol = 6,byrow = TRUE))
  #   col310c <- cbind(col310b[,1:5],rep(NA,nrow(col310b)),rep(NA,nrow(col310b)),rep(NA,nrow(col310b)),col310b[,6])
  #   names(col310c) <- names9
  #   
  #   col311 <- unlist(strsplit(g3_11,":"))
  #   col311b <- data.frame(matrix(col311,ncol = 6,byrow = TRUE))
  #   col311c <- cbind(col311b[,1:5],rep(NA,nrow(col311b)),rep(NA,nrow(col311b)),rep(NA,nrow(col311b)),col311b[,6:8])
  #   names(col311c) <- names9
  #   
  #   # GROUP 4
  #   g4_10 <- as.character(vcffile$V10[vcffile$group==4])
  #   g4_11 <- as.character(vcffile$V11[vcffile$group==4])
  #   id_4 <- vcffile$ID[vcffile$group==4]
  #   
  #   col410 <- unlist(strsplit(g4_10,":"))
  #   col410b <- data.frame(matrix(col410,ncol = 10,byrow = TRUE))
  #   col410c <- cbind(col410b[,1:5],rep(NA,nrow(col410b)),col410b[,6:10])
  #   names(col410c) <- names9
  #   
  #   col411 <- unlist(strsplit(g4_11,":"))
  #   col411b <- data.frame(matrix(col411,ncol = 10,byrow = TRUE))
  #   col411c <- cbind(col411b[,1:5],rep(NA,nrow(col411b)),col411b[,6:10])
  #   names(col411c) <- names9
  #   
  #   id_t <- c(id_1,id_2,id_3,id_4)
  #   dfcol <- rbind(col10c,col210c,col310c,col410c)
  #   dfcol3 <- rbind(col11c,col211c,col311c,col411c)
  #   dfcol2 <- cbind(id_t,dfcol,dfcol3)
  #   dfcol2 <- dfcol2[order(dfcol2$id_t),]
  #   names(dfcol2)[13:23]<- paste0(names(dfcol2)[2:12],"_tum")
  #   
  #   ad1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AD),",")),ncol = 2,byrow = TRUE))
  #   names(ad1) <- c("AD_ref","AD_alt")
  #   gss1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$QSS),",")),ncol = 2,byrow = TRUE))
  #   names(gss1) <- c("QSS_ref","QSS_alt")
  #   
  #   ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AD_tum),",")),ncol = 2,byrow = TRUE))
  #   names(ad1b) <- c("AD_ref_tum","AD_alt_tum")
  #   gss1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$QSS_tum),",")),ncol = 2,byrow = TRUE))
  #   names(gss1b) <- c("QSS_ref_tum","QSS_alt_tum")
  # }
  
  # GROUP 3
  g3_10 <- as.character(vcffile$V10[vcffile$group==3])
  g3_11 <- as.character(vcffile$V11[vcffile$group==3])
  id_3 <- vcffile$ID[vcffile$group==3]
  
  col310 <- unlist(strsplit(g3_10,":"))
  col310b <- data.frame(matrix(col310,ncol = 8,byrow = TRUE))
  col310c <- cbind(col310b[,1:5],rep(NA,nrow(col310b)),rep(NA,nrow(col310b)),rep(NA,nrow(col310b)),col310b[,6:8])
  names(col310c) <- names9
  
  col311 <- unlist(strsplit(g3_11,":"))
  col311b <- data.frame(matrix(col311,ncol = 8,byrow = TRUE))
  col311c <- cbind(col311b[,1:5],rep(NA,nrow(col311b)),rep(NA,nrow(col311b)),rep(NA,nrow(col311b)),col311b[,6:8])
  names(col311c) <- names9
  
  # GROUP 4
  g4_10 <- as.character(vcffile$V10[vcffile$group==4])
  g4_11 <- as.character(vcffile$V11[vcffile$group==4])
  id_4 <- vcffile$ID[vcffile$group==4]
  
  col410 <- unlist(strsplit(g4_10,":"))
  col410b <- data.frame(matrix(col410,ncol = 10,byrow = TRUE))
  col410c <- cbind(col410b[,1:5],rep(NA,nrow(col410b)),col410b[,6:10])
  names(col410c) <- names9
  
  col411 <- unlist(strsplit(g4_11,":"))
  col411b <- data.frame(matrix(col411,ncol = 10,byrow = TRUE))
  col411c <- cbind(col411b[,1:5],rep(NA,nrow(col411b)),col411b[,6:10])
  names(col411c) <- names9
  
  id_t <- c(id_1,id_2,id_3,id_4)
  dfcol <- rbind(col10c,col210c,col310c,col410c)
  dfcol3 <- rbind(col11c,col211c,col311c,col411c)
  dfcol2 <- cbind(id_t,dfcol,dfcol3)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  names(dfcol2)[13:23]<- paste0(names(dfcol2)[2:12],"_tum")
  
  ad1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AD),",")),ncol = 2,byrow = TRUE))
  names(ad1) <- c("AD_ref","AD_alt")
  gss1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$QSS),",")),ncol = 2,byrow = TRUE))
  names(gss1) <- c("QSS_ref","QSS_alt")
  
  ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AD_tum),",")),ncol = 2,byrow = TRUE))
  names(ad1b) <- c("AD_ref_tum","AD_alt_tum")
  gss1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$QSS_tum),",")),ncol = 2,byrow = TRUE))
  names(gss1b) <- c("QSS_ref_tum","QSS_alt_tum")
  
  # FINAL FORMAT 
  formatfin <- cbind(dfcol2[,c(1,2)],ad1,dfcol2[,c(4:8)],gss1,dfcol2[,c(11,13)],ad1b,dfcol2[,c(15:20)],gss1b,dfcol2[,c(22,23)])
  
  ### INFO
  f8 <- vcffile$V8
  
  library(tidyr)
  s1 <- separate(data = vcffile, col = V8, into = paste0("a",1:48), sep = ";")
  s2 <- s1[,c(8:55)]
  s2[grep("RPA",s2$a6),names(s2)] <- s2[grep("RPA",s2$a6),names(s2)[-c(6:8)]]
  
  snames <- unlist(strsplit(as.character(s2[1,]),"="))
  snames2 <- snames[seq(1,95,by = 2)]
  
  for(h in c(1:41,43:47))
  {
    print(h)
    s2[,h] <- gsub(paste0(snames2[h],"="),"",s2[,h])
    names(s2)[h] <- snames2[h]
  }
  names(s2)[42]<- "GERP"
  
  dftot <- cbind(vcffile[,1:7],s2,formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")
  
  for(q in c(2,58:63,65:67,69:74,77:80))
  {
    dftot[,q] <- as.numeric(as.character(dftot[,q]))
  }
  
  dftot$Cob_DNA_nor <- dftot$AD_ref+dftot$AD_alt
  dftot$Cob_DNA_tum <- dftot$AD_ref_tum+dftot$AD_alt_tum
  
  dftot <- dftot[,!is.na(names(dftot))]
  
  return(dftot)
  
}

mut41_2_excel <- function(vcffile)
{
  # vcffile <- mut41
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  pgt <- grep("PGT",vcffile$V9)
  
  vcffile$pgt <- ifelse(rownames(vcffile)%in%pgt,"Pgt","Other")
  
  vcffile$group <- 0
  vcffile$group[vcffile$pgt == "Pgt"] <- 1
  vcffile$ID <- 1:(nrow(vcffile))
  
  table(vcffile$group)
  
  # Extract names
  u1 <- unique(vcffile$V9[vcffile$pgt == "Pgt"])
  u1 <- u1[!is.na(u1)] 
  names9 <- strsplit(u1,":")[[1]]
  len9 <- length(names9)
  
  # GROUP 1 
  g1_10 <- as.character(vcffile$V10[vcffile$group==1])
  g1_11 <- as.character(vcffile$V11[vcffile$group==1])
  id_1 <- vcffile$ID[vcffile$group==1]
  
  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = 9,byrow = TRUE))
  col10c <- col10b
  names(col10c) <- names9
  
  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = 9,byrow = TRUE))
  col11c <- col11b
  names(col11c) <- paste0(names9,"_2")
  
  # GROUP 0
  g2_10 <- as.character(vcffile$V10[vcffile$group==0])
  g2_11 <- as.character(vcffile$V11[vcffile$group==0])
  id_2 <- vcffile$ID[vcffile$group==0]
  
  col210 <- unlist(strsplit(g2_10,":"))
  col210b <- data.frame(matrix(col210,ncol = 6,byrow = TRUE))
  col210c <- cbind(col210b,rep(NA,nrow(col210b)),rep(NA,nrow(col210b)),rep(NA,nrow(col210b)))
  names(col210c) <- names9
  
  col211 <- unlist(strsplit(g2_11,":"))
  col211b <- data.frame(matrix(col211,ncol = 6,byrow = TRUE))
  col211c <- cbind(col211b,rep(NA,nrow(col211b)),rep(NA,nrow(col211b)),rep(NA,nrow(col211b)))
  names(col211c) <- paste0(names9,"_2")
  
  id_t <- c(id_1,id_2)
  dfcol <- rbind(col10c,col210c)
  dfcolb <- rbind(col11c,col211c)
  dfcol2 <- cbind(id_t,dfcol,dfcolb)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  
  ###
  
  l <- strsplit(as.character(dfcol2$AD),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  ad1 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(ad1) <- c("AD_ref","AD_alt","AD_3","AD_4")
  
  l <- strsplit(as.character(dfcol2$F1R2),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  F1R2 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(F1R2) <- c("F1R2_ref","F1R2_alt","F1R2_3","F1R2_4")
  
  l <- strsplit(as.character(dfcol2$F2R1),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  F2R1 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(F2R1) <- c("F1R2_ref","F1R2_alt","F1R2_3","F1R2_4")
  
  l <- strsplit(as.character(dfcol2$AF),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  AF <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(AF) <- c("AF_ref","AF_alt","AF_3","AF_4")
  
  ###
  
  l <- strsplit(as.character(dfcol2$AD_2),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  ad1_2 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(ad1_2) <- c("AD_ref2","AD_alt2","AD_32","AD_42")
  
  l <- strsplit(as.character(dfcol2$F1R2_2),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  F1R2_2 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(F1R2_2) <- c("F1R2_ref2","F1R2_alt2","F1R2_32","F1R2_42")
  
  l <- strsplit(as.character(dfcol2$F2R1_2),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  F2R1_2 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(F2R1_2) <- c("F1R2_ref2","F1R2_alt2","F1R2_32","F1R2_42")
  
  l <- strsplit(as.character(dfcol2$AF_2),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  AF_2 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(AF_2) <- c("AF_ref2","AF_alt2","AF_32","AF_42")
  
  # FINAL FORMAT 
  formatfin <- cbind(dfcol2[,c(1,2)],ad1,AF,dfcol2[,5],F1R2,F2R1,dfcol2[,c(8,11)],ad1_2,AF_2,dfcol2[,14],F1R2_2,F2R1_2,dfcol2[,c(17:19)])
  names(formatfin)[c(11,30)] <- c("DP","DP_2")
  
  ### INFO
  f8 <- vcffile$V8
  
  library(tidyr)
  s1 <- separate(data = vcffile, col = V8, into = paste0("a",1:48), sep = ";")
  s2 <- s1[,c(8:55)]
  snames <- unlist(strsplit(as.character(s2[1,]),"="))
  snames2 <- snames[seq(1,95,by = 2)]
  
  # ### break
  # vapply(strsplit(as.character(f8[1]),";"), function(x) length(x))
  # 
  # le <- NA
  # for(h in 1:length(f8))
  # {
  #   ss <- strsplit(as.character(f8[h]),";")
  #   le <- c(le,length(ss[[1]]))
  # }
  # ###
  
  for(h in c(1:41,43:47))
  {
    print(h)
    s2[,h] <- gsub(paste0(snames2[h],"="),"",s2[,h])
    names(s2)[h] <- snames2[h]
  }
  names(s2)[42]<- "GERP"
  
  dftot <- cbind(vcffile[,1:7],s2,formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")
  
  return(dftot)
  
}

mut41_2_excel_v2 <- function(vcffile, txtfile)
{
  # vcffile <- mut41
  # txtfile <- mut41txt
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  pgt <- grep("PGT",vcffile$V9)
  
  vcffile$pgt <- ifelse(rownames(vcffile)%in%pgt,"Pgt","Other")
  
  vcffile$group <- 0
  vcffile$group[vcffile$pgt == "Pgt"] <- 1
  vcffile$ID <- 1:(nrow(vcffile))
  
  table(vcffile$group)
  
  # Extract names
  u1 <- unique(vcffile$V9[vcffile$pgt == "Pgt"])
  u1 <- u1[!is.na(u1)]
  names9 <- strsplit(u1,":")[[1]]
  len9 <- length(names9)
  
  # GROUP 1
  g1_10 <- as.character(vcffile$V10[vcffile$group==1])
  g1_11 <- as.character(vcffile$V11[vcffile$group==1])
  id_1 <- vcffile$ID[vcffile$group==1]
  
  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = 9,byrow = TRUE))
  col10c <- col10b
  names(col10c) <- names9
  
  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = 9,byrow = TRUE))
  col11c <- col11b
  names(col11c) <- paste0(names9,"_2")
  
  # GROUP 0
  g2_10 <- as.character(vcffile$V10[vcffile$group==0])
  g2_11 <- as.character(vcffile$V11[vcffile$group==0])
  id_2 <- vcffile$ID[vcffile$group==0]
  
  col210 <- unlist(strsplit(g2_10,":"))
  col210b <- data.frame(matrix(col210,ncol = 6,byrow = TRUE))
  col210c <- cbind(col210b,rep(NA,nrow(col210b)),rep(NA,nrow(col210b)),rep(NA,nrow(col210b)))
  names(col210c) <- names9
  
  col211 <- unlist(strsplit(g2_11,":"))
  col211b <- data.frame(matrix(col211,ncol = 6,byrow = TRUE))
  col211c <- cbind(col211b,rep(NA,nrow(col211b)),rep(NA,nrow(col211b)),rep(NA,nrow(col211b)))
  names(col211c) <- paste0(names9,"_2")
  
  id_t <- c(id_1,id_2)
  dfcol <- rbind(col10c,col210c)
  dfcolb <- rbind(col11c,col211c)
  dfcol2 <- cbind(id_t,dfcol,dfcolb)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  
  ###
  
  l <- strsplit(as.character(dfcol2$AD),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  ad1 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(ad1) <- c("AD_ref","AD_alt","AD_3","AD_4")
  
  l <- strsplit(as.character(dfcol2$F1R2),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  F1R2 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(F1R2) <- c("F1R2_ref","F1R2_alt","F1R2_3","F1R2_4")
  
  l <- strsplit(as.character(dfcol2$F2R1),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  F2R1 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(F2R1) <- c("F2R1_ref","F2R1_alt","F2R1_3","F2R1_4")
  
  l <- strsplit(as.character(dfcol2$AF),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  AF <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(AF) <- c("AF_ref","AF_alt","AF_3","AF_4")
  
  ###
  
  l <- strsplit(as.character(dfcol2$AD_2),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  ad1_2 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(ad1_2) <- c("AD_ref2","AD_alt2","AD_32","AD_42")
  
  l <- strsplit(as.character(dfcol2$F1R2_2),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  F1R2_2 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(F1R2_2) <- c("F1R2_ref2","F1R2_alt2","F1R2_32","F1R2_42")
  
  l <- strsplit(as.character(dfcol2$F2R1_2),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  F2R1_2 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(F2R1_2) <- c("F2R1_ref2","F2R1_alt2","F2R1_32","F2R1_42")
  
  l <- strsplit(as.character(dfcol2$AF_2),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  AF_2 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(AF_2) <- c("AF_ref2","AF_alt2","AF_32","AF_42")
  
  # FINAL FORMAT
  formatfin <- cbind(dfcol2[,c(1,2)],ad1,AF,dfcol2[,5],F1R2,F2R1,dfcol2[,c(8,11)],ad1_2,AF_2,dfcol2[,14],F1R2_2,F2R1_2,dfcol2[,c(17:19)])
  names(formatfin)[c(11,30)] <- c("DP","DP_2")
  
  
  dftot <- cbind(vcffile[,1:7],formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")
  
  mut41_txt2 <- txtfile[,1:10]
  mut41_txt2$merge <- paste0(mut41_txt2$Chr,mut41_txt2$Start)
  dftot$merge <- paste0(dftot$Chr,dftot$Pos)
  
  dftot2 <- merge(mut41_txt2,dftot, by = c("merge"))
  dftotx <- merge(mut41_txt2,dftot, by = c("merge"), all.x = TRUE)
  dftoty <- merge(mut41_txt2,dftot, by = c("merge"), all.y = TRUE)
  
  return(list(dftot2,dftotx,dftoty))
  
}

# cuidado con GT, no quitar o cambiar codigo despues de cambiar
strelka_2_excel <- function(vcffile)
{
  # vcffile <- strelka
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  vcffile$ID <- 1:(nrow(vcffile))
  
  # Extract names
  u1 <- unique(vcffile$V9)
  u1 <- u1[!is.na(u1)] 
  names9 <- strsplit(u1,":")[[1]]
  # names9 <- names9[-1]
  # names9 <- names9[-1]
  len9 <- length(names9)
  
  # GROUP 1 
  g1_10 <- as.character(vcffile$V10)
  nchar(g1_10)
  
  g1_11 <- as.character(vcffile$V11)
  id_1 <- vcffile$ID
  
  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = len9,byrow = TRUE))
  col10c <- col10b
  names(col10c) <- names9
  
  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = len9,byrow = TRUE))
  col11c <- col11b
  names(col11c) <- names9
  
  
  id_t <- id_1
  dfcol <- rbind(col10c)
  dfcol3 <- rbind(col11c)
  dfcol2 <- cbind(id_t,dfcol,dfcol3)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  nnc <- ncol(dfcol2)
  nnc2 <- (nnc-1) / 2 
  nnc3 <- nnc2 * 2 + 1
  names(dfcol2)[(nnc2+2):nnc3]<- paste0(names(dfcol2)[2:(nnc2+1)],"_2")
  
  ad1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AU),",")),ncol = 2,byrow = TRUE))
  names(ad1) <- c("AU_T1","AU_T2")
  gss1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$CU),",")),ncol = 2,byrow = TRUE))
  names(gss1) <- c("CU_T1","CU_T2")
  gss2 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$GU),",")),ncol = 2,byrow = TRUE))
  names(gss2) <- c("GU_T1","GU_T2")
  gss3 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TU),",")),ncol = 2,byrow = TRUE))
  names(gss3) <- c("TU_T1","TU_T2")
  
  ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AU_2),",")),ncol = 2,byrow = TRUE))
  names(ad1b) <- c("AU_T1_2","AU_T2_2")
  gss1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$CU_2),",")),ncol = 2,byrow = TRUE))
  names(gss1b) <- c("CU_T1_2","CU_T2_2")
  gss2b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$GU_2),",")),ncol = 2,byrow = TRUE))
  names(gss2b) <- c("GU_T1_2","GU_T2_2")
  gss3b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TU),",")),ncol = 2,byrow = TRUE))
  names(gss3b) <- c("TU_T1_2","TU_T2_2")
  
  # FINAL FORMAT 
  formatfin <- cbind(dfcol2[,c(1:5)],ad1,gss1,gss2,gss3,dfcol2[,c(10:13)],ad1b,gss1b,gss2b,gss3b)
  
  ### INFO
  f8 <- vcffile$V8
  
  library(tidyr)
  s1 <- separate(data = vcffile, col = V8, into = paste0("a",1:38), sep = ";")
  s2 <- s1[,c(8:27)]
  snames <- unlist(strsplit(as.character(s2[1,]),"="))
  snames2 <- snames[c(1,seq(2,38,by = 2))]
  
  for(h in c(2:20))
  {
    print(h)
    s2[,h] <- gsub(paste0(snames2[h],"="),"",s2[,h])
    names(s2)[h] <- snames2[h]
  }
  
  dftot <- cbind(vcffile[,1:7],s2,formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")
  
  # plot(as.numeric(dftot$AU_T1),as.numeric(dftot$AU_T2))
  # cor(as.numeric(dftot$AU_T1),as.numeric(dftot$AU_T2))
  # cor.test(as.numeric(dftot$AU_T1),as.numeric(dftot$AU_T2))
  # 
  # plot(as.numeric(dftot$GU_T1),as.numeric(dftot$GU_T2))
  # cor(as.numeric(dftot$GU_T1),as.numeric(dftot$GU_T2))
  # cor.test(as.numeric(dftot$GU_T1),as.numeric(dftot$GU_T2))
  
  for(y in 29:ncol(dftot))
  {
    dftot[,y] <- as.numeric(as.character(dftot[,y]))
  }
  
  
  dftot$Ref_T1 <- 0
  dftot$Alt_T1 <- 0
  dftot$Ref_T2 <- 0
  dftot$Alt_T2 <- 0
  
  for(q in 1:nrow(dftot))
  {
    print(paste0(q,"/",nrow(dftot)))
    dftot$Ref_T1[q] <- ifelse(dftot$Ref[q]=="A",dftot[q,"AU_T1"],
                              ifelse(dftot$Ref[q]=="C",dftot[q,"CU_T1"],
                                     ifelse(dftot$Ref[q]=="G",dftot[q,"GU_T1"],
                                            ifelse(dftot$Ref[q]=="T",dftot[q,"TU_T1"],NA))))
    
    dftot$Ref_T2[q]  <- ifelse(dftot$Ref[q]=="A",dftot[q,"AU_T2"],
                               ifelse(dftot$Ref[q]=="C",dftot[q,"CU_T2"],
                                      ifelse(dftot$Ref[q]=="G",dftot[q,"GU_T2"],
                                             ifelse(dftot$Ref[q]=="T",dftot[q,"TU_T2"],NA))))
    
    dftot$Alt_T1[q]  <- ifelse(dftot$Alt[q]=="A",dftot[q,"AU_T1"],
                               ifelse(dftot$Alt[q]=="C",dftot[q,"CU_T1"],
                                      ifelse(dftot$Alt[q]=="G",dftot[q,"GU_T1"],
                                             ifelse(dftot$Alt[q]=="T",dftot[q,"TU_T1"],NA))))
    
    dftot$Alt_T2[q]  <- ifelse(dftot$Alt[q]=="A",dftot[q,"AU_T2"],
                               ifelse(dftot$Alt[q]=="C",dftot[q,"CU_T2"],
                                      ifelse(dftot$Alt[q]=="G",dftot[q,"GU_T2"],
                                             ifelse(dftot$Alt[q]=="T",dftot[q,"TU_T2"],NA))))
  }
  
  
  
  
  dftot$Ref_T1_2 <- 0
  dftot$Alt_T1_2 <- 0
  dftot$Ref_T2_2 <- 0
  dftot$Alt_T2_2 <- 0
  
  for(q in 1:nrow(dftot))
  {
    print(paste0(q,"/",nrow(dftot)))
    dftot$Ref_T1_2[q] <- ifelse(dftot$Ref[q]=="A",dftot[q,"AU_T1_2"],
                                ifelse(dftot$Ref[q]=="C",dftot[q,"CU_T1_2"],
                                       ifelse(dftot$Ref[q]=="G",dftot[q,"GU_T1_2"],
                                              ifelse(dftot$Ref[q]=="T",dftot[q,"TU_T1_2"],NA))))
    
    dftot$Ref_T2_2[q]  <- ifelse(dftot$Ref[q]=="A",dftot[q,"AU_T2_2"],
                                 ifelse(dftot$Ref[q]=="C",dftot[q,"CU_T2_2"],
                                        ifelse(dftot$Ref[q]=="G",dftot[q,"GU_T2_2"],
                                               ifelse(dftot$Ref[q]=="T",dftot[q,"TU_T2_2"],NA))))
    
    dftot$Alt_T1_2[q]  <- ifelse(dftot$Alt[q]=="A",dftot[q,"AU_T1_2"],
                                 ifelse(dftot$Alt[q]=="C",dftot[q,"CU_T1_2"],
                                        ifelse(dftot$Alt[q]=="G",dftot[q,"GU_T1_2"],
                                               ifelse(dftot$Alt[q]=="T",dftot[q,"TU_T1_2"],NA))))
    
    dftot$Alt_T2_2[q]  <- ifelse(dftot$Alt[q]=="A",dftot[q,"AU_T2_2"],
                                 ifelse(dftot$Alt[q]=="C",dftot[q,"CU_T2_2"],
                                        ifelse(dftot$Alt[q]=="G",dftot[q,"GU_T2_2"],
                                               ifelse(dftot$Alt[q]=="T",dftot[q,"TU_T2_2"],NA))))
  }
  
  
  
  dftote <- dftot
  dftot[dftot$NT=="hom",c("Ref_T1","Ref_T2","Ref_T1_2","Ref_T2_2")] <- dftote[dftote$NT=="hom",c("Alt_T1","Alt_T2","Alt_T1_2","Alt_T2_2")]
  dftot[dftot$NT=="hom",c("Alt_T1","Alt_T2","Alt_T1_2","Alt_T2_2")] <- dftote[dftote$NT=="hom",c("Ref_T1","Ref_T2","Ref_T1_2","Ref_T2_2")]
  
  
  dftot$VAF_T1 <- dftot$Alt_T1_2/(dftot$Alt_T1_2+dftot$Ref_T1_2)
  dftot$VAF_T2 <- dftot$Alt_T2_2/(dftot$Alt_T2_2+dftot$Ref_T2_2)
  
  # ff <- dftot$Ref_T1-dftot$Alt_T1
  
  # aggregate(ff,by=list(dftot$NT),FUN=function(x) c(min(x),mean(x),max(x)))
  return(dftot)
  
}

strelka_2_excel_snvs_fast <- function(vcffile)
{
  # vcffile <- strelka
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  vcffile$ID <- 1:(nrow(vcffile))
  
  # Extract names
  u1 <- unique(vcffile$V9)
  u1 <- u1[!is.na(u1)] 
  names9 <- strsplit(u1,":")[[1]]
  # names9 <- names9[-1]
  # names9 <- names9[-1]
  len9 <- length(names9)
  
  # GROUP 1 
  g1_10 <- as.character(vcffile$V10)
  nchar(g1_10)
  
  g1_11 <- as.character(vcffile$V11)
  id_1 <- vcffile$ID
  
  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = len9,byrow = TRUE))
  col10c <- col10b
  names(col10c) <- names9
  
  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = len9,byrow = TRUE))
  col11c <- col11b
  names(col11c) <- names9
  
  
  id_t <- id_1
  dfcol <- rbind(col10c)
  dfcol3 <- rbind(col11c)
  dfcol2 <- cbind(id_t,dfcol,dfcol3)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  nnc <- ncol(dfcol2)
  nnc2 <- (nnc-1) / 2 
  nnc3 <- nnc2 * 2 + 1
  names(dfcol2)[(nnc2+2):nnc3]<- paste0(names(dfcol2)[2:(nnc2+1)],"_2")
  
  ad1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AU),",")),ncol = 2,byrow = TRUE))
  names(ad1) <- c("AU_T1_normal","AU_T2_normal")
  gss1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$CU),",")),ncol = 2,byrow = TRUE))
  names(gss1) <- c("CU_T1_normal","CU_T2_normal")
  gss2 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$GU),",")),ncol = 2,byrow = TRUE))
  names(gss2) <- c("GU_T1_normal","GU_T2_normal")
  gss3 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TU),",")),ncol = 2,byrow = TRUE))
  names(gss3) <- c("TU_T1_normal","TU_T2_normal")
  
  ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AU_2),",")),ncol = 2,byrow = TRUE))
  names(ad1b) <- c("AU_T1_tumor","AU_T2_tumor")
  gss1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$CU_2),",")),ncol = 2,byrow = TRUE))
  names(gss1b) <- c("CU_T1_tumor","CU_T2_tumor")
  gss2b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$GU_2),",")),ncol = 2,byrow = TRUE))
  names(gss2b) <- c("GU_T1_tumor","GU_T2_tumor")
  gss3b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TU_2),",")),ncol = 2,byrow = TRUE))
  names(gss3b) <- c("TU_T1_tumor","TU_T2_tumor")
  
  # FINAL FORMAT 
  formatfin <- cbind(dfcol2[,c(1:5)],ad1,gss1,gss2,gss3,dfcol2[,c(10:13)],ad1b,gss1b,gss2b,gss3b)
  
  ### INFO
  f8 <- vcffile$V8
  
  library(tidyr)
  s1 <- separate(data = vcffile, col = V8, into = paste0("a",1:38), sep = ";")
  s2 <- s1[,c(8:27)]
  snames <- unlist(strsplit(as.character(s2[1,]),"="))
  snames2 <- snames[c(1,seq(2,38,by = 2))]
  
  for(h in c(2:20))
  {
    print(h)
    s2[,h] <- gsub(paste0(snames2[h],"="),"",s2[,h])
    names(s2)[h] <- snames2[h]
  }
  
  dftot <- cbind(vcffile[,1:7],s2,formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")
  
  # plot(as.numeric(dftot$AU_T1),as.numeric(dftot$AU_T2))
  # cor(as.numeric(dftot$AU_T1),as.numeric(dftot$AU_T2))
  # cor.test(as.numeric(dftot$AU_T1),as.numeric(dftot$AU_T2))
  # 
  # plot(as.numeric(dftot$GU_T1),as.numeric(dftot$GU_T2))
  # cor(as.numeric(dftot$GU_T1),as.numeric(dftot$GU_T2))
  # cor.test(as.numeric(dftot$GU_T1),as.numeric(dftot$GU_T2))
  
  for(y in 29:ncol(dftot))
  {
    dftot[,y] <- as.numeric(as.character(dftot[,y]))
  }
  
  
  dftot$Ref_T1_normal <- 0
  dftot$Alt_T1_normal <- 0
  dftot$Ref_T2_normal <- 0
  dftot$Alt_T2_normal <- 0
  
  dftot$Ref_T1_normal <- ifelse(dftot$Ref=="A",dftot[,"AU_T1_normal"],
                                ifelse(dftot$Ref=="C",dftot[,"CU_T1_normal"],
                                       ifelse(dftot$Ref=="G",dftot[,"GU_T1_normal"],
                                              ifelse(dftot$Ref=="T",dftot[,"TU_T1_normal"],NA))))
  
  dftot$Ref_T2_normal  <- ifelse(dftot$Ref=="A",dftot[,"AU_T2_normal"],
                                 ifelse(dftot$Ref=="C",dftot[,"CU_T2_normal"],
                                        ifelse(dftot$Ref=="G",dftot[,"GU_T2_normal"],
                                               ifelse(dftot$Ref=="T",dftot[,"TU_T2_normal"],NA))))
  
  dftot$Alt_T1_normal  <- ifelse(dftot$Alt=="A",dftot[,"AU_T1_normal"],
                                 ifelse(dftot$Alt=="C",dftot[,"CU_T1_normal"],
                                        ifelse(dftot$Alt=="G",dftot[,"GU_T1_normal"],
                                               ifelse(dftot$Alt=="T",dftot[,"TU_T1_normal"],NA))))
  
  dftot$Alt_T2_normal  <- ifelse(dftot$Alt=="A",dftot[,"AU_T2_normal"],
                                 ifelse(dftot$Alt=="C",dftot[,"CU_T2_normal"],
                                        ifelse(dftot$Alt=="G",dftot[,"GU_T2_normal"],
                                               ifelse(dftot$Alt=="T",dftot[,"TU_T2_normal"],NA))))
  
  
  dftot$Ref_T1_tumor <- 0
  dftot$Alt_T1_tumor <- 0
  dftot$Ref_T2_tumor <- 0
  dftot$Alt_T2_tumor <- 0
  
  dftot$Ref_T1_tumor <- ifelse(dftot$Ref=="A",dftot[,"AU_T1_tumor"],
                               ifelse(dftot$Ref=="C",dftot[,"CU_T1_tumor"],
                                      ifelse(dftot$Ref=="G",dftot[,"GU_T1_tumor"],
                                             ifelse(dftot$Ref=="T",dftot[,"TU_T1_tumor"],NA))))
  
  dftot$Ref_T2_tumor  <- ifelse(dftot$Ref=="A",dftot[,"AU_T2_tumor"],
                                ifelse(dftot$Ref=="C",dftot[,"CU_T2_tumor"],
                                       ifelse(dftot$Ref=="G",dftot[,"GU_T2_tumor"],
                                              ifelse(dftot$Ref=="T",dftot[,"TU_T2_tumor"],NA))))
  
  dftot$Alt_T1_tumor  <- ifelse(dftot$Alt=="A",dftot[,"AU_T1_tumor"],
                                ifelse(dftot$Alt=="C",dftot[,"CU_T1_tumor"],
                                       ifelse(dftot$Alt=="G",dftot[,"GU_T1_tumor"],
                                              ifelse(dftot$Alt=="T",dftot[,"TU_T1_tumor"],NA))))
  
  dftot$Alt_T2_tumor  <- ifelse(dftot$Alt=="A",dftot[,"AU_T2_tumor"],
                                ifelse(dftot$Alt=="C",dftot[,"CU_T2_tumor"],
                                       ifelse(dftot$Alt=="G",dftot[,"GU_T2_tumor"],
                                              ifelse(dftot$Alt=="T",dftot[,"TU_T2_tumor"],NA))))
  
  
  
  dftote <- dftot
  dftot[dftot$NT=="hom",c("Ref_T1_normal","Ref_T2_normal","Ref_T1_tumor","Ref_T2_tumor")] <- dftote[dftote$NT=="hom",c("Alt_T1_normal","Alt_T2_normal","Alt_T1_tumor","Alt_T2_tumor")]
  dftot[dftot$NT=="hom",c("Alt_T1_normal","Alt_T2_normal","Alt_T1_tumor","Alt_T2_tumor")] <- dftote[dftote$NT=="hom",c("Ref_T1_normal","Ref_T2_normal","Ref_T1_tumor","Ref_T2_tumor")]
  
  
  dftot$VAF_T1 <- dftot$Alt_T1_tumor/(dftot$Alt_T1_tumor+dftot$Ref_T1_tumor)
  dftot$VAF_T2 <- dftot$Alt_T2_tumor/(dftot$Alt_T2_tumor+dftot$Ref_T2_tumor)
  
  # ff <- dftot$Ref_T1-dftot$Alt_T1
  
  # aggregate(ff,by=list(dftot$NT),FUN=function(x) c(min(x),mean(x),max(x)))
  
  return(dftot)
  
}

strelka_indels_2_excel <- function(vcffile)
{
  # vcffile <- strelka_indels
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  vcffile$ID <- 1:(nrow(vcffile))
  
  # Extract names
  u1 <- unique(vcffile$V9)
  u1 <- u1[!is.na(u1)] 
  names9 <- strsplit(u1,":")[[1]]
  # names9 <- names9[-1]
  len9 <- length(names9)
  
  # GROUP 1 
  g1_10 <- as.character(vcffile$V10)
  g1_11 <- as.character(vcffile$V11)
  id_1 <- vcffile$ID
  
  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = len9,byrow = TRUE))
  col10c <- col10b
  names(col10c) <- names9
  
  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = len9,byrow = TRUE))
  col11c <- col11b
  names(col11c) <- names9
  
  
  id_t <- id_1
  dfcol <- rbind(col10c)
  dfcol3 <- rbind(col11c)
  dfcol2 <- cbind(id_t,dfcol,dfcol3)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  nnc <- ncol(dfcol2)
  nnc2 <- (nnc-1) / 2 
  nnc3 <- nnc2 * 2 + 1
  names(dfcol2)[(nnc2+2):nnc3]<- paste0(names(dfcol2)[2:(nnc2+1)],"_2")
  
  ad1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TAR),",")),ncol = 2,byrow = TRUE))
  names(ad1) <- c("TAR_T1","TAR_T2")
  gss1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TIR),",")),ncol = 2,byrow = TRUE))
  names(gss1) <- c("TIR_T1","TIR_T2")
  gss2 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TOR),",")),ncol = 2,byrow = TRUE))
  names(gss2) <- c("TOR_T1","TOR_T2")
  
  ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TAR_2),",")),ncol = 2,byrow = TRUE))
  names(ad1b) <- c("TAR_T1_2","TAR_T2_2")
  gss1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TIR_2),",")),ncol = 2,byrow = TRUE))
  names(gss1b) <- c("TIR_T1_2","TIR_T2_2")
  gss2b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$TOR_2),",")),ncol = 2,byrow = TRUE))
  names(gss2b) <- c("TOR_T1_2","TOR_T2_2")
  
  
  # FINAL FORMAT 
  formatfin <- cbind(dfcol2[,c(1:3)],ad1,gss1,gss2,dfcol2[,c(7:12)],ad1b,gss1b,gss2b)
  
  ### INFO
  f8 <- vcffile$V8
  
  library(tidyr)
  s1 <- separate(data = vcffile, col = V8, into = paste0("a",1:21), sep = ";")
  s2 <- s1[,c(8:30)]
  snames <- unlist(strsplit(as.character(s2[1,]),"="))
  snames2 <- snames[c(1,seq(2,40,by = 2))]
  
  for(h in c(2:23))
  {
    print(h)
    s2[,h] <- gsub(paste0(snames2[h],"="),"",s2[,h])
    names(s2)[h] <- snames2[h]
  }
  
  dftot <- cbind(vcffile[,1:7],s2,formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")
  
  dftot$TIR_T1_2 <- as.numeric(as.character(dftot$TIR_T1_2))
  dftot$TAR_T1_2 <- as.numeric(as.character(dftot$TAR_T1_2))
  
  dftot$vaf <- dftot$TIR_T1_2/(dftot$TIR_T1_2+dftot$TAR_T1_2)
  
  naw <- which(is.na(names(dftot)))
  if(length(naw)>0)
  {
    dftot <- dftot[,-naw]
  }
  
  pos1 <- which(names(dftot)=="DP")
  
  for(g in pos1:ncol(dftot))
  {
    dftot[,g] <- as.numeric(as.character(dftot[,g]))
  }
  
  return(dftot)
  
}

somaticsniper_2_excel <- function(vcffile)
{
  # vcffile <- somaticsniper
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  vcffile$ID <- 1:(nrow(vcffile))
  
  # Extract names
  u1 <- unique(vcffile$V9)
  u1 <- u1[!is.na(u1)] 
  names9 <- strsplit(u1,":")[[1]]
  len9 <- length(names9)
  
  # GROUP 1 
  g1_10 <- as.character(vcffile$V10)
  g1_11 <- as.character(vcffile$V11)
  id_1 <- vcffile$ID
  
  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = len9,byrow = TRUE))
  col10c <- col10b
  names(col10c) <- names9
  
  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = len9,byrow = TRUE))
  col11c <- col11b
  names(col11c) <- names9
  
  
  id_t <- id_1
  dfcol <- rbind(col10c)
  dfcol3 <- rbind(col11c)
  dfcol2 <- cbind(id_t,dfcol,dfcol3)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  nnc <- ncol(dfcol2)
  nnc2 <- (nnc-1) / 2 
  nnc3 <- nnc2 * 2 + 1
  names(dfcol2)[(nnc2+2):nnc3]<- paste0(names(dfcol2)[2:(nnc2+1)],"_2")
  
  ad1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$DP4),",")),ncol = 4,byrow = TRUE))
  names(ad1) <- c("Norm_DP_ref_forw","Norm_DP_ref_rev","Norm_DP_alt_forw","Norm_DP_alt_rev")
  gss1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$BCOUNT),",")),ncol = 4,byrow = TRUE))
  names(gss1) <- c("NormCOUNT_A","NormCOUNT_C","NormCOUNT_G","NormCOUNT_T")
  
  ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$DP4_2),",")),ncol = 4,byrow = TRUE))
  names(ad1b) <- c("Tum_DP_ref_forw","Tum_DP_ref_rev","Tum_DP_alt_forw","Tum_DP_alt_rev")
  gss1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$BCOUNT_2),",")),ncol = 4,byrow = TRUE))
  names(gss1b) <- c("TumCOUNT_A","TumCOUNT_C","TumCOUNT_G","TumCOUNT_T")
  
  for(t in 1:4)
  {
    ad1[,t] <- as.numeric(as.character(ad1[,t]))
    ad1b[,t] <- as.numeric(as.character(ad1b[,t]))
  }
  
  alt1 <- ad1$Norm_DP_alt_forw+ad1$Norm_DP_alt_rev
  alt2 <- ad1b$Tum_DP_alt_forw+ad1b$Tum_DP_alt_rev
  
  dfcol2$DP <- as.numeric(as.character(dfcol2$DP))
  dfcol2$DP_2 <- as.numeric(as.character(dfcol2$DP_2))
  
  
  vaf1 <- alt1/dfcol2$DP
  vaf2 <- alt2/dfcol2$DP_2
  
  # FINAL FORMAT 
  formatfin <- cbind(dfcol2[,c(1:4)],ad1,alt1,vaf1,gss1,dfcol2[,c(7:17)],ad1b,alt2,vaf2,gss1b,dfcol2[,c(20:27)])
  
  ### INFO
  f8 <- vcffile$V8
  
  library(tidyr)
  s1 <- separate(data = vcffile, col = V8, into = paste0("a",1:42), sep = ";")
  s2 <- s1[,c(8:49)]
  snames <- unlist(strsplit(as.character(s2[1,]),"="))
  snames2 <- snames[c(seq(2,83,by = 2))]
  
  for(h in c(1:35,37:41))
  {
    print(h)
    s2[,h+1] <- gsub(paste0(snames2[h],"="),"",s2[,h+1])
    names(s2)[h+1] <- snames2[h]
  }
  
  names(s2)[37] <- "GERP"
  # s2[,37] <- gsub(paste0("GERP++_RS=.","",s2[,h+1]))
  dftot <- cbind(vcffile[,1:7],s2,formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")
  
  
  nakendu <- sapply(dftot, function(e) {
    len1 <- length(e)
    pr90 <- len1 * 90 / 100
    pr90 < sum(is.na(e))
  })
  dftot <- dftot[,!nakendu]
  
  return(dftot)
  
}

varscan_2_excel <- function(vcffile)
{
  # vcffile <- varscan
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  vcffile$ID <- 1:(nrow(vcffile))
  
  # Extract names
  u1 <- unique(vcffile$V9)
  u1 <- u1[!is.na(u1)] 
  names9 <- strsplit(u1,":")[[1]]
  len9 <- length(names9)
  
  # GROUP 1 
  g1_10 <- as.character(vcffile$V10)
  g1_11 <- as.character(vcffile$V11)
  id_1 <- vcffile$ID
  
  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = len9,byrow = TRUE))
  col10c <- col10b
  names(col10c) <- names9
  
  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = len9,byrow = TRUE))
  col11c <- col11b
  names(col11c) <- names9
  
  
  id_t <- id_1
  dfcol <- rbind(col10c)
  dfcol3 <- rbind(col11c)
  dfcol2 <- cbind(id_t,dfcol,dfcol3)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  nnc <- ncol(dfcol2)
  nnc2 <- (nnc-1) / 2 
  nnc3 <- nnc2 * 2 + 1
  names(dfcol2)[(nnc2+2):nnc3]<- paste0("Tum_",names(dfcol2)[2:(nnc2+1)])
  
  ad1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$DP4),",")),ncol = 4,byrow = TRUE))
  names(ad1) <- c("Nor_Ref_ST_p","Nor_Ref_ST_n","Nor_Alt_ST_p","Nor_Alt_ST_n")
  # gss1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$BCOUNT),",")),ncol = 4,byrow = TRUE))
  # names(gss1) <- c("BCOUNT_ref","BCOUNT_alt","BCOUNT_3","BCOUNT_4")
  
  ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$Tum_DP4),",")),ncol = 4,byrow = TRUE))
  names(ad1b) <- c("Tum_Ref_ST_p","Tum_Ref_ST_n","Tum_Alt_ST_p","Tum_Alt_ST_n")
  # gss1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$BCOUNT_2),",")),ncol = 4,byrow = TRUE))
  # names(gss1b) <- c("BCOUNT_ref_2","BCOUNT_alt_2","BCOUNT_2_3","BCOUNT_2_4")
  
  
  # FINAL FORMAT 
  formatfin <- cbind(dfcol2[,c(1:7)],ad1,dfcol2[,c(9:14)],ad1b)
  
  ### INFO
  f8 <- vcffile$V8
  
  library(tidyr)
  s1 <- separate(data = vcffile, col = V8, into = paste0("a",1:47), sep = ";")
  s2 <- s1[,c(8,10:54)]
  # s2 <- s2[,-2]
  snames <- unlist(strsplit(as.character(s2[1,]),"="))
  snames2 <- snames[c(seq(1,93,by = 2))]
  
  azkeNA <- which(is.na(snames2))[1]-1
  for(h in c(1:azkeNA))
  {
    print(h)
    if(length(grep("GERP\\+\\+_RS", snames2[h]))>0) 
    {
      s2[,h] <- gsub("GERP\\+\\+_RS","GERP_RS",s2[,h])
      snames2[h] <- gsub("GERP\\+\\+_RS","GERP_RS",snames2[h])
    }
    
    s2[,h] <- gsub(paste0(snames2[h],"="),"",s2[,h])
    names(s2)[h] <- snames2[h]
  }
  
  kendu <- sapply(s2,function(x) 
  {
    nr <- length(x)
    sr <- sum(is.na(x))
    nr!=sr
  })
  
  tt <- table(kendu)
  if(length(tt)>1) s2 <- s2[,kendu]
  
  names(s2)[which(names(s2)=="DP")] <- "DP_Tum_Nor"
  
  # names(s2)[42] <- "GERP"
  # s2[,37] <- gsub(paste0("GERP++_RS=.","",s2[,h+1]))
  dftot <- cbind(vcffile[,1:7],s2,formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")
  
  names(dftot)[which(names(dftot)=="RD")] <- "Nor_Ref_D"
  names(dftot)[which(names(dftot)=="DP")] <- "Nor_DP"
  names(dftot)[which(names(dftot)=="AD")] <- "Nor_Alt_D"
  
  names(dftot)[which(names(dftot)=="Tum_RD")] <- "Tum_Ref_D"
  names(dftot)[which(names(dftot)=="Tum_AD")] <- "Tum_Alt_D"
  
  w_names <- which(names(dftot)%in%c("Pos","DP_Tum_Nor","SSC","GPV","SPV","id_t","Nor_DP","Nor_Ref_D","Nor_Alt_D","Nor_Ref_ST_p" ,"Nor_Ref_ST_n", "Nor_Alt_ST_p", "Nor_Alt_ST_n",
                                     "Tum_DP","Tum_Ref_D","Tum_Alt_D","Tum_Ref_ST_p","Tum_Ref_ST_n","Tum_Alt_ST_p","Tum_Alt_ST_n"))
  
  for(z in w_names)
  {
    dftot[,z] <- as.numeric(as.character(dftot[,z]))
  }
  
  dftot$Tum_Ratio_Alt_st <- dftot$Tum_Alt_ST_p/dftot$Tum_Alt_ST_n
  
  dftot$FREQ <- as.character(dftot$FREQ)
  dftot$FREQ <- as.numeric(gsub("%","",dftot$FREQ))
  dftot$Tum_FREQ <- as.character(dftot$Tum_FREQ)
  dftot$Tum_FREQ <- as.numeric(gsub("%","",dftot$Tum_FREQ))
  
  dftot$Tum_FREQ <- dftot$Tum_FREQ/100
  dftot$FREQ <- dftot$FREQ/100
  
  names(dftot)[which(names(dftot)=="FREQ")] <- "Nor_VAF"
  names(dftot)[which(names(dftot)=="Tum_FREQ")] <- "Tum_VAF"
  
  
  return(dftot)
  
}

radia_2_excel <- function(vcffile, txtfile)
{
  # vcffile <- radia
  # txtfile <- radiatxt
  
  # separate format
  vcffile$V9 <- as.character(vcffile$V9)
  vcffile$ID <- 1:(nrow(vcffile))
  
  # Extract names
  u1 <- unique(vcffile$V9)
  u1 <- u1[!is.na(u1)] 
  names9 <- strsplit(u1,":")[[1]]
  len9 <- length(names9)
  
  # GROUP 1 
  g1_10 <- as.character(vcffile$V10)
  g1_11 <- as.character(vcffile$V11)
  g1_12 <- as.character(vcffile$V12)
  id_1 <- vcffile$ID
  
  col10 <- unlist(strsplit(g1_10,":"))
  col10b <- data.frame(matrix(col10,ncol = len9,byrow = TRUE))
  col10c <- col10b
  names(col10c) <- names9
  
  col11 <- unlist(strsplit(g1_11,":"))
  col11b <- data.frame(matrix(col11,ncol = len9,byrow = TRUE))
  col11c <- col11b
  names(col11c) <- names9
  
  col12 <- unlist(strsplit(g1_12,":"))
  col12b <- data.frame(matrix(col12,ncol = len9,byrow = TRUE))
  col12c <- col12b
  names(col12c) <- names9
  
  id_t <- id_1
  dfcol <- rbind(col10c)
  dfcol3 <- rbind(col11c)
  dfcol4 <- rbind(col12c)
  dfcol2 <- cbind(id_t,dfcol,dfcol3,dfcol4)
  dfcol2 <- dfcol2[order(dfcol2$id_t),]
  nnc <- ncol(dfcol2)
  nnc2 <- (nnc-1) / 3 
  nnc3 <- nnc2 * 2 + 1
  names(dfcol2)[(nnc2+2):nnc3]<- paste0(names(dfcol2)[2:(nnc2+1)],"_2")
  nnc2b <- nnc3 + 1
  nnc3b <- nnc2b + nnc2 - 1
  names(dfcol2)[(nnc2b):nnc3b]<- paste0(names(dfcol2)[2:(nnc2+1)],"_3")
  
  l <- strsplit(as.character(dfcol2$AD),",")
  n.obs <- sapply(l, length)
  seq.max <- seq_len(max(n.obs))
  ad1 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  names(ad1) <- c("AD_ref","AD_alt","AD_3","AD_4")
  
  # l <- strsplit(as.character(dfcol2$AF),",")
  # n.obs <- sapply(l, length)
  # seq.max <- seq_len(max(n.obs))
  # ad1 <- data.frame(t(sapply(l, "[", i = seq.max))[,1:4])
  # names(ad1) <- c("AD_ref","AD_alt","AD_3","AD_4")
  
  # af1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AF),",")),ncol = 4,byrow = TRUE))
  af1 <- data.frame(t(sapply(strsplit(as.character(dfcol2$AF),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$AF),","), length)))))[,1:4])
  names(af1) <- c("AF_ref","AF_alt","AF_3","AF_4")
  # dp1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$DP4),",")),ncol = 4,byrow = TRUE))
  dp1 <- data.frame(t(sapply(strsplit(as.character(dfcol2$DP4),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$DP4),","), length)))))[,1:4])
  names(dp1) <- c("DP4_ref","DP4_alt","DP4_3","DP4_4")
  # mq01 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$MQ0),",")),ncol = 4,byrow = TRUE))
  mq01 <- data.frame(t(sapply(strsplit(as.character(dfcol2$MQ0),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$MQ0),","), length)))))[,1:4])
  names(mq01) <- c("MQ0_ref","MQ0_alt","MQ0_3","MQ0_4")
  # mmq1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$MMQ),",")),ncol = 4,byrow = TRUE))
  mmq1 <- data.frame(t(sapply(strsplit(as.character(dfcol2$MMQ),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$MMQ),","), length)))))[,1:4])
  names(mmq1) <- c("MMQ_ref","MMQ_alt","MMQ_3","MMQ_4")
  # mqa1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$MQA),",")),ncol = 4,byrow = TRUE))
  mqa1 <- data.frame(t(sapply(strsplit(as.character(dfcol2$MQA),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$MQA),","), length)))))[,1:4])
  names(mqa1) <- c("MQA_ref","MQA_alt","MQA_3","MQA_4")
  # bq1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$BQ),",")),ncol = 4,byrow = TRUE))
  bq1 <- data.frame(t(sapply(strsplit(as.character(dfcol2$BQ),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$BQ),","), length)))))[,1:4])
  names(bq1) <- c("BQ_ref","BQ_alt","BQ_3","BQ_4")
  # sb1 <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$SB),",")),ncol = 4,byrow = TRUE))
  sb1 <- data.frame(t(sapply(strsplit(as.character(dfcol2$SB),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$SB),","), length)))))[,1:4])
  names(sb1) <- c("SB_ref","SB_alt","SB_3","SB_4")
  
  # ad1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AD_2),",")),ncol = 4,byrow = TRUE))
  ad1b <- data.frame(t(sapply(strsplit(as.character(dfcol2$AD_2),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$AD_2),","), length)))))[,1:4])
  names(ad1b) <- c("AD_ref_2","AD_alt_2","AD_3_2","AD_4_2")
  # af1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AF_2),",")),ncol = 4,byrow = TRUE))
  af1b <- data.frame(t(sapply(strsplit(as.character(dfcol2$AF_2),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$AF_2),","), length)))))[,1:4])
  names(af1b) <- c("AF_ref_2","AF_alt_2","AF_3_2","AF_4_2")
  # dp1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$DP4_2),",")),ncol = 4,byrow = TRUE))
  dp1b <- data.frame(t(sapply(strsplit(as.character(dfcol2$DP4_2),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$DP4_2),","), length)))))[,1:4])
  names(dp1b) <- c("DP4_ref_2","DP4_alt_2","DP4_3_2","DP4_4_2")
  # mq01b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$MQ0_2),",")),ncol = 4,byrow = TRUE))
  mq01b <- data.frame(t(sapply(strsplit(as.character(dfcol2$MQ0_2),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$MQ0_2),","), length)))))[,1:4])
  names(mq01b) <- c("MQ0_ref_2","MQ0_alt_2","MQ0_3_2","MQ0_4_2")
  # mmq1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$MMQ_2),",")),ncol = 4,byrow = TRUE))
  mmq1b <- data.frame(t(sapply(strsplit(as.character(dfcol2$MMQ_2),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$MMQ_2),","), length)))))[,1:4])
  names(mmq1b) <- c("MMQ_ref_2","MMQ_alt_2","MMQ_3_2","MMQ_4_2")
  # mqa1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$MQA_2),",")),ncol = 4,byrow = TRUE))
  mqa1b <- data.frame(t(sapply(strsplit(as.character(dfcol2$MQA_2),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$MQA_2),","), length)))))[,1:4])
  names(mqa1b) <- c("MQA_ref_2","MQA_alt_2","MQA_3_2","MQA_4_2")
  # bq1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$BQ_2),",")),ncol = 4,byrow = TRUE))
  bq1b <- data.frame(t(sapply(strsplit(as.character(dfcol2$BQ_2),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$BQ_2),","), length)))))[,1:4])
  names(bq1b) <- c("BQ_ref_2","BQ_alt_2","BQ_3_2","BQ_4_2")
  # sb1b <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$SB_2),",")),ncol = 4,byrow = TRUE))
  sb1b <- data.frame(t(sapply(strsplit(as.character(dfcol2$SB_2),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$SB_2),","), length)))))[,1:4])
  names(sb1b) <- c("SB_ref_2","SB_alt_2","SB_3_2","SB_4_2")
  
  # ad1c <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AD_3),",")),ncol = 4,byrow = TRUE))
  ad1c <- data.frame(t(sapply(strsplit(as.character(dfcol2$AD_3),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$AD_3),","), length)))))[,1:4])
  names(ad1c) <- c("AD_ref_3","AD_alt_3","AD_3_3","AD_4_3")
  # af1c <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$AF_3),",")),ncol = 4,byrow = TRUE))
  af1c <- data.frame(t(sapply(strsplit(as.character(dfcol2$AF_3),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$AF_3),","), length)))))[,1:4])
  names(af1c) <- c("AF_ref_3","AF_alt_3","AF_3_3","AF_4_3")
  # dp1c <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$DP4_3),",")),ncol = 4,byrow = TRUE))
  dp1c <- data.frame(t(sapply(strsplit(as.character(dfcol2$DP4_3),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$DP4_3),","), length)))))[,1:4])
  names(dp1c) <- c("DP4_ref_3","DP4_alt_3","DP4_3_3","DP4_4_3")
  # mq01c <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$MQ0_3),",")),ncol = 4,byrow = TRUE))
  mq01c <- data.frame(t(sapply(strsplit(as.character(dfcol2$MQ0_3),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$MQ0_3),","), length)))))[,1:4])
  names(mq01c) <- c("MQ0_ref_3","MQ0_alt_3","MQ0_3_3","MQ0_4_3")
  # mmq1c <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$MMQ_3),",")),ncol = 4,byrow = TRUE))
  mmq1c <- data.frame(t(sapply(strsplit(as.character(dfcol2$MMQ_3),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$MMQ_3),","), length)))))[,1:4])
  names(mmq1c) <- c("MMQ_ref_3","MMQ_alt_3","MMQ_3_3","MMQ_4_3")
  # mqa1c <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$MQA_3),",")),ncol = 4,byrow = TRUE))
  mqa1c <- data.frame(t(sapply(strsplit(as.character(dfcol2$MQA_3),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$MQA_3),","), length)))))[,1:4])
  names(mqa1c) <- c("MQA_ref_3","MQA_alt_3","MQA_3_3","MQA_4_3")
  # bq1c <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$BQ_3),",")),ncol = 4,byrow = TRUE))
  bq1c <- data.frame(t(sapply(strsplit(as.character(dfcol2$BQ_3),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$BQ_3),","), length)))))[,1:4])
  names(bq1c) <- c("BQ_ref_3","BQ_alt_3","BQ_3_3","BQ_4_3")
  # sb1c <- data.frame(matrix(unlist(strsplit(as.character(dfcol2$SB_3),",")),ncol = 4,byrow = TRUE))
  sb1c <- data.frame(t(sapply(strsplit(as.character(dfcol2$SB_3),","), "[", i = seq_len(max(sapply(strsplit(as.character(dfcol2$SB_3),","), length)))))[,1:4])
  names(sb1c) <- c("SB_ref_3","SB_alt_3","SB_3_3","SB_4_3")
  
  
  # FINAL FORMAT 
  formatfin <- cbind(dfcol2[,c(1:3)],ad1,af1,
                     dfcol2[,c(6:7)],dp1,
                     dfcol2[,c(9:10)],mq01,mmq1,mqa1,bq1,sb1,
                     
                     dfcol2[,c(16:18)],ad1b,af1b,
                     dfcol2[,c(21:22)],dp1b,
                     dfcol2[,c(24:25)],mq01b,mmq1b,mqa1b,bq1b,sb1b,
                     
                     dfcol2[,c(31:33)],ad1c,af1c,
                     dfcol2[,c(36:37)],dp1c,
                     dfcol2[,c(39:40)],mq01c,mmq1c,mqa1c,bq1c,sb1c,
                     dfcol2[,c(46)])
  names(formatfin)[ncol(formatfin)] <- "MMP_3"
  
  dftot <- cbind(vcffile[,1:7],formatfin)
  names(dftot)[1:7] <- c("Chr","Pos","V3","Ref","Alt","V6","Filt")
  
  mut41_txt2 <- txtfile[,1:10]
  mut41_txt2$merge <- paste0(mut41_txt2$Chr,mut41_txt2$Start)
  dftot$merge <- paste0(dftot$Chr,dftot$Pos)
  
  dftot2 <- merge(mut41_txt2,dftot, by = c("merge"))
  dftotx <- merge(mut41_txt2,dftot, by = c("merge"), all.x = TRUE)
  dftoty <- merge(mut41_txt2,dftot, by = c("merge"), all.y = TRUE)
  
  return(dftot2)
  
}




# functions

f1_load_missings <- function(raw_data = data, n_lags = n_lags, 
                             relevant_variables = NULL, 
                             ozup = ozup, ozlow = ozlow,
                             int_var, lag_var2, date_var)
{
  
  print(int_var)
  print(lag_var2)
  
  
  
  
  lags_var <- c(int_var,lag_var2)
  
  processed_data = raw_data[,lags_var]
  
  nb_lag = n_lags
  for (j in 1:ncol(processed_data)){
    for (l in 1:nb_lag){
      new_col_name = paste(colnames(processed_data)[j],"_lag_",toString(l),sep = "")
      processed_data = cbind(processed_data, c(rep(NA,l),processed_data[[j]][1:(nrow(processed_data)-l)]))
      colnames(processed_data)[ncol(processed_data)] = new_col_name
    }
  }
  # 
  processed_data <- cbind(raw_data[,-which(names(raw_data)%in%lags_var)],processed_data)
  processed_data$date <- processed_data[,date_var]
  
  #------------------------------------------------------------------------#
  # add explicit covariates for day of week, month, year (coded as integers)
  month = factor(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                 levels = c("Jan","Feb","Mar","Apr","May","Jun", "Jul","Aug","Sep","Oct","Nov","Dec"))
  processed_data$month = month[as.numeric(format(processed_data$date,"%m"))]
  processed_data$year = as.numeric(format(processed_data$date,"%Y"))
  # construct ozone high/low indicator variable
  ozone_threshold_upper = ozup
  ozone_threshold_lower = ozlow
  # ozone_threshold = 26.42824 # WARNING: this value was arbitrarily set to the empirical 0.75-quantile !!!
  # processed_data$is_highO3 = processed_data$o3 > ozone_threshold
  processed_data$is_highO3 = rep(NA,nrow(processed_data))
  processed_data$is_highO3[processed_data[,int_var] > ozone_threshold_upper] = TRUE
  processed_data$is_highO3[processed_data[,int_var] < ozone_threshold_lower] = FALSE
  
  # print("3")
  # print(table(processed_data$is_highO3))
  
  # construct treatment assigment variable
  lag_o3 = n_lags # number of preceding days with low ozone levels (must be at least 1)
  processed_data$is_treated = rep(NA,nrow(processed_data))
  for (l in (lag_o3+1):nrow(processed_data)){
    low_ozone_on_previous_days = all(!processed_data$is_highO3[(l-lag_o3):(l-1)])
    is_admissible = low_ozone_on_previous_days & (!is.na(low_ozone_on_previous_days))
    # is_admissible = all(!processed_data$is_highO3[(l-lag_o3):(l-1)])
    if(is_admissible){
      processed_data$is_treated[l] = processed_data$is_highO3[l]
    }
  }
  # remove the first lag_o3 days (for which the treatment assigment is not well-defined)
  lag_max = max(lag_o3,nb_lag)
  processed_data = processed_data[(lag_max+1):nrow(processed_data),]
  # remove the days for which assignment is undefined
  matching_data = processed_data[!is.na(processed_data$is_treated),]
  N = nrow(matching_data)
  #--------------------------------------------------------------------------------------------------------------#
  #--------------------------------------------------------------------------------------------------------------#
  
  #--------- explore treated units ---------#
  treated_units = subset(matching_data,is_treated)
  control_units = subset(matching_data,!is_treated)
  
  # print(dim(matching_data))
  # print("4")
  
  return(list(treated_units = treated_units, control_units = control_units, matching_data = matching_data, processed_data = processed_data))
}


f2_thresholds <- function(matching_data, n_lags = n_lags, treated_units, control_units, thresholds2, processed_data = processed_data, 
                          showmatched = TRUE, showmatched_detail = TRUE, remove = FALSE, how_many_remove = 0)
{
  scaling =  rep(list(1),ncol(matching_data))
  names(scaling) = colnames(matching_data)
  
  # browser()
  
  discrepancies = discrepancyMatrix(treated_units, control_units, thresholds2, scaling)
  
  # print("14411")
  cat("Number of prospective matched treated units =",sum(rowSums(discrepancies<Inf)>0),"\n")
  # Perform matching
  rownames(discrepancies) = format(matching_data$date[which(matching_data$is_treated)],"%Y-%m-%d")
  # print("1412")
  colnames(discrepancies) = format(matching_data$date[which(!matching_data$is_treated)],"%Y-%m-%d")
  # print("14")
  rownames(matching_data) = matching_data$date
  
  # --------------------------------------------------------------------------------------------------------------#
  # print("142")
  matched_groups = fullmatch(discrepancies,data=matching_data,remove.unmatchables = TRUE,mean.controls=1)
  # Get list of groups
  # print("143")
  groups_labels = unique(matched_groups[!is.na(matched_groups)])
  groups_list = list()
  # print("144")
  for (i in 1:length(groups_labels)){
    IDs = names(matched_groups)[(matched_groups==groups_labels[i])]
    groups_list[[i]] = as.Date(IDs[!is.na(IDs)])
  }
  # Print diagnostics
  # print(matched_groups,grouped=TRUE)
  
  N_treated = nrow(treated_units)
  N_control = nrow(control_units)
  
  adj = (discrepancies<Inf)
  edges_mat = which(adj,arr.ind = TRUE)
  weights = 1/(1+sapply(1:nrow(edges_mat),function(i)discrepancies[edges_mat[i,1],edges_mat[i,2]]))
  edges_mat[,"col"] = edges_mat[,"col"] + N_treated
  edges_vector = c(t(edges_mat))
  BG = make_bipartite_graph(c(rep(TRUE,N_treated),rep(FALSE,N_control)), edges = edges_vector)
  MBM = maximum.bipartite.matching(BG,weights = weights)
  pairs_list = list()
  N_matched = 0
  for (i in 1:N_treated){
    if (!is.na(MBM$matching[i])){
      N_matched = N_matched + 1
      pairs_list[[N_matched]] = c(treated_units$date[i],control_units$date[MBM$matching[i]-N_treated])
    }
  }
  cat("Number of matched treated units =", N_matched,"\n")
  # print("16")
  total <- processed_data[1,]
  relevant_fields = colnames(matching_data)[which(unlist(thresholds2)<Inf)]
  for (i in 1:N_matched){
    t_i = toString(subset(matching_data,date%in%(pairs_list[[i]])&is_treated)$date)
    c_i = toString(subset(matching_data,date%in%(pairs_list[[i]])&!is_treated)$date)
    if(showmatched)
    {
      cat("\n-------------------- Matched pair",i,":",t_i,",",c_i,"--------------------\n")
      if(showmatched_detail)
      {
        print(subset(matching_data,date%in%(pairs_list[[i]]))[relevant_fields])
      }
      cat("...........\n")
    }
    
    for (j in 1:length(pairs_list[[i]])){
      # print(subset(processed_data,date%in%(pairs_list[[i]][j]-0:lag_o3))[relevant_fields])
      # print(n_lags)
      total <- rbind(total,subset(processed_data,date%in%(pairs_list[[i]][j]-0:n_lags)))
    }
  }
  total <- total[-1,]
  if(remove)
  {
    total <- total[-c(1:how_many_remove),]
  }
  
  return(list(N_matched = N_matched, pairs_list = pairs_list, total = data.frame(total), thresholds = thresholds2))
}


plot_sen <- function(min1,max1,tarte,punto,
                     city = var_city, n_lags = n_lags, 
                     relevant_variables, thresholds2 = thresholds, mort = FALSE)
{
  sek <- seq(min1,max1,tarte)
  um <- punto
  up1 <- punto+sek
  low1 <- punto-sek
  lend <- length(sek)
  
  min2 <- min(low1)
  max2 <- max(up1)
  
  df1 <- data.frame(matrix(NA,lend,6))
  names(df1) <- c("Down","Up","Matches", "Mortcont", "MortTreat", "MortRatio")
  
  for(ii in 1:lend)
  {
    print(paste0("ii - ", ii))
    res1 <- res1 <- f1_load_missings(city = city, n_lags = n_lags, 
                                     relevant_variables = NULL, ozup=up1[ii], ozlow=low1[ii])
    
    # names(res1)
    treated_units = res1$treated_units
    control_units = res1$control_units
    matching_data = res1$matching_data
    processed_data = res1$processed_data
    
    # print("res2aurretik")
    
    res2 <- f2_thresholds (matching_data, n_lags = n_lags,treated_units = treated_units, control_units = control_units, thresholds2 = NULL, scaling2 = NULL, processed_data = processed_data)
    
    # print("res2ondoren")
    
    zen <- res2$N_matched
    
    # print(class(zen))
    
    if(zen[1]<5) break
    df1[ii,1]<-low1[ii]
    df1[ii,2]<-up1[ii]
    df1[ii,3]<-zen
    df1[ii,4]<-sum(control_units$death)
    df1[ii,5]<-sum(treated_units$death)
  }
  
  # print("loop kanpoan")
  
  df1 <- df1[!is.na(df1[,1]),1:6]
  df1$MortRatio <- df1$MortTreat / df1$Mortcont
  df1$id <- 1:dim(df1)[1]
  yylim <- c(min(df1$Down),max(df1$Up))
  max3 <- max(df1$Matches)
  
  par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
  plot(df1$id,df1$Up,type="l",ylim=yylim,ylab="O3 concentration",xlab="",
       col=2,pch=15,bg="grey",lty=3,lwd=2,main=paste0("Num of matches:",max3, " - Threshold limit-",punto)) # first plot
  lines(df1$id,df1$Down,type="l",col=2,lty=3,lwd=2) # first plot
  abline(h=seq(min2,max2,1),col="grey")
  
  if(mort)
  {
    par(new = TRUE)
    plot(df1$id,df1$MortRatio, type = "l", lwd=3,
         axes = FALSE, xlab = "", ylab = "",col="green")
    axis(side=4, at = pretty(range(df1$MortRatio)))
    mtext("Mortality proportion", side=4, line=3)
  } else {
    par(new = TRUE)
    plot(df1$id,df1$Matches, type = "l", lwd=3,
         axes = FALSE, xlab = "", ylab = "",col="blue")
    axis(side=4, at = pretty(range(df1$Matches)))
    mtext("Number of matches", side=4, line=3)
  }
  
  return(df1)
}

plot_other_pol <- function(tot, var = c("o3", "tmpd", "no2", "so2"), zenb = 1, abl = c(29,30))
{
  
  par(mfrow=c(2,2))
  
  varint <- var[1]
  ylim2 <- range(tot[1:(zenb*2),varint])
  ylim2[1] <- ifelse(ylim2[1]>0, 0, ylim2[1])
  # print(ylim2)
  plot(1:zenb, tot[1:zenb,varint],type="b", col = 2, ylab = varint, lwd = 2, ylim = ylim2, main = varint, xlab = "Days")
  lines(1:zenb,tot[(zenb+1):(zenb*2),varint],col=4,type="b",lwd=2)
  abline(h=abl[1],lty=3)
  abline(h=abl[2],lty=3)
  
  varint <- var[2]
  ylim2 <- range(tot[1:(zenb*2),varint])
  ylim2[1] <- ifelse(ylim2[1]>0, 0, ylim2[1])
  plot(1:zenb,tot[1:zenb,varint],type="b",col=2, ylab = varint, lwd = 2, ylim = ylim2, main = varint, xlab = "Days")
  lines(1:zenb,tot[(zenb+1):(zenb*2),varint],col=4,type="b",lwd=2)
  # abline(h=abl,lty=3)
  
  varint <- var[3]
  ylim2 <- range(tot[1:(zenb*2),varint])
  ylim2[1] <- ifelse(ylim2[1]>0, 0, ylim2[1])
  plot(1:zenb,tot[1:zenb,varint],type="b",col=2, ylab = varint, lwd = 2, ylim = ylim2, main = varint, xlab = "Days")
  lines(1:zenb,tot[(zenb+1):(zenb*2),varint],col=4,type="b",lwd=2)
  
  varint <- var[4]
  ylim2 <- range(tot[1:(zenb*2),varint])
  ylim2[1] <- ifelse(ylim2[1]>0, 0, ylim2[1])
  plot(1:zenb,tot[1:zenb,varint],type="b",col=2, ylab = varint, lwd = 2, ylim = ylim2, main = varint, xlab = "Days")
  lines(1:zenb,tot[(zenb+1):(zenb*2),varint],col=4,type="b",lwd=2)
  
  par(mfrow=c(1,1))
}






library(Rcpp)
# pairwise difference between real-univariate covariate of treated VS control group
cppFunction('NumericMatrix pairDistCpp(NumericVector treated, NumericVector control) {
            NumericMatrix D(treated.size(), control.size());
            for (int t = 0; t < treated.size(); t++) {
            for (int c = 0; c < control.size(); c++) {
            D(t,c) = abs(treated[t] - control[c]);
            }
            }
            return D;
            }')
# pairwise absolute difference between real-univariate covariates of treated VS control group
cppFunction('NumericMatrix pairAbsDistCpp(NumericVector treated, NumericVector control) {
            NumericMatrix D(treated.size(), control.size());
            for (int t = 0; t < treated.size(); t++) {
            for (int c = 0; c < control.size(); c++) {
            D(t,c) = abs(treated[t] - control[c]);
            }
            }
            return D;
            }')
# pairwise difference between factor-valued (i.e. bounded integer-valued) covariates 
# (e.g. day of the week, month, ...) of treated VS control group, assuming the facotr levels are cyclic
# and only the shortest difference modulo nb_levels matters.
pairModuloDist = function(factors_treated, factors_control, nb_levels) {
  return (pmin(pairDistCpp(as.integer(factors_treated),as.integer(factors_control))%%nb_levels,
               t(pairDistCpp(as.integer(factors_control),as.integer(factors_treated))%%nb_levels)))
}
# pairwise difference between covariates of treated VS control group
# Inputs: treated/control are of covariate vectors (one entry per unit, for a given covariate)
# Outputs: pairwise difference matrix
pairdifference = function(treated, control){
  if(is.factor(treated[1])){
    # if factor-valued, use shortest difference modulo number of levels
    return (pairModuloDist(treated, control, length(levels(treated[1]))))
  } else {
    return (pairAbsDistCpp(treated,control))
  }
}
# pairwise discrepancy between treated VS control group
# Inputs: treated/control are lists or dataframes (with one column per covariate, one row per unit), 
# thresholds is a LIST of values for each covariate to be matched (match is admissible
# if and only if the differences between covariates are all less than the associated thresholds), 
# standard_deviations is an optional LIST of values to standardize/reweight the differences
discrepancyMatrix = function(treated, control, thresholds, scaling = NULL){
  nb_covariates = ncol(treated)
  # matrix of pairwise discrepancies (computed as standardized L1-distance)
  D = matrix(0, nrow = nrow(treated), ncol = nrow(control))
  # keep track of pairs that are non-admissible matches
  non_admissible = matrix(FALSE, nrow = nrow(treated), ncol = nrow(control))
  for (i in 1:nb_covariates){
    # only compute the distances for covariates that are matched on (i.e. finite thresholds)
    # print(names(treated[i]))
    if (thresholds[[i]]<Inf){
      # print("dentro inf")
      differences = pairdifference(treated[[i]], control[[i]])
      # print(differences)
      if(is.null(scaling))
      {
        D = D + differences
      } else {
        D = D + differences*scaling[[i]]
      }
      # print(D)
      # The user is responsible for inputing complete data (i.e. impute missing data beforhand if needed).
      # In the undesirable case where some covariates that are matched on are NA, we exclude the corresponding
      # unit from the matching (default behavior for convenience, but NOT for statistical validity, especially
      # if the missing-data mechanism is non-ignorable)
      differences[is.na(differences)] = Inf
      non_admissible = non_admissible | (differences > thresholds[[i]])
    }
  }
  D = D/nb_covariates # "standardize" the discrepancies (just for convenience, doesn't change the matching at all)
  D[non_admissible] = Inf # give infinite penalty to non-admissible pairs
  return (D)
}

# missings per year
# name of the date column shold be date
missing_per_year <- function(dataset)
{
  day1 <- as.Date(substr(dataset$date[1],1,10), format="%Y-%m-%d")
  day_last<-  as.Date(substr(dataset$date[nrow(dataset)],1,10), format="%Y-%m-%d")
  
  all_days <- data.frame(seq(day1, day_last, "days"))
  names(all_days) <- "AllDays"
  all_years2 <- substr(all_days[,1],1,4)
  all_days2 <- data.frame(all_years2, all_days)
  
  all_years <- substr(seq(day1, day_last, "years"),1,4)
  
  
  dates <- data.frame(as.Date(substr(dataset$date,1,10), format="%Y-%m-%d"))
  dates2 <- data.frame(dates,dates)
  names(dates2) <- c("AllDays", "WithDays")
  
  m <- merge(all_days2, dates2, all.x = TRUE, by= "AllDays")
  m3 <- NULL
  for( y in all_years)
  {
    m2 <- m[m$all_years2==y,]
    m3 <- c(m3,sum(is.na(m2$WithDays)))
  }
  print(m3)
  dfm <- data.frame(all_years,m3)
  names(dfm) <- c("Year", "Missing days")
  return(dfm)
}


# plot patterns

plotone <- function(tot, lags0, varint, abl = c(25,32), int = "o3")
{
  dat <- data.frame(rep(1:lags0,2),
                    c(tot[1:lags0,varint],tot[(lags0+1):(lags0*2),varint]),
                    c(rep("Treatment",4), rep("Control",4)))
  
  names(dat) <- c("Lags","Concentration","Group")
  
  htot <- NA
  lowestp <- NA
  for(h in 1:4)
  {
    htot <- c(htot, round(dat[h,2]-dat[h+lags0,2],2))
    lowestp <- c(lowestp, round(ifelse(dat[h,2]>dat[h+lags0,2],dat[h+lags0,2]-dat[h+lags0,2]*0.2,dat[h,2]-dat[h,2]*0.2),1))
  }
  dat$dif <- htot[-1]
  dat$low <- lowestp[-1]
  
  ylim2 <- range(tot[1:(lags0*2),varint])
  ylim2[1] <- ifelse(ylim2[1]>0, 0, ylim2[1])
  
  
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # The palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # To use for fills, add
  scale_fill_manual(values=cbPalette)
  
  # To use for line and point colors, add
  scale_colour_manual(values=cbPalette)
  
  p <- ggplot(dat, aes(Lags, Concentration, color = Group)) +
    geom_line(size = 2, alpha = .8) +
    theme_minimal() +
    ggtitle(varint) +
    labs(x = "Lags", y = "Concentration") + ylim(ylim2) +
    annotate("text", x = dat$Lags[1:lags0], y = dat$low[1:lags0], label = dat$dif[1:lags0]) + theme(legend.position="none")  + 
    scale_color_manual(values=cbPalette)
  
  if(varint==int)
  {
    # browser()
    p <- p + geom_hline(yintercept=abl[1], linetype="dashed", color = "black") +
      geom_hline(yintercept=abl[2], linetype="dashed", color = "black") 
  }
  
  p
}
