#install.packages("tidyr")
#install.packages("ggplot2")
#install.packages("corrplot")
library(openxlsx)
library(tibble)
library(tidyr)
library(ggplot2)
library(dplyr)
library(corrplot)
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "Campionate")
mygenere <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/percentuali_generi1.xlsx")
rownames(mygenere) <- mygenere$nome_ID
mygenere$nome_ID <- NULL
#cambio in numerico i valori
mygenere <- mutate_all(mygenere, function(x) as.numeric(as.character(x)))
#calcolo la media
mygenere$media <- rowMeans(mygenere)
#ordino i valori medi in ordine decrescente
mygenere <- mygenere[order(mygenere$media, decreasing=TRUE),]
#prendo solo le prime 25 righe
mygenere <- head(mygenere, 25)
#elimino la colonna media
mygenere$media <- NULL
#elimino le parentesi quadre e l'ID
row.names(mygenere) <- gsub ("-", "_", row.names(mygenere))
row.names(mygenere) <- gsub (" ", "_", row.names(mygenere))
row.names(mygenere) <- gsub ("\\[", "", row.names(mygenere))
row.names(mygenere) <- gsub ("\\]", "", row.names(mygenere))
row.names(mygenere) <- gsub("_[^_]*$", "", row.names(mygenere))
#cambio le righe con le colonne e viceversa
tmygenere <- data.frame(t(mygenere))
#controllo se tmygenere è un dataframe
is.data.frame(tmygenere)
#cerco le prime 4 colonne e le prime 4 righe
tmygenere[1:4, 1:4]
metagenere <- merge (mymetadata, tmygenere, by.y = "row.names", by.x = "Sample_name")
metagenere[1:4, 1:4]
head(metagenere$Spathaspora)
mygenere[1:4, 1:4]
#lista delle proprieta
myprop<-c("Ca","Cu","Fe","K","Mg","Mn","Na","P","S","Zn","pH","Conductivity.uS/cm", "NPOC.mg/g.suolo", "TN.mg/g.suolo", "nmol.ATP.g-1.suolo")
#lista delle specie
mynames<-row.names(mygenere)
#funzione per generare tutti i possibili confronti tra specie e proprieta
cdf<-data.frame(expand.grid(mynames,myprop),stringsAsFactors=F)
names(cdf)<-c("nome", "proprieta")
#trasformo in caratteri
cdf$nome <- as.character(cdf$nome)
cdf$proprieta <- as.character(cdf$proprieta)
#loop per generare tutte le possibili correlazioni
for (i in 1:nrow(cdf)) {
  #test di correlazione
  mycor <- cor.test(metagenere[,cdf$nome[i]], metagenere[,cdf$proprieta[i]], method = "spearman")
  cdf$rho[i] <- round(mycor$estimate, 10)
  cdf$pvalue[i] <- round(mycor$p.value, 10)
}
#funzione per aggiustare il pvalue
cdf$fdr <- p.adjust(cdf$pvalue, method="fdr")
#ordina il pvalue
write.xlsx (cdf, "C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/correlazioni_generiT1.xlsx", row.names=F)
#genera un dataframe necessario per il grafico successivo
corrho<-corp<-matrix(NA,nrow=length(mynames),ncol=length(myprop),dimnames=list(mynames,myprop))
#riempio il dataframe con rho e pvalue usando un loop 
for (n in 1:nrow(cdf)) {
  corrho[cdf$nome[n], cdf$proprieta[n]] <- cdf$rho[n]
  corp [cdf$nome[n], cdf$proprieta[n]] <- cdf$pvalue[n]
}
corrplot(corrho, p.mat=corp,
         diag = TRUE,
         sig.level = c(0.05),
         type = "full",
         tl.col = "black",
         tl.cex = 0.6,
         insig = "label_sig",
         pch.cex = 0.6)
