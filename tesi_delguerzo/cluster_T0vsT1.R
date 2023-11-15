library(openxlsx)
library(data.table)
library(vegan)
library(ggplot2)
library(ggpubr)
library(colorhcplot)
library(tibble)
library(tidyr)
library(dplyr)
library(pheatmap)
mygenere <- read.xlsx ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/percentuali_generi.xlsx")
mygenere1 <- read.xlsx ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/percentuali_generi1.xlsx")
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "T0")
mymetadataT1 <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "Campionate")
rownames(mygenere) <- mygenere$nome_ID
mygenere$nome_ID <- NULL
colnames(mygenere) = paste0(colnames(mygenere), " ", "T0")
rownames(mygenere1) <- mygenere1$nome_ID
mygenere1$nome_ID <- NULL
#MANURE
#seleziono solo manure
mygenere_manure <- mygenere[1:6]
mygenere1_manure <- mygenere1[1:12]
manure <- merge (mygenere_manure, mygenere1_manure, by="row.names", all = TRUE)
manure [is.na(manure)]=0
rownames(manure) <- manure$Row.names
manure$Row.names <- NULL
#rendo numerico
manure <- mutate_all(manure, function(x) as.numeric(as.character(x)))
#seleziono una soglia
keep <- rowSums(manure)>= 0.05
manure <- manure [keep,]
#traspongo
manure_trasp <- as.data.frame(t(manure))
#calcolo la distanza
dist_genere <- dist(manure_trasp)
hcgenere <- hclust(dist_genere)
plot(hcgenere)
fact_manure <- c(rep("T0", 6), rep("T1", 12))
fact_manure <- as.factor(fact_manure)
colorhcplot(hcgenere, fact_manure, color = c("blue", "darkgreen"))
#PEAT
#seleziono solo peat
mygenere_peat <- mygenere[7:12]
mygenere1_peat <- mygenere1[25:36]
peat <- merge (mygenere_peat, mygenere1_peat, by="row.names", all = TRUE)
peat [is.na(peat)]=0
rownames(peat) <- peat$Row.names
peat$Row.names <- NULL
#rendo numerico
peat <- mutate_all(peat, function(x) as.numeric(as.character(x)))
#seleziono una soglia
keep <- rowSums(peat)>= 0.05
peat <- peat [keep,]
#traspongo
peat_trasp <- as.data.frame(t(peat))
#calcolo la distanza
dist_genere <- dist(peat_trasp)
hcgenere <- hclust(dist_genere)
plot(hcgenere)
fact_peat <- c(rep("T0", 6), rep("T1", 12))
fact_peat <- as.factor(fact_peat)
colorhcplot(hcgenere, fact_peat, color = c("blue", "darkgreen"))