#install.packages ("colorhcplot")
#install.packages("pheatmap")
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
#SUOLO T0
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "T0")
rownames(mygenere) <- mygenere$nome_ID
mygenere$nome_ID <- NULL
#rendo numerico
mygenere <- mutate_all(mygenere, function(x) as.numeric(as.character(x)))
#seleziono una soglia
keep <- rowSums(mygenere)>= 0.05
mygenere <- mygenere [keep,]
#traspongo
mygenere_trasp <- as.data.frame(t(mygenere))
#calcolo la distanza
dist_genere <- dist(mygenere_trasp)
hcgenere <- hclust(dist_genere)
plot(hcgenere)
fact <- as.factor(mymetadata$Soil)
colorhcplot(hcgenere, fact, color = c("black", "red"))

#SUOLO T1
mymetadataT1 <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "Campionate")
ionomics <- mymetadataT1[,8:17]
sc_ionomics <- scale(ionomics)
dist_ionomics <- dist(sc_ionomics)
hcionomics <- hclust(dist_ionomics)
plot(hcionomics)
fact <- as.factor(mymetadataT1$Soil)
colorhcplot(hcionomics, fact, color = c("black", "red", "gold"))

#STERILIZZAZIONE T0
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "T0")
rownames(mygenere) <- mygenere$nome_ID
mygenere$nome_ID <- NULL
#rendo numerico
mygenere <- mutate_all(mygenere, function(x) as.numeric(as.character(x)))
#seleziono una soglia
keep <- rowSums(mygenere)>= 0.05
mygenere <- mygenere [keep,]
#traspongo
mygenere_trasp <- as.data.frame(t(mygenere))
#calcolo la distanza
dist_genere <- dist(mygenere_trasp)
hcgenere <- hclust(dist_genere)
plot(hcgenere)
fact <- as.factor(mymetadata$Ster_soil)
colorhcplot(hcgenere, fact, color = c("red", "green"))

#STERILIZZAZIONE T1
mymetadataT1 <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "Campionate")
ionomics <- mymetadataT1[,8:17]
sc_ionomics <- scale(ionomics)
dist_ionomics <- dist(sc_ionomics)
hcionomics <- hclust(dist_ionomics)
plot(hcionomics)
fact <- as.factor(mymetadataT1$Ster_soil)
colorhcplot(hcionomics, fact, color = c("red", "green"))