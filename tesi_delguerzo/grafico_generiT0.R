#install.packages("tidyr")
#install.packages("ggplot2")
library(openxlsx)
library(tibble)
library(tidyr)
library(ggplot2)
library(dplyr)
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "T0")
mygenere <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/percentuali_generi.xlsx")
#elimino il nome della prima colonna in modo che la colonna 1 dell'excel corrisponde alla prima riga di myfragment
rownames(mygenere) <- mygenere$nome_ID
mygenere$nome_ID <- NULL
#cambio in numerico i valori
mygenere <- mutate_all(mygenere, function(x) as.numeric(as.character(x)))
#calcolo la media
mygenere$media <- rowMeans(mygenere)
#ordino i valori medi in ordine decrescente
mygenere <- mygenere[order(mygenere$media, decreasing=TRUE),]
#prendo solo le prime 10 righe
mygenere <- head(mygenere, 10)
#elimino la colonna media
mygenere$media <- NULL
mygenere <- tibble::rownames_to_column(mygenere,"nome_ID")
mygenere$nome_ID <- gsub("_[^_]*$", "", mygenere$nome_ID)
#comando per organizzare i dati in modo che ad ogni colonna corrisponda una variabile, necessario per ggplot2
mygenere2 <- mygenere %>%
  pivot_longer(names_to = "Campioni",
               values_to = "percentuali",
               cols=-nome_ID)
positions <- c(1:12)
#costruisco il grafico con ggplot
ggplot(mygenere2, aes(x = Campioni, y=percentuali, fill=nome_ID))+
  geom_bar(position = "stack", stat = "identity")+
  ggtitle ("Stacked bar plot dei generi")+
  ylab ("Abbondanza relativa (%)")+
  scale_x_discrete(limits=positions)