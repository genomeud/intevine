#install.packages("tidyr")
#install.packages("ggplot2")
library(openxlsx)
library(tibble)
library(tidyr)
library(ggplot2)
library(dplyr)
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "T0")
myphylum <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/percentuali_phyla.xlsx")
#elimino il nome della prima colonna in modo che la colonna 1 dell'excel corrisponde alla prima riga di myfragment
rownames(myphylum) <- myphylum$nome_ID
myphylum$nome_ID <- NULL
#cambio in numerico i valori
myphylum <- mutate_all(myphylum, function(x) as.numeric(as.character(x)))
#calcolo la media
myphylum$media <- rowMeans(myphylum)
#ordino i valori medi in ordine decrescente
myphylum <- myphylum[order(myphylum$media, decreasing=TRUE),]
#prendo solo le prime 10 righe
myphylum <- head(myphylum, 10)
#elimino la colonna media 
myphylum$media <- NULL
myphylum <- tibble::rownames_to_column(myphylum,"nome_ID")
myphylum$nome_ID <- gsub("_[^_]*$", "", myphylum$nome_ID)
#comando per organizzare i dati in modo che ad ogni colonna corrisponda una variabile, necessario per ggplot2
myphylum2 <- myphylum %>%
  pivot_longer(names_to = "Campioni",
               values_to = "percentuali",
               cols=-nome_ID)

#costruisco il grafico con ggplot
positions <- c(1:12)
ggplot(myphylum2, aes(x=Campioni, y=percentuali, fill=nome_ID))+
  geom_bar(position = "stack", stat = "identity")+
  ggtitle ("Stacked bar plot dei phyla")+
  ylab ("Abbondanza relativa (%)")+
  scale_x_discrete(limits=positions)