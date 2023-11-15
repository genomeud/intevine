options(scipen=999)
library(tibble)

#SPECIE
#carico dati file excel
library(openxlsx)
myspecie <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/conteggio_specie1.xlsx")
head(myspecie)
#carico dati frammenti 
library(data.table)
myfragment <- fread("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/count_fragments_ITS.txt", data.table = F)
head(myfragment)
#elimino il nome della prima colonna in modo che la colonna 1 dell'excel corrisponde alla prima riga di myfragment
rownames(myspecie) <- myspecie$nome_ID
myspecie$nome_ID <- NULL
#controllo che i due file abbiano lo stesso ordine
names(myspecie)==rownames(myfragment)
#faccio un ciclo for per reiterare l'operazione
for (i in 1:ncol(myspecie)) {
  myspecie[,i] <- (myspecie[,i]/myfragment[i,2]*100)
}
myspecie <- format(myspecie, scientific=F)
myspecie <- tibble::rownames_to_column(myspecie,"nome_ID")
write.xlsx(myspecie, paste0 ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/", "percentuali_specie1.xlsx"), row.names=F)

#GENERE
#carico dati file excel
library(openxlsx)
mygenere <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/conteggio_generi1.xlsx")
head(mygenere)
#carico dati frammenti 
library(data.table)
myfragment <- fread("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/count_fragments_ITS.txt", data.table = F)
head(myfragment)
#elimino il nome della prima colonna in modo che la colonna 1 dell'excel corrisponde alla prima riga di myfragment
rownames(mygenere) <- mygenere$nome_ID
mygenere$nome_ID <- NULL
#controllo che i due file abbiano lo stesso ordine
names(mygenere)==rownames(myfragment)
#faccio un ciclo for per reiterare l'operazione
for (i in 1:ncol(mygenere)) {
  mygenere[,i] <- (mygenere[,i]/myfragment[i,2]*100)
}
mygenere <- format(mygenere, scientific=F)
mygenere <- tibble::rownames_to_column(mygenere,"nome_ID")
write.xlsx(mygenere, paste0 ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/", "percentuali_generi1.xlsx"), row.names=F)

#FAMIGLIA
#carico dati file excel
library(openxlsx)
myfamiglia <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/conteggio_famiglie1.xlsx")
head(myfamiglia)
#carico dati frammenti 
library(data.table)
myfragment <- fread("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/count_fragments_ITS.txt", data.table = F)
head(myfragment)
#elimino il nome della prima colonna in modo che la colonna 1 dell'excel corrisponde alla prima riga di myfragment
rownames(myfamiglia) <- myfamiglia$nome_ID
myfamiglia$nome_ID <- NULL
#controllo che i due file abbiano lo stesso ordine
names(myfamiglia)==rownames(myfragment)
#faccio un ciclo for per reiterare l'operazione
for (i in 1:ncol(myfamiglia)) {
  myfamiglia[,i] <- (myfamiglia[,i]/myfragment[i,2]*100)
}
myfamiglia <- format(myfamiglia, scientific=F)
myfamiglia <- tibble::rownames_to_column(myfamiglia,"nome_ID")
write.xlsx(myfamiglia, paste0 ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/", "percentuali_famiglie1.xlsx"), row.names=F)

#PHYLUM
#carico dati file excel
library(openxlsx)
myphylum <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/conteggio_phyla1.xlsx")
head(myphylum)
#carico dati frammenti 
library(data.table)
myfragment <- fread("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/count_fragments_ITS.txt", data.table = F)
head(myfragment)
#elimino il nome della prima colonna in modo che la colonna 1 dell'excel corrisponde alla prima riga di myfragment
rownames(myphylum) <- myphylum$nome_ID
myphylum$nome_ID <- NULL
#controllo che i due file abbiano lo stesso ordine
names(myphylum)==rownames(myfragment)
#faccio un ciclo for per reiterare l'operazione
for (i in 1:ncol(myphylum)) {
  myphylum[,i] <- (myphylum[,i]/myfragment[i,2]*100)
}
myphylum <- format(myphylum, scientific=F)
myphylum <- tibble::rownames_to_column(myphylum,"nome_ID")
write.xlsx(myphylum, paste0 ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/", "percentuali_phyla1.xlsx"), row.names=F)
