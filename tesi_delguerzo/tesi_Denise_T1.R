  #T1
library(data.table)
library(openxlsx)

#SPECIE
#Carico i dati
mydata <- fread("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/T1/T1/107805_ID2646_1-1-A01-ITS_S260_L001_report.txt", data.table = F)
head(mydata)
#Cambio nome colonne
names(mydata) <- c("percentuale", "1", "frammenti", "codice", "taxon_ID", "nome")
#Seleziono le righe che hanno come codice S
myspecie <- subset(mydata, codice == "S")
#Sommo le percentuali
sum(myspecie$percentuale)
#Elimino alcune colonne
myspecie$percentuale <- myspecie$frammenti <- myspecie$codice <- NULL
head(myspecie)
#Creo nuova colonna unendo nome e taxon ID
myspecie$nome_ID <- paste(myspecie$nome, myspecie$taxon_ID, sep = "_")
myspecie$taxon_ID <- myspecie$nome <- NULL
#Cambio gli spazi con _
myspecie$nome_ID <- gsub(" ", "_", myspecie$nome_ID)
#Cambio ordine colonne
myspecie <- myspecie[,c("nome_ID", "1")]
listfile <- list.files("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/T1/T1/")
#Ciclo for che percorre tutti i file caricati di T1
for (ciclo in 2:length(listfile)) {
  mydata <- fread(paste0("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/T1/T1/", listfile[ciclo]), data.table = F)
  #Cambio nome colonne
  names(mydata) <- c("percentuale", ciclo, "frammenti", "codice", "taxon_ID", "nome")
  #Seleziono le righe che hanno come codice S
  myciclo <- subset(mydata, codice == "S")
  #Elimino alcune colonne
  myciclo$percentuale <- myciclo$frammenti <- myciclo$codice <- NULL
  #Creo nuova colonna unendo nome e taxon ID
  myciclo$nome_ID <- paste(myciclo$nome, myciclo$taxon_ID, sep = "_")
  myciclo$taxon_ID <- myciclo$nome <- NULL
  #Cambio gli spazi con _
  myciclo$nome_ID <- gsub(" ", "_", myciclo$nome_ID)
  myciclo <- myciclo[,c("nome_ID", ciclo)]
  myspecie <- merge(myspecie, myciclo, by = "nome_ID", all=T)
}
#Sostituisco NA con 0
myspecie [is.na(myspecie)]=0
write.xlsx(myspecie, paste0 ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/", "conteggio_specie1.xlsx"), row.names=F)

#GENERE
#Carico i dati
mydata <- fread("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/T1/T1/107805_ID2646_1-1-A01-ITS_S260_L001_report.txt", data.table = F)
head(mydata)
#Cambio nome colonne
names(mydata) <- c("percentuale", "1", "frammenti", "codice", "taxon_ID", "nome")
#Seleziono le righe che hanno come codice G
mygenere <- subset(mydata, codice == "G")
#Elimino alcune colonne
mygenere$percentuale <- mygenere$frammenti <- mygenere$codice <- NULL
head(mygenere)
#Creo nuova colonna unendo nome e taxon ID
mygenere$nome_ID <- paste(mygenere$nome, mygenere$taxon_ID, sep = "_")
mygenere$taxon_ID <- mygenere$nome <- NULL
#Cambio ordine colonne
mygenere <- mygenere[,c("nome_ID", "1")]
listfile <- list.files("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/T1/T1/")
#Ciclo for che percorre tutti i file caricati di T1
for (ciclo in 2:length(listfile)) {
  mydata <- fread(paste0("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/T1/T1/", listfile[ciclo]), data.table = F)
  #Cambio nome colonne
  names(mydata) <- c("percentuale", ciclo, "frammenti", "codice", "taxon_ID", "nome")
  #Seleziono le righe che hanno come codice G
  myciclo <- subset(mydata, codice == "G")
  #Elimino alcune colonne
  myciclo$percentuale <- myciclo$frammenti <- myciclo$codice <- NULL
  #Creo nuova colonna unendo nome e taxon ID
  myciclo$nome_ID <- paste(myciclo$nome, myciclo$taxon_ID, sep = "_")
  myciclo$taxon_ID <- myciclo$nome <- NULL
  #Cambio gli spazi con _
  myciclo$nome_ID <- gsub(" ", "_", myciclo$nome_ID)
  myciclo <- myciclo[,c("nome_ID", ciclo)]
  mygenere <- merge(mygenere, myciclo, by = "nome_ID", all=T)
}
#Sostituisco NA con 0
mygenere [is.na(mygenere)]=0
write.xlsx(mygenere, paste0 ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/", "conteggio_generi1.xlsx"), row.names=F)

#FAMIGLIA
#Carico i dati
mydata <- fread("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/T1/T1/107805_ID2646_1-1-A01-ITS_S260_L001_report.txt", data.table = F)
head(mydata)
#Cambio nome colonne
names(mydata) <- c("percentuale", "1", "frammenti", "codice", "taxon_ID", "nome")
#Seleziono le righe che hanno come codice F
myfamiglia <- subset(mydata, codice == "F")
#Elimino alcune colonne
myfamiglia$percentuale <- myfamiglia$frammenti <- myfamiglia$codice <- NULL
#Creo nuova colonna unendo nome e taxon ID
myfamiglia$nome_ID <- paste(myfamiglia$nome, myfamiglia$taxon_ID, sep = "_")
myfamiglia$taxon_ID <- myfamiglia$nome <- NULL
#Cambio ordine colonne
myfamiglia <- myfamiglia[,c("nome_ID", "1")]
listfile <- list.files("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/T1/T1/")
#Ciclo for che percorre tutti i file caricati di T1
for (ciclo in 2:length(listfile)) {
  mydata <- fread(paste0("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/T1/T1/", listfile[ciclo]), data.table = F)
  #Cambio nome colonne
  names(mydata) <- c("percentuale", ciclo, "frammenti", "codice", "taxon_ID", "nome")
  #Seleziono le righe che hanno come codice F
  myciclo <- subset(mydata, codice == "F")
  #Elimino alcune colonne
  myciclo$percentuale <- myciclo$frammenti <- myciclo$codice <- NULL
  #Creo nuova colonna unendo nome e taxon ID
  myciclo$nome_ID <- paste(myciclo$nome, myciclo$taxon_ID, sep = "_")
  myciclo$taxon_ID <- myciclo$nome <- NULL
  #Cambio gli spazi con _
  myciclo$nome_ID <- gsub(" ", "_", myciclo$nome_ID)
  myciclo <- myciclo[,c("nome_ID", ciclo)]
  myfamiglia <- merge(myfamiglia, myciclo, by = "nome_ID", all=T)
}
#Sostituisco NA con 0
myfamiglia [is.na(myfamiglia)]=0
write.xlsx(myfamiglia, paste0 ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/", "conteggio_famiglie1.xlsx"), row.names=F)

#PHYLUM
#Carico i dati
mydata <- fread("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/T1/T1/107805_ID2646_1-1-A01-ITS_S260_L001_report.txt", data.table = F)
head(mydata)
#Cambio nome colonne
names(mydata) <- c("percentuale", "1", "frammenti", "codice", "taxon_ID", "nome")
#Seleziono le righe che hanno come codice P
myphylum <- subset(mydata, codice == "P")
#Elimino alcune colonne
myphylum$percentuale <- myphylum$frammenti <- myphylum$codice <- NULL
#Creo nuova colonna unendo nome e taxon ID
myphylum$nome_ID <- paste(myphylum$nome, myphylum$taxon_ID, sep = "_")
myphylum$taxon_ID <- myphylum$nome <- NULL
#Cambio ordine colonne
myphylum <- myphylum[,c("nome_ID", "1")]
listfile <- list.files("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/T1/T1/")
#Ciclo for che percorre tutti i file caricati di T1
for (ciclo in 2:length(listfile)) {
  mydata <- fread(paste0("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/T1/T1/", listfile[ciclo]), data.table = F)
  #Cambio nome colonne
  names(mydata) <- c("percentuale", ciclo, "frammenti", "codice", "taxon_ID", "nome")
  #Seleziono le righe che hanno come codice P
  myciclo <- subset(mydata, codice == "P")
  #Elimino alcune colonne
  myciclo$percentuale <- myciclo$frammenti <- myciclo$codice <- NULL
  #Creo nuova colonna unendo nome e taxon ID
  myciclo$nome_ID <- paste(myciclo$nome, myciclo$taxon_ID, sep = "_")
  myciclo$taxon_ID <- myciclo$nome <- NULL
  #Cambio gli spazi con _
  myciclo$nome_ID <- gsub(" ", "_", myciclo$nome_ID)
  myciclo <- myciclo[,c("nome_ID", ciclo)]
  myphylum <- merge(myphylum, myciclo, by = "nome_ID", all=T)
}
#Sostituisco NA con 0
myphylum [is.na(myphylum)]=0
write.xlsx(myphylum, paste0 ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/", "conteggio_phyla1.xlsx"), row.names=F)
