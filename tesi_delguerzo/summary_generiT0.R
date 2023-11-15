#install.packages("dplyr")
library(dplyr)
library(openxlsx)
mygenere <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/percentuali_generi.xlsx")
is.numeric(mygenere[1,2])
#elimino il nome della prima colonna in modo che la colonna 1 dell'excel corrisponde alla prima riga di myfragment
rownames(mygenere) <- mygenere$nome_ID
mygenere$nome_ID <- NULL
#funzione dplyr che permette di fare un'operazione su ogni numero del dataframe convertendolo a numerico
mygenere <- mutate_all(mygenere, function(x) as.numeric(as.character(x)))
summary(mygenere)
