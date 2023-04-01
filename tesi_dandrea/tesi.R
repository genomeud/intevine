#Carico le librerie
library(multcompView)

grapevine <- read.csv2("grapevine.csv", encoding = "latin1", header = T, sep = ";")


# Rinomino le righe con le etichette della prima colonna e poi la elimino
row.names(grapevine) <- grapevine[,1]
grapevine[,1] <- NULL

# Trasformo Soil, ster_soil e Ster_root in factor
grapevine$Suolo <- as.factor(grapevine$Suolo)
grapevine$Suolo_autoclav <- as.factor(grapevine$Suolo_autoclav)
grapevine$Termoterapia <- as.factor(grapevine$Termoterapia)

# Mediana di ciascun parametro
# Suolo
med_soil <- round(aggregate(grapevine[,5:10], by = list(grapevine$Suolo), FUN = "median")[,2:7],2)
rownames(med_soil) <- c("Sabbia", "Stallatico", "Torba")
# median(grapevine[which(grapevine$Suolo == "Sabbia"),]$Conucibilità.µS.cm)

# Ster soil
med_stersoil <- round(aggregate(grapevine[,5:10], by = list(grapevine$Suolo_autoclav), FUN = "median")[,2:7],2)
rownames(med_stersoil) <- c("S_non autocl", "S_autoclavato")
# median(grapevine[which(grapevine$Suolo_autoclav == "Si"),]$Conucibilità.µS.cm)

#Ster root
med_sterroot <- round(aggregate(grapevine[,5:10], by = list(grapevine$Termoterapia), FUN = "median")[,2:7],2)
rownames(med_sterroot) <- c("No termoter", "Termoterapia")


# STRIPCHART E TEST DI WILCOXON

#
tri.to.squ<-function(x)
{
  rn<-row.names(x)
  cn<-colnames(x)
  an<-unique(c(cn,rn))
  myval<-x[!is.na(x)]
  mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
  for(ext in 1:length(cn))
  {
    for(int in 1:length(rn))
    {
      if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
      mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
      mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
    }
    
  }
  return(mymat)
}

# PARAMETRI GRAFICI
# font.axis = 4 (font degli assi) grassetto corsivo
# font.lab = 3 (font delle etichette degli assi) corsivo
# (pch = 22, col = "blue", bg = "orange") bg riempimento, col = contorno


# 1)SUOLO  
par(mfrow = c(2,3))

col_Str <- c("goldenrod2","lightsalmon3","yellowgreen") # Sabbia, stallatico, torba
# cloripunti <- col_Str[grapevine$Soil]

library(multcompView)

# pH
w_pH <- pairwise.wilcox.test(grapevine[,5], grapevine[,2], p.adjust.method = "none", paired = F)
mymat_pH <- tri.to.squ(w_pH$p.value)
myletters <- multcompLetters(mymat_pH, compare = "<=", threshold = 0.05, Letters = letters)
stripchart(grapevine$pH ~ grapevine$Suolo, vertical = T, method = "jitter", jitter = 0.2, pch = 19 , cex = 2, cex.main = 1.2, 
           col = col_Str, las = 1, frame = F, font.axis = 3, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5, font.lab = 3, font.main = 4, 
           ylab = "pH", xlab = "Suolo", main = "pH\n(Suolo)")
med_S_pH <- aggregate(grapevine$pH, by = list(grapevine$Suolo), FUN = "median")
text(c(1,2,3),med_S_pH[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)
legend("topright", legend = levels(grapevine$Suolo), bty = "n",fill = col_Str, text.font = 3)

# Conductivity 
w_cond <- pairwise.wilcox.test(grapevine[,6], grapevine[,2], p.adjust.method = "none", paired = F)
mymat_c <- tri.to.squ(w_cond$p.value)
myletters <- multcompLetters(mymat_c, compare = "<=", threshold = 0.05, Letters = letters)

stripchart(grapevine$Conucibilità.µS.cm ~ grapevine$Suolo, vertical = T, method = "jitter", jitter = 0.2, pch = 19 , cex = 2, 
           cex.main = 1.2, col = col_Str,las = 1, ylab = "Conucibilità µS/cm", xlab = "Suolo", main = "Conducibilità\n(Suolo)", 
           font.axis = 3, font.lab = 3, font.main = 4, frame = F, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)

med_S_C <- aggregate(grapevine$Conucibilità.µS.cm, by = list(grapevine$Suolo), FUN = "median")
text(c(1,2,3),med_S_C[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# DOC
w_DOC <- pairwise.wilcox.test(grapevine[,7], grapevine[,2], p.adjust.method = "none", paired = F)
mymat_DOC <- tri.to.squ(w_DOC$p.value)
myletters <- multcompLetters(mymat_DOC, compare = "<=", threshold = 0.05, Letters = letters)
stripchart(grapevine$DOC.mg.g.suolo ~ grapevine$Suolo, vertical = T, method = "jitter",jitter = 0.2, pch = 19 , cex = 2, 
           cex.main = 1.2, col = col_Str, las = 1, frame = F, ylab = "DOC mg/g suolo", xlab = "Suolo", main = "DOC\n(Suolo)",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
med_S_DOC <- aggregate(grapevine$DOC.mg.g.suolo, by = list(grapevine$Suolo), FUN = "median")
text(c(1,2,3),med_S_DOC[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# DN
w_DN <- pairwise.wilcox.test(grapevine[,8], grapevine[,2], p.adjust.method = "none", paired = F)
mymat_DN <- tri.to.squ(w_DN$p.value)
myletters <- multcompLetters(mymat_DN, compare = "<=", threshold = 0.05, Letters = letters)
stripchart(grapevine$DN.mg.g.suolo ~ grapevine$Suolo, vertical = T, method = "jitter", jitter = 0.2, pch = 19 , cex = 2, cex.main = 1.2,
           col = col_Str, las = 1, frame = F, ylab = "DN mg/g suolo", xlab = "Suolo", main = "DN\n(Suolo)",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
med_S_DN <- aggregate(grapevine$DN.mg.g.suolo, by = list(grapevine$Suolo), FUN = "median")
text(c(1,2,3),med_S_DN[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# ATP
w_ATP <- pairwise.wilcox.test(grapevine[,9], grapevine[,2], p.adjust.method = "none", paired = F)
mymat_ATP <- tri.to.squ(w_ATP$p.value)
myletters <- multcompLetters(mymat_ATP, compare = "<=", threshold = 0.05, Letters = letters)
stripchart(grapevine$µgATP.g.suolo ~ grapevine$Suolo, vertical = T, method = "jitter", jitter = 0.2, pch = 19 , cex = 2, cex.main = 1.2, 
           col = col_Str, las = 1, frame = F, ylab = "µgATP/g suolo", xlab = "Suolo", main = "ATP\n(Suolo)",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
med_S_ATP <- aggregate(grapevine$µgATP.g.suolo, by = list(grapevine$Suolo), FUN = "median")
text(c(1,2,3),med_S_ATP[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# DNA
w_DNA <- pairwise.wilcox.test(grapevine[,10], grapevine[,2], p.adjust.method = "none", paired = F)
mymat_DNA <- tri.to.squ(w_DNA$p.value)
myletters <- multcompLetters(mymat_DNA, compare = "<=", threshold = 0.05, Letters = letters)
stripchart(grapevine$DNA.ng.µL ~ grapevine$Suolo, vertical = T, method = "jitter",jitter = 0.2, pch = 19 , cex = 2, cex.main = 1.2, 
           col = col_Str, las = 1, frame = F, ylab = "DNA ng/µL", xlab = "Suolo", main = "DNA\n(Suolo)", 
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
med_S_DNA <- aggregate(grapevine$DNA.ng.µL, by = list(grapevine$Suolo), FUN = "median")
text(c(1,2,3),med_S_DNA[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)


# 2) SUDDIVISIONE SUOLO AUTOCLAVATO 
par(mfrow = c(2,3))
colors <- c("orange1", "olivedrab3")

# pH
w_pH <- pairwise.wilcox.test(grapevine[,5], grapevine[,3], p.adjust.method = "none", paired = F)
mymat_pH <- tri.to.squ(w_pH$p.value)
myletters <- multcompLetters(mymat_pH, compare = "<=", threshold = 0.05, Letters = letters)

stripchart(grapevine$pH  ~ grapevine$Suolo_autoclav, vertical = T, method = "jitter", jitter = 0.2, pch = 19, col = colors, cex = 2,
           las = 1, frame = F, ylab = "pH", xlab = "Suolo autoclavato", main ="pH", 
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
med_pH <- aggregate(grapevine[,5], by = list(grapevine$Suolo_autoclav), FUN= "median")
text(c(1,2,3),med_pH[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# Conductivity
w_cond <- pairwise.wilcox.test(grapevine[,6], grapevine[,3], p.adjust.method = "none", paired = F)
mymat_c <- tri.to.squ(w_cond$p.value)
myletters <- multcompLetters(mymat_c, compare = "<=", threshold = 0.05, Letters = letters)

stripchart(grapevine$Conucibilità.µS.cm  ~ grapevine$Suolo_autoclav, vertical = T, method = "jitter", jitter = 0.2, pch = 19, col = colors, cex = 2,
           las = 1, frame = F, ylab = "Conducibilità µS/cm", xlab = "Suolo autoclavato", main ="Conducibilità",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
med_C <- aggregate(grapevine[,6], by = list(grapevine$Suolo_autoclav), FUN= "median")
text(c(1,2,3),med_C[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# DOC
w_DOC <- pairwise.wilcox.test(grapevine[,7], grapevine[,3], p.adjust.method = "none", paired = F)
mymat_DOC <- tri.to.squ(w_DOC$p.value)
myletters <- multcompLetters(mymat_DOC, compare = "<=", threshold = 0.05, Letters = letters)

stripchart(grapevine$DOC.mg.g.suolo  ~ grapevine$Suolo_autoclav, vertical = T, method = "jitter", jitter = 0.2, pch = 19, col = colors, cex = 2,
           las = 1, frame = F, ylab = "DOC mg/g suolo", xlab = "Suolo autoclavato", main ="DOC",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
med_DOC <- aggregate(grapevine[,7], by = list(grapevine$Suolo_autoclav), FUN= "median")
text(c(1,2,3),med_DOC[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# DN
w_DN <- pairwise.wilcox.test(grapevine[,8], grapevine[,3], p.adjust.method = "none", paired = F)
mymat_DN <- tri.to.squ(w_DN$p.value)
myletters <- multcompLetters(mymat_DN, compare = "<=", threshold = 0.05, Letters = letters)

stripchart(grapevine$DN.mg.g.suolo  ~ grapevine$Suolo_autoclav, vertical = T, method = "jitter", jitter = 0.2, pch = 19, col = colors, cex = 2,
           las = 1, frame = F, ylab = "DN mg/g suolo", xlab = "Suolo autoclavato", main ="DN",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
med_DN <- aggregate(grapevine[,8], by = list(grapevine$Suolo_autoclav), FUN= "median")
text(c(1,2,3),med_DN[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# ATP
w_ATP <- pairwise.wilcox.test(grapevine[,9], grapevine[,3], p.adjust.method = "none", paired = F)
mymat_ATP <- tri.to.squ(w_ATP$p.value)
myletters <- multcompLetters(mymat_ATP, compare = "<=", threshold = 0.05, Letters = letters)

stripchart(grapevine$µgATP.g.suolo  ~ grapevine$Suolo_autoclav, vertical = T, method = "jitter", jitter = 0.2, pch = 19, col = colors, cex = 2,
           las = 1, frame = F, ylab = "µgATP/g suolo", xlab = "Suolo autoclavato", main ="ATP",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
med_ATP <- aggregate(grapevine[,9], by = list(grapevine$Suolo_autoclav), FUN= "median")
text(c(1,2,3),med_ATP[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)
label <- c("S.non autocl.", "Autocl.")
legend("topright", legend = label, fill = colors, bty = "n")

# DNA
w_DNA <- pairwise.wilcox.test(grapevine[,10], grapevine[,3], p.adjust.method = "none", paired = F)
mymat_DNA <- tri.to.squ(w_DNA$p.value)
myletters <- multcompLetters(mymat_DNA, compare = "<=", threshold = 0.05, Letters = letters) 
stripchart(grapevine$DNA.ng.µL  ~ grapevine$Suolo_autoclav, vertical = T, method = "jitter", jitter = 0.2, pch = 19, col = colors, cex =2,
           las = 1, frame = F, ylab = "DNA ng/µL", xlab = "Suolo autoclavato", main ="DNA",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
med_DNA <- aggregate(grapevine[,10], by = list(grapevine$Suolo_autoclav), FUN= "median")
text(c(1,2,3),med_DNA[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# 3) SUDDIVISIONE RADICE STERILIZZATA
par(mfrow = c(2,3))

colors <- c("orange1", "olivedrab3")

# pH
w_pH <- pairwise.wilcox.test(grapevine[,5], grapevine[,4], p.adjust.method = "none", paired = F)
mymat_pH <- tri.to.squ(w_pH$p.value)
myletters <- multcompLetters(mymat_pH, compare = "<=", threshold = 0.05, Letters = letters)

stripchart(grapevine$pH  ~ grapevine$Termoterapia, vertical = T, method = "jitter", jitter = 0.2, pch = 19, col = colors, cex = 2,
           las = 1, frame = F, ylab = "pH", xlab = "Termoterapia radice", main ="pH",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
medd_pH <- aggregate(grapevine[,5], by = list(grapevine$Termoterapia), FUN= "median")
text(c(1,2,3),medd_pH[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# Conductivity
w_cond <- pairwise.wilcox.test(grapevine[,6], grapevine[,4], p.adjust.method = "none", paired = F)
mymat_c <- tri.to.squ(w_cond$p.value)
myletters <- multcompLetters(mymat_c, compare = "<=", threshold = 0.05, Letters = letters)

stripchart(grapevine$Conucibilità.µS.cm  ~ grapevine$Termoterapia, vertical = T, method = "jitter", jitter = 0.2, pch = 19, col = colors, cex = 2,
           las = 1, frame = F, ylab = "Conducibilità µS/cm", xlab = "Termoterapia radice", main ="Conducibilità",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
medd_C <- aggregate(grapevine[,6], by = list(grapevine$Termoterapia), FUN= "median")
text(c(1,2,3),medd_C[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# DOC
w_DOC <- pairwise.wilcox.test(grapevine[,7], grapevine[,4], p.adjust.method = "none", paired = F)
mymat_DOC <- tri.to.squ(w_DOC$p.value)
myletters <- multcompLetters(mymat_DOC, compare = "<=", threshold = 0.05, Letters = letters) 

stripchart(grapevine$DOC.mg.g.suolo  ~ grapevine$Termoterapia, vertical = T, method = "jitter", jitter = 0.2, pch = 19, col = colors, cex = 2,
           las = 1, frame = F, ylab = "DOC mg/g suolo", xlab = "Termoterapia radice", main ="DOC",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
medd_DOC <- aggregate(grapevine[,7], by = list(grapevine$Termoterapia), FUN= "median")
text(c(1,2,3),medd_DOC[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# DN
w_DN <- pairwise.wilcox.test(grapevine[,8], grapevine[,4], p.adjust.method = "none", paired = F)
mymat_DN <- tri.to.squ(w_DN$p.value)
myletters <- multcompLetters(mymat_DN, compare = "<=", threshold = 0.05, Letters = letters)

stripchart(grapevine$DN.mg.g.suolo  ~ grapevine$Termoterapia, vertical = T, method = "jitter",jitter = 0.2, pch = 19, col = colors, cex = 2,
           las = 1, frame = F, ylab = "DN mg/g suolo", xlab = "Termoterapia radice", main ="DN",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
medd_DN <- aggregate(grapevine[,8], by = list(grapevine$Termoterapia), FUN= "median")
text(c(1,2,3),medd_DN[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# ATP
w_ATP <- pairwise.wilcox.test(grapevine[,9], grapevine[,4], p.adjust.method = "none", paired = F)
mymat_ATP <- tri.to.squ(w_ATP$p.value)
myletters <- multcompLetters(mymat_ATP, compare = "<=", threshold = 0.05, Letters = letters)

stripchart(grapevine$µgATP.g.suolo  ~ grapevine$Termoterapia, vertical = T, method = "jitter", jitter = 0.2, pch = 19, col = colors, cex = 2,
           las = 1, frame = F, ylab = "µgATP/g suolo", xlab = "Termoterapia radice", main ="ATP",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
medd_ATP <- aggregate(grapevine[,9], by = list(grapevine$Termoterapia), FUN= "median")
text(c(1,2,3),medd_ATP[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)

# DNA
w_DNA <- pairwise.wilcox.test(grapevine[,10], grapevine[,4], p.adjust.method = "none", paired = F)
mymat_DNA <- tri.to.squ(w_DNA$p.value)
myletters <- multcompLetters(mymat_DNA, compare = "<=", threshold = 0.05, Letters = letters)

stripchart(grapevine$DNA.ng.µL  ~ grapevine$Termoterapia, vertical = T, method = "jitter", jitter = 0.2, pch = 19, col = colors, cex = 2,
           las = 1, frame = F, ylab = "DNA ng/µl", xlab = "Termoterapia radice", main ="DNA",
           font.axis = 3, font.lab = 3, font.main = 4, cex.axis = 1.2,cex.lab = 1.2,cex.main=1.5)
medd_DNA <- aggregate(grapevine[,10], by = list(grapevine$Termoterapia), FUN= "median")
text(c(1,2,3),medd_DNA[,2],c(myletters$Letters[1],myletters$Letters[2],myletters$Letters[3]), cex = 1.7)
label <- c("R.non\ntermotr.", "R.termotr.")
legend("topleft", legend = label, fill = colors, bty = "n")


# CORRELOGRAMMA
library(Hmisc)
library(corrplot)
subset_cor <- grapevine[, 7:10] # creo df con i soli valori numerici da DOC a DNA

# se applico rcorr direttamete su subset -->  'list' object cannot be coerced to type 'double'
sc_subset_cor <- scale(subset_cor)
cormat2 <- rcorr(sc_subset_cor, type = "spearman")

library(RColorBrewer)

#trace(corrplot, edit=TRUE) # add +0.25 to [sig.locs])+0.25 to row 446

corrplot(cormat2$r, p.mat = cormat2$P ,type = "lower", tl.col = "black", col=colorRampPalette(c("olivedrab1","lightsalmon2"))(100), 
         diag = F, insig = "label_sig", pch.cex = 1.5, addCoef.col = "black",tl.srt=45, cl.align.text ="l", cl.ratio = 0.2,
         number.font = 3, addgrid.col = colorRampPalette(c("olivedrab1","lightsalmon2"))(100))
cormat2$P

## CLUSTERING

# Standardizzazione
sc_grape <- scale(grapevine[,5:10])

# Matrice di distanza
grape_dist <- dist(sc_grape)

# Clustering gerarchico 
hc_grape <- hclust(grape_dist, method = "complete")

# Dendrogramma con colore delle etichette e dei rami diverso in base al tipo di suolo 
library(colorhcplot)

colorhcplot(hc_grape,grapevine$Suolo, colors = col_Str, main = "Cluster analisi\n(Tipologia di suolo)")
# Suolo ster
colorhcplot(hc_grape,grapevine$Suolo_autoclav, colors = c("orange1", "olivedrab3"), main = "Cluster analisi\n(Suolo autoclavato)")
# radice ster
colorhcplot(hc_grape,grapevine$Termoterapia, colors = c("orange1", "olivedrab3"), main = "Cluster analisi\n(Termoterapia della radice)")


## DNA reads

library(data.table)
tdat <- data.table()
fdat <- data.table()
mydata <- dir("." , pattern = "_report.txt") # crea un vettore con tutti gli elementi che nella cartella corrente terminano con "...txt"


for (k in 1:length(mydata)) {
  tdat <- fread(mydata[k], data.table = F)  # df del file contenuti in kraken
  tdat <- tdat[which(tdat$V4== "G"),c("V6","V1")]
  names(tdat) <- c("Genus", mydata[k])
  if (k == 1) {
    fdat <- tdat
  } else {
    fdat <- merge(fdat, tdat, "Genus", sort=F, all = T)
  }
}

# fdat è un df con 1006 righe(generi)*37 colonne (la prima colonna è il nome del genere)

# la tabella a disposizione è il risultato di kraken 2
# colonna 1 -> % di read assegnate ad un taxon all'interno di un campione
# colonna 4 -> iniziale del livello tassonomico, S = specie, G = genere
# colonna 5 -> taxon ID di NCBI (codice univoco utilizzabile sulla banca)
# colonna 6 -> nome colloquiale del taxon

row.names(fdat) <- fdat[,1] # la prima colonna -genere- diventa il nome delle righe

fdat[,1] <- NULL # elimino la prima colonna

# Porre i Nas uguale a zero
fdat[is.na(fdat)] <- 0 # inserisce 0 al posto dei valori NA. evito di perdere dei campioni in cui potrebbero  essere presenti dei generi

#aggiungo colonna con totale riga (somma delle percentuali di ciscun genere) NB 1: per le righe, 2 per colonne
fdat$total <- apply(fdat, 1, sum) 

# ordino il dataframe sulla base dei generi più abbondanti
fdat <-fdat[order(fdat$total, decreasing = T),]

# se total è 0.36 significa che ciascun campione ha almeno 0.01% di reads riferite ad un dato genere (oppure 0 e alcuni campioni >0.01%)
# fdat[which(fdat$total<=0.36),]
# Mantengo solo i generi che hanno valore magg/uguale a 0.36 nella colonna total
fdat <- fdat[fdat$total >= 0.36,]

# Abbondanza percentuale media dei primi 10 generi
fdat$mean_abundance <- fdat$total/36
fdat[1:10, 37:38]
fdat$mean_abundance <- NULL

# Trasposta
# elimino la colonna total 
fdat$total <- NULL 
tfdat <- t(fdat) # ottengo un oggetto tipo matrice. 36 righe (campioni)*430 colonne (generi)

tfdat <- as.data.frame(tfdat) # trasformo in df

# labels è un foglio excel creato per fare un merge con tfdat, per ottenere un dendrogramma con le etichette "LS,KT.." come nomi riga
# e colorazione diversa in base al tipo di suolo
labels <- read.csv("labels1.csv", header = F, sep = ";")

tfdat$V1 <- rownames(tfdat) # aggiungo la colonna V1 (ID) che servirà per eseguire il merge

unione <- merge(tfdat, labels, "V1", sort = F, all = T) # i generi vanno da colonna 2 a 431 
# v1 = prima colonna con ID2.., da 432 a 434 etichette (V2), colori (V3) e tipo terreno (V4)

grapevine$V2 <- rownames(grapevine)
gv <- grapevine[, c(3:4,11)] # df 3 colonne (etichette, ster_suolo e ster_rad)*36 righe
unione <- merge(unione, gv, "V2", sort = F, all = T)

row.names(unione) <-  unione$V2 # setto i giusti nomi di riga e poi elimino la colonna V2

grapevine$V2 <- NULL
unione[, 1:2] <- NULL # i generi vanno da colonna 1 a 430
unione$V4 <- as.factor(unione$V4)

# BARPLOT 10 GENERI PRINCIPALI

palette1 <- c("red","orange", "olivedrab1","green", "limegreen",
              "turquoise", "steelblue2", "slateblue1", "orchid" ,"violetred" )
barplot(t(unione[,1:10]), col = palette1, las = 2, cex.names = 0.75, ylim = c(0, 25), ylab = "Relative abundance (%)", font.lab = 4)
legend("topleft", legend = colnames(unione[,1:10]), fill = palette1, ncol = 3, bty = "n", text.font = 4, cex = 0.85) # elimina il bordo della legenda


# CLUSTER ANALISI
d_unione <- dist(unione[,1:430]) # i generi vanno da colonna 1 a 430
hc_unione <- hclust(d_unione)

# utilizzo funzione colorhcplot per poter asseganre i colori in base alla tipologia di suolo
library(colorhcplot)
unione$V4 <- as.factor(unione$V4)
colorhcplot(hc_unione, unione$V4, colors = col_Str, lab.cex = 1, 
            main = "Cluster analisi \n(Tipologia di suolo)")

# Trattamento suolo
colorhcplot(hc_unione, unione$Suolo_autoclav, colors = c("orange1", "olivedrab3"), main = "Cluster analisi \n(Suolo autoclavato)")

# Trattamento radice
colorhcplot(hc_unione, unione$Termoterapia, colors = c("orange1", "olivedrab3"), main = "Cluster analisi \n(Termoterapia della radice)")


# CORRELAZIONE TRA I GENERI E I PARAMETRI CHIMICO-FISICI DEL SUOLO
# definisco un unico df che nelle colonne presenta i 10 generi più abbondanti e i parametri chimici, le righe sono i diversi campioni
# al df creato è possibile applicare RCORR, funzione che fornisce come output sia la matrice con i coeff di corr, che la matrice
# con i rispettivi p-value

# creo un sottogruppo che contiene i 10 generi più abbondanti (meglio osservare i primi 10)
sub_set <- unione[,1:10]
sub_set$names <- rownames(sub_set) # aggiungo colonna names per eseguire il merge

grape <- grapevine[,5:10] 
grape$names <- rownames(grape) # aggiunta colonna names che serve per il merge

final <- merge(sub_set, grape, "names", sort = F, all = T) # final è un df con 36 righe(campioni) e 57 colonne(50 generi, 
# 6 parametri chimici, colonna names)
rownames(final)<- final[,1] # colonna names diventa nome di righe e poi elimino
final[,1]<- NULL # df 36 righe * 16 colonne, ultime 8 colonne sono paramteri chimici e tratt terreno e rad

# Standardizzo i valori del df final per renderli confrontabili 
sc_final <- scale(final)
library(Hmisc)
S_corr <- rcorr(sc_final, type = "spearman")

# ottengo la matrice con i coeff di corr
rho <- S_corr$r
rho <- rho[1:10,11:16] # oggetto di classe matrix, con 10 righe (Generi) e sei colonne (parametri chimici)

# ottengo la matrice con i p-value
pv <- S_corr$P 
pv <- pv[1:10,11:16]


## CORRELOGRAMMA
corrplot(rho, p.mat = pv, col = colorRampPalette(c("olivedrab3","orange1"))(100), tl.cex = 1, 
         tl.col = "black",insig = "label_sig", pch.cex = 1.5,addgrid.col=colorRampPalette(c("olivedrab3","orange1"))(100), 
         addCoef.col = "black",  number.font = 3, number.cex = 0.9)



## PHYLA
p_tdat <- data.table()
p_fdat <- data.table()
my_pdata <- dir("." , pattern = "_report.txt")


for (k in 1:length(my_pdata)) {
  p_tdat <- fread(my_pdata[k], data.table = F)  
  p_tdat <- p_tdat[which(p_tdat$V4== "P"),c("V6","V1")]
  names(p_tdat) <- c("Phylum", my_pdata[k])
  if (k == 1) {
    p_fdat <- p_tdat
  } else {
    p_fdat <- merge(p_fdat, p_tdat, "Phylum", sort=F, all = T)
  }
}


row.names(p_fdat) <- p_fdat[,1] # la prima colonna -Phylum- diventa il nome delle righe

p_fdat[,1] <- NULL # elimino la prima colonna

# Porre i Nas uguale a zero
p_fdat[is.na(p_fdat)] <- 0 # inserisce 0 al posto dei valori NA. evito di perdere dei campioni in cui potrebbero  essere presenti dei generi

#aggiungo colonna con totale riga (somma delle percentuali di ciscun phylum) 1 sta per righe, 2 per colonne
p_fdat$total <- apply(p_fdat, 1, sum) 

p_fdat <-p_fdat[order(p_fdat$total, decreasing = T),]

p_fdat$mean <- p_fdat$total/36
p_fdat[1:10, 37:38]
p_fdat$mean <- NULL
p_fdat$total <- NULL

t_pfdat <- t(p_fdat) # ottengo un oggetto tipo matrice. 36 righe (campioni)*40 colonne (phyla)

t_pfdat <- as.data.frame(t_pfdat) # trasformo in df

t_pfdat$V1 <- rownames(t_pfdat) # aggiungo la colonna V1 (ID) che servirà per eseguire il merge

uni_one <- merge(t_pfdat, labels, "V1", sort = F, all = T) # i phyla vanno da colonna 2 a 41
# v1 = prima colonna con ID2.., da 42 a 44 etichette (V2), colori (V3) e tipo terreno (V4)

rownames(uni_one) <- uni_one$V2
uni_one$V2 <- NULL
uni_one$V1 <- NULL

# BARPLOT 10 PHYLA PRINCIPALI

barplot(t(uni_one[1:36,1:10]), col = palette1, las = 2, cex.names = 0.75, ylim = c(0, 130),  ylab = "Relative abundance (%)", font.lab = 4)
legend("topleft", legend = colnames(uni_one[,1:10]), fill = palette1, ncol = 3, bty = "n", text.font = 4, cex = 0.7) # elimina il bordo della legenda

