library(openxlsx)
library(data.table)
library(vegan)
library(ggplot2)
library(ggpubr)
mygenere <- read.xlsx ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/conteggio_generi1.xlsx")
rownames(mygenere) <- mygenere$nome_ID
mygenere$nome_ID <- NULL
mytra <- as.data.frame(t(mygenere))
#media della abbondanza totale del campione che userem come soglia per escludere i generi meno abbondanti
totalmedia <- summary (rowSums(mytra)) [4]
soglia <- totalmedia/1000
total <- colSums(mytra)
mytra <- mytra[,total > soglia]
#Indice di Shannon
shannon <- diversity(mytra, index = "shannon")
#Indice di Simpson
simpson <- diversity(mytra, index = "simpson")
#indice di Chao
chao <- estimateR(mytra)
write.table (shannon, ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/indice_Shannon1.txt"), quote = F)
write.table (simpson, ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/indice_Simpson1.txt"), quote = F)
write.table (chao, ("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/tabelle_excel/indice_Chao1.txt"), quote = F)

#SHANNON
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "Campionate")
#trasformo il vettore shannon in un dataframe
shannonframe <- data.frame (shannon)
#unisco mymetadata e shannonframe
shannonmeta <- merge(shannonframe, mymetadata, by.x = "row.names", by.y = "Sample_name")
#wilcoxon test
shannonwilcox <- pairwise.wilcox.test(shannonmeta$shannon, shannonmeta$Soil)
my_comparisons <- list(
  c("Peat", "Manure"),
  c("Peat", "Sand"),
  c("Manure", "Sand"))
ggplot(data = shannonmeta, aes(x=Soil, y=shannon, fill=Soil))+
  stat_compare_means(comparisons = my_comparisons)+
  geom_boxplot()+
  scale_fill_manual(values = c("black", "red", "yellow"))

#CHAO
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "Campionate")
chaoframe <- chao[row.names(chao)=="S.chao1", ,drop=F]
chaoframe <- as.data.frame(t(chaoframe))
names(chaoframe) <- c("chao")
#unisco mymetadata e chaoframe
chaometa <- merge(chaoframe, mymetadata, by.x = "row.names", by.y = "Sample_name")
#wilcoxon test
chaowilcox <- pairwise.wilcox.test(chaometa$chao, chaometa$Soil)
my_comparisons <- list(
  c("Peat", "Manure"),
  c("Peat", "Sand"),
  c("Manure", "Sand"))
ggplot(data = chaometa, aes(x=Soil, y=chao, fill=Soil))+
  stat_compare_means(comparisons = my_comparisons)+
  geom_boxplot()+
  scale_fill_manual(values = c("black", "red", "yellow"))

#SIMPSON
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "Campionate")
#trasformo il vettore simpson in un dataframe
simpsonframe <- data.frame (simpson)
#unisco mymetadata e simpsonframe
simpsonmeta <- merge(simpsonframe, mymetadata, by.x = "row.names", by.y = "Sample_name")
#wilcoxon test
simpsonwilcox <- pairwise.wilcox.test(simpsonmeta$simpson, simpsonmeta$Soil)
my_comparisons <- list(
  c("Peat", "Manure"),
  c("Peat", "Sand"),
  c("Manure", "Sand"))
ggplot(data = simpsonmeta, aes(x=Soil, y=simpson, fill=Soil))+
  stat_compare_means(comparisons = my_comparisons)+
  geom_boxplot()+
  scale_fill_manual(values = c("black", "red", "yellow"))

#RICHNESS
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "Campionate")
richnessframe <- chao[row.names(chao)=="S.obs", ,drop=F]
richnessframe <- as.data.frame(t(richnessframe))
names(richnessframe) <- c("richness")
#unisco mymetadata e richnessframe
richnessmeta <- merge(richnessframe, mymetadata, by.x = "row.names", by.y = "Sample_name")
#wilcoxon test
richnesswilcox <- pairwise.wilcox.test(richnessmeta$richness, richnessmeta$Soil)
my_comparisons <- list(
  c("Peat", "Manure"),
  c("Peat", "Sand"),
  c("Manure", "Sand"))
ggplot(data = richnessmeta, aes(x=Soil, y=richness, fill=Soil))+
  stat_compare_means(comparisons = my_comparisons)+
  geom_boxplot()+
  scale_fill_manual(values = c("black", "red", "yellow"))

#STERILIZATION
#SHANNON_STERILIZZAZIONE
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "Campionate")
#trasformo il vettore shannon in un dataframe
shannonframe1 <- data.frame (shannon)
#unisco mymetadata e shannonframe
shannonmeta1 <- merge(shannonframe1, mymetadata, by.x = "row.names", by.y = "Sample_name")
#wilcoxon test
shannonwilcox1 <- pairwise.wilcox.test(shannonmeta1$shannon, shannonmeta1$Ster_soil)
ggplot(data = shannonmeta1, aes(x=Ster_soil, y=shannon, fill=Ster_soil))+
  stat_compare_means(comparisons = list(c("Yes", "No")))+
  geom_boxplot()+
  scale_fill_manual(values = c("red", "green"))

#CHAO_STERILIZZAZIONE
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "Campionate")
chaoframe1 <- chao[row.names(chao)=="S.chao1", ,drop=F]
chaoframe1 <- as.data.frame(t(chaoframe1))
names(chaoframe1) <- c("chao")
#unisco mymetadata e chaoframe
chaometa1 <- merge(chaoframe1, mymetadata, by.x = "row.names", by.y = "Sample_name")
#wilcoxon test
chaowilcox1 <- pairwise.wilcox.test(chaometa1$chao, chaometa1$Ster_soil)
ggplot(data = chaometa1, aes(x=Ster_soil, y=chao, fill=Ster_soil))+
  stat_compare_means(comparisons = list(c("Yes", "No")))+
  geom_boxplot()+
  scale_fill_manual(values = c("red", "green"))

#SIMPSON_STERILIZZAZIONE
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "Campionate")
#trasformo il vettore simpson in un dataframe
simpsonframe1 <- data.frame (simpson)
#unisco mymetadata e simpsonframe
simpsonmeta1 <- merge(simpsonframe1, mymetadata, by.x = "row.names", by.y = "Sample_name")
#wilcoxon test
simpsonwilcox1 <- pairwise.wilcox.test(simpsonmeta1$simpson, simpsonmeta1$Ster_Soil)
ggplot(data = simpsonmeta1, aes(x=Ster_soil, y=simpson, fill=Ster_soil))+
  stat_compare_means(comparisons = list(c("Yes", "No")))+
  geom_boxplot()+
  scale_fill_manual(values = c("red", "green"))

#RICHNESS_STERILIZZAZIONE
mymetadata <- read.xlsx("C:/Users/denis/OneDrive/Desktop/DENISE/UNIVERSITY SAN/TESI/map_grapevine_final.xlsx", sheet = "Campionate")
richnessframe1 <- chao[row.names(chao)=="S.obs", ,drop=F]
richnessframe1 <- as.data.frame(t(richnessframe1))
names(richnessframe1) <- c("richness")
#unisco mymetadata e richnessframe
richnessmeta1 <- merge(richnessframe1, mymetadata, by.x = "row.names", by.y = "Sample_name")
#wilcoxon test
richnesswilcox1 <- pairwise.wilcox.test(richnessmeta1$richness, richnessmeta1$Ster_soil)
ggplot(data = richnessmeta, aes(x=Ster_soil, y=richness, fill=Ster_soil))+
  stat_compare_means(comparisons = list(c("Yes", "No")))+
  geom_boxplot()+
  scale_fill_manual(values = c("red", "green"))