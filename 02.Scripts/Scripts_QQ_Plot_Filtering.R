
source("//192.168.153.238/biodata10/rstudio/carolina/FSA/READ_GWAS_SNP_Remover.R")
source("//192.168.153.238/biodata10/rstudio/carolina/FSA/qqPlod_CSD.R")

load("GWAS_Analysis_NUE_low_high_N_Discovery_populaiton.Rdata")

genofile <- "Genotypic_data_discovery_population_220_genotypes_ACGT_Format.csv.gz"
phenofile <- "Phenotypic_data_NUE_Low_High_N_Discovery.csv"



scores <- NUE@scores$NUE_High_N
scores <- scores %>% rownames_to_column("SNP") %>% 
  separate(SNP, into = c("A","B","r","a","c"), remove = F, sep = "_") %>% 
  mutate(Chromosome = ifelse(is.na(c), A, paste0(A,"_",B)),
         Position = as.numeric(ifelse(is.na(c), B, r)))

scores <- scores[,!colnames(scores) %in% c("A","B","c","r","a")]

scores <- scores[,c(1,13,14,2:12)]
qqCDS(scores, modelos = "general", Q1 = 0.05, Q3 = 0.95)

###########################
NUE <- set.threshold(data = NUE, method = "FDR",level = 0.05, n.permute = 1000,n.core = 1)

qtl <- get.QTL(data = NUE, traits = NULL, models = NULL, bp.window = NULL)

# Para saber cuales SNPS vas a remover
datos_UEN_nivel_1 <- subset(qtl, Trait == "NUE_High_N" & Model == "general")

snp5 <- datos_UEN_nivel_1$Marker[datos_UEN_nivel_1$Score >=5]

data_SNP3 <- GWAS_CSD(ploidy=10,
                       pheno.file=phenofile,
                       geno.file=genofile,
                       format="ACGT",
                       n.traits=2,
                       delim=",", snp_remover = snp5)

data_SNP2<- set.K(data_SNP3,LOCO=FALSE,n.core=1)
params <- set.params(geno.freq = 1 - 10/220,
                     fixed=c("Q1", "Q2", "Q3", "Q4","Q5"),
                     fixed.type=c("numeric", "numeric","numeric","numeric", "numeric"))

data_UEN_aereo_nivel_1_mayor_5 <- GWASpoly(data = data_SNP2,
                                           models=c("general"),
                                           traits=c("NUE_High_N"),
                                           params=params,
                                           n.core=1)


################
scores <- data_BLUPS_UEN_AEREO@scores$Blup_bk_fijo_UEN_areo_NIVEL_1

scores = as.data.frame(data_BLUPS_UEN_AEREO@scores$Blup_bk_fijo_UEN_areo_NIVEL_1) %>% rownames_to_column("SNP")
scores <- scores[,c("SNP","general")]

scores <- scores %>% separate("SNP", into = c("A","B","C"), remove = F, sep = "_")
scores$Chromosome <- ifelse(is.na(scores$C),scores$A, paste0(scores$A,"_",scores$B))
scores$Position <- as.numeric(ifelse(is.na(scores$C),scores$B, scores$C))
scores <- scores[,c(1, 6, 7, 5 )]

## Cambiar por cada una de las pruebas del score
sc4 = as.data.frame(data_UEN_aereo_nivel_1_mayor_5@scores$Blup_bk_fijo_UEN_areo_NIVEL_1) %>% rownames_to_column("SNP")
colnames(sc4)[2] <- "Gen_mayor_5"
scores <- left_join(scores, sc4)

sc5 = as.data.frame(data_UEN_aereo_nivel_1_mayor_626@scores$Blup_bk_fijo_UEN_areo_NIVEL_1) %>% rownames_to_column("SNP")
colnames(sc5)[2] <- "Gen_mayor_626"
scores <- left_join(scores, sc5)

sc6 = as.data.frame(data_UEN_aereo_nivel_1_mayor_55@scores$Blup_bk_fijo_UEN_areo_NIVEL_1) %>% rownames_to_column("SNP")
colnames(sc6)[2] <- "Gen_mayor_55"
scores <- left_join(scores, sc6)

qqCDS(scores, modelos =c("general", "Gen_mayor_626",  "Gen_mayor_55"), Q1 = 0.05, Q3 = 0.95)

########################### modelo 1-DOM-ALT

scores <- data_BLUPS_UEN_AEREO@scores$Blup_bk_fijo_UEN_areo_NIVEL_1
scores <- scores %>% rownames_to_column("SNP")
scores$Marker <- scores$SNP
scores <- scores %>% separate(Marker, into = c("A","B","c"))
scores$Chromosome <- ifelse(is.na(scores$c), scores$A, paste0(scores$A,"_",scores$B))
scores$Position <- ifelse(is.na(scores$c), scores$B, scores$c)
scores <- scores[,!colnames(scores) %in% c("A","B","c")]
#scores <- scores[,c(1,13,2:12)]
scores <- scores[,c(1,13,14,2:12)]
qqCDS(scores, modelos = "1-dom-alt", Q1 = 0.05, Q3 = 0.95)

###

datos_UEN_nivel_1 <- subset(qtl_data_BLUPS_UEN_AEREO_fdr, qtl_data_BLUPS_UEN_AEREO_fdr$Trait == "Blup_bk_fijo_UEN_areo_NIVEL_1" & qtl_data_BLUPS_UEN_AEREO_fdr$Model == "1-dom-alt")

snp2 <- datos_UEN_nivel_1$Marker[datos_UEN_nivel_1$Score >=4.8]

data_SNP3 <- GWAS_CSD (ploidy=10,
                       pheno.file=phenofile,
                       geno.file=genofile,
                       format="ACGT",
                       n.traits=2,
                       delim=",", snp_remover = snp2)

data_SNP2<- set.K(data_SNP3,LOCO=FALSE,n.core=1)
params <- set.params(geno.freq = 1 - 10/220,
                     fixed=c("Q1", "Q2", "Q3", "Q4"),
                     fixed.type=c("numeric", "numeric","numeric", "numeric"))

data_UEN_aereo_nivel_1_1domalt_48 <- GWASpoly(data = data_SNP2,
                                              models=c("1-dom"),
                                              traits=c("Blup_bk_fijo_UEN_areo_NIVEL_1"),
                                              params=params,
                                              n.core=1)


#setwd("G:/.shortcut-targets-by-id/1U9YXpkWxprY7XjilTWGYnZdMbnHhp06L/02.Articulos/sacarosa_paper/nitrogeno pruebas 2023/UEN_2023")

save(data_UEN_aereo_nivel_1_1domalt_48, file = "GWASPOLY_data_UEN_aereo_nivel_1_1domalt_48.RData")
load("GWASPOLY_data_UEN_aereo_nivel_1_1domalt_48.RData")
load("GWASPOLY_data_UEN_aereo_nivel_1_1domalt_5.RData")

scores = as.data.frame(data_BLUPS_UEN_AEREO@scores$Blup_bk_fijo_UEN_areo_NIVEL_1) %>% rownames_to_column("SNP")
scores <- scores[,c("SNP","1-dom-alt")]

scores <- scores %>% separate("SNP", into = c("A","B","C"), remove = F, sep = "_")
scores$Chromosome <- ifelse(is.na(scores$C),scores$A, paste0(scores$A,"_",scores$B))
scores$Position <- as.numeric(ifelse(is.na(scores$C),scores$B, scores$C))
scores <- scores[,c(1, 6, 7, 5 )]

## Cambiar por cada una de las pruebas del score
sc4 = as.data.frame(data_UEN_aereo_nivel_1_1domalt_5@scores$Blup_bk_fijo_UEN_areo_NIVEL_1) %>% rownames_to_column("SNP")
colnames(sc4)[2] <- "Gen_mayor_5"
scores <- left_join(scores, sc4)

sc5 = as.data.frame(data_UEN_aereo_nivel_1_1domalt_48@scores$Blup_bk_fijo_UEN_areo_NIVEL_1) %>% rownames_to_column("SNP")
colnames(sc5)[2] <- "Gen_mayor48"
scores <- left_join(scores, sc5)

qqCDS(scores, modelos =c("1-dom-alt", "Gen_mayor_5", "Gen_mayor48"), Q1 = 0.05, Q3 = 0.95)


#######################################################################################
############################      nivel 2          ####################################
#######################################################################################

scores <- data_BLUPS_UEN_AEREO@scores$Blup_bk_fijo_UEN_areo_N_2
scores <- scores %>% rownames_to_column("SNP")
scores$Marker <- scores$SNP
scores <- scores %>% separate(Marker, into = c("A","B","c"))
scores$Chromosome <- ifelse(is.na(scores$c), scores$A, paste0(scores$A,"_",scores$B))
scores$Position <- ifelse(is.na(scores$c), scores$B, scores$c)
scores <- scores[,!colnames(scores) %in% c("A","B","c")]
#scores <- scores[,c(1,13,2:12)]
scores <- scores[,c(1,13,14,2:12)]
qqCDS(scores, modelos = NULL, Q1 = 0.05, Q3 = 0.95)


datos_UEN_nivel_2 <- subset(qtl_data_BLUPS_UEN_AEREO_fdr, qtl_data_BLUPS_UEN_AEREO_fdr$Trait == "Blup_bk_fijo_UEN_areo_N_2" & qtl_data_BLUPS_UEN_AEREO_fdr$Model == "general")

snp5 <- datos_UEN_nivel_2$Marker[datos_UEN_nivel_2$Score >=5.5]

data_SNP3 <- GWAS_CSD (ploidy=10,
                       pheno.file=phenofile,
                       geno.file=genofile,
                       format="ACGT",
                       n.traits=2,
                       delim=",", snp_remover = snp5)

data_SNP2<- set.K(data_SNP3,LOCO=FALSE,n.core=1)
params <- set.params(geno.freq = 1 - 10/220,
                     fixed=c("Q1", "Q2", "Q3", "Q4"),
                     fixed.type=c("numeric", "numeric","numeric", "numeric"))

data_UEN_aereo_nivel_2_mayor_55 <- GWASpoly(data = data_SNP2,
                                            models=c("general"),
                                            traits=c("Blup_bk_fijo_UEN_areo_N_2"),
                                            params=params,
                                            n.core=1)


#setwd("G:/.shortcut-targets-by-id/1U9YXpkWxprY7XjilTWGYnZdMbnHhp06L/02.Articulos/sacarosa_paper/nitrogeno pruebas 2023/UEN_2023")

save(data_UEN_aereo_nivel_2_mayor_55, file = "GWASPOLY_data_UEN_aereo_nivel_2_mayor_55.RData")



scores = as.data.frame(data_BLUPS_UEN_AEREO@scores$Blup_bk_fijo_UEN_areo_N_2) %>% rownames_to_column("SNP")
scores <- scores[,c("SNP","general")]

scores <- scores %>% separate("SNP", into = c("A","B","C"), remove = F, sep = "_")
scores$Chromosome <- ifelse(is.na(scores$C),scores$A, paste0(scores$A,"_",scores$B))
scores$Position <- as.numeric(ifelse(is.na(scores$C),scores$B, scores$C))
scores <- scores[,c(1, 6, 7, 5 )]

## Cambiar por cada una de las pruebas del score
sc4 = as.data.frame(data_UEN_aereo_nivel_2_mayor_55@scores$Blup_bk_fijo_UEN_areo_N_2) %>% rownames_to_column("SNP")
colnames(sc4)[2] <- "Gen_mayor_55"
scores <- left_join(scores, sc4)

qqCDS(scores, modelos =c("general", "Gen_mayor_55"), Q1 = 0.05, Q3 = 0.95)

##############################################################################################
##############################################################################################


##################################################################
##################################################################
##################################################################
data_BLUPS_UEN_AEREO_fdr<- set.threshold(data_BLUPS_UEN_AEREO, "FDR")

qtl_data_BLUPS_UEN_AEREO_fdr_nivel_2<- get.QTL(data_BLUPS_UEN_AEREO_fdr, traits ="Blup_bk_fijo_UEN_areo_N_2")
qtl_data_BLUPS_UEN_AEREO_fdr_nivel_1<- get.QTL(data_BLUPS_UEN_AEREO_fdr, traits ="Blup_bk_fijo_UEN_areo_NIVEL_1")

r2_Nivel_2_UEN_aereo <- fit.QTL_CSD(data=data_BLUPS_UEN_AEREO_fdr,trait="Blup_bk_fijo_UEN_areo_N_2",
                                    qtl=qtl_data_BLUPS_UEN_AEREO_fdr_nivel_1[,c("Marker","Model")],
                                    fixed=data.frame(Effect=c("Q1", "Q2","Q3", "Q4"),
                                                     Type=c("numeric", "numeric","numeric","numeric")))

library(openxlsx)
write.csv(r2_Nivel_2_UEN_aereo, "r2_UEN_aereo_2023_log_fit.QTL_CSD_NIVEL_1.csv", sep=",")              
