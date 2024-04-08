library(GWASpoly)
library(tidyverse)

source("//192.168.153.238/biodata10/rstudio/carolina/FSA/READ_GWAS_SNP_Remover.R")
source("//192.168.153.238/biodata10/rstudio/carolina/FSA/qqPlod_CSD.R")

load("GWAS_Analysis_NUE_low_high_N_Discovery_populaiton.Rdata")

genofile <- "Genotypic_data_discovery_population_220_genotypes_ACGT_Format.csv.gz"
phenofile <- "Phenotypic_data_NUE_Low_High_N_Discovery.csv"


NUE <- set.threshold(data = NUE, method = "FDR",level = 0.05, n.permute = 1000,n.core = 1)

qtl <- get.QTL(data = NUE, traits = NULL, models = NULL, bp.window = NULL)

# Filter using qqplot
# This step is realize for each genetic model

datos_UEN_high <- subset(qtl, Trait == "NUE_High_N" & Model == "general")

snp_remove <- datos_UEN_nivel_1$Marker[datos_UEN_nivel_1$Score >=5]


#########################################################################
############################# --REMOVE SNPs-- ###########################
#########################################################################

#Using the function for remover SNPs according to score
# this function is realized in different scores within genetic model, until eliminated 
#all SNPs no significatives

data_remove <- GWAS_CSD(ploidy=10,
                       pheno.file=phenofile,
                       geno.file=genofile,
                       format="ACGT",
                       n.traits=2,
                       delim=",", snp_remover = snp_remove)

data_remove<- set.K(data_SNP3,LOCO=FALSE,n.core=1)
params <- set.params(geno.freq = 1 - 10/220,
                     fixed=c("Q1", "Q2", "Q3", "Q4","Q5"),
                     fixed.type=c("numeric", "numeric","numeric","numeric", "numeric"))

data_for_score <- GWASpoly(data = data_remove,
                                           models=c("general"),
                                           traits=c("NUE_High_N"),
                                           params=params,
                                           n.core=1)


# QQplot generation

scores <- NUE@scores$NUE_High_N

scores = as.data.frame(NUE@scores$NUE_High_N) %>% rownames_to_column("SNP")
scores <- scores[,c("SNP","general")]

scores <- scores %>% separate("SNP", into = c("A","B","C"), remove = F, sep = "_")
scores$Chromosome <- ifelse(is.na(scores$C),scores$A, paste0(scores$A,"_",scores$B))
scores$Position <- as.numeric(ifelse(is.na(scores$C),scores$B, scores$C))
scores <- scores[,c(1, 6, 7, 5 )]

## Change for each of the score tests
sc4 = as.data.frame(data_for_score@scores$NUE_High_N) %>% rownames_to_column("SNP")
colnames(sc4)[2] <- "score evaluated"
scores <- left_join(scores, sc4)


qqCDS(scores, modelos =c("general", "score evaluated"), Q1 = 0.05, Q3 = 0.95)

# Repeat REMOVE SNPs until markers with false associations are eliminated in each model.


