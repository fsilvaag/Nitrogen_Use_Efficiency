library(GWASpoly)
library(data.table)
library(tidyverse)

setwd("G:/.shortcut-targets-by-id/1U9YXpkWxprY7XjilTWGYnZdMbnHhp06L/02.Articulos/sacarosa_paper/nitrogeno pruebas 2023/UEN_2023")


genofile <- "Genotypic_data_discovery_population_220_genotypes_ACGT_Format.csv.gz"
phenofile <- "Phenotypic_data_NUE_Low_High_N_Discovery.csv"

data <- read.GWASpoly(ploidy=10,pheno.file=phenofile,geno.file=genofile,
                      format="ACGT",n.traits=2,delim=",")

data.original <- set.K(data,LOCO=FALSE,n.core=1)

params <- set.params(geno.freq = 1 - 10/220,fixed=c("Q1", "Q2","Q3", "Q4","Q5"),
                     fixed.type=c("numeric", "numeric","numeric", "numeric","numeric"))

NUE <- GWASpoly(data.original,models=c("general","1-dom","2-dom","3-dom","4-dom","5-dom"),
                                     traits=c("UEN_High_N","UEN_Low_N"),
                                params=params,n.core=1)

NUE <- set.threshold(data = NUE, method = "FDR",level = 0.05, n.permute = 1000,n.core = 1)

qtl <- get.QTL(data = NUE, traits = NULL, models = NULL, bp.window = NULL)

library(openxlsx)
write.csv(qtl, "qtl_data_BLUPS_UEN_AEREO_fdr.csv")
save(NUE, file = "GWAS_Analysis_NUE_low_high_N_Discovery_populaiton.Rdata")