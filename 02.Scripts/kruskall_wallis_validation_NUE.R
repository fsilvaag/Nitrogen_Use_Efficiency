library(nlcor)
rm(list = ls())
d = data.table::fread(".../Phenotypic_data_NUE_Low_High_N_Validation.csv", data.table = F)
d1 = data.table::fread(".../Markers for validation.csv", data.table = F)

d = d[,colnames(d) %in% c(d1$Marker,"model","nivel 1","nivel 2")]
d1 = subset(d1, d1$Marker %in% colnames(d))
colnames(d) <- gsub(" ","", colnames(d))

d$nivel2 = ifelse(is.na(d$nivel2), 0 , d$nivel2)

nm = unique(d1$Marker)

kwt1 = NULL
for (g in 1:length(nm)) {
  kwt = data.frame(Marker = character(), Nivel1 = numeric(), Nivel2 = numeric(), 
                   Nivel1ncl = numeric(), Nivel2ncl = numeric(),Modelo = character())
  gh = subset(d1, d1$Marker == nm[g])
  
  if(gh$modelo == "general") {
    asd = d[,c("model","nivel1","nivel2", nm[g])]
    colnames(asd) <- gsub(" ","", colnames(asd))
    aj = unique(asd[,nm[g]])
    if (length(aj) !=1 ) {
      x = model.matrix(~ 0 + as.factor(asd[,paste0(nm[g])]))
      colnames(x) = gsub(paste(c("as.factor","\\(","\\)","asd","\\[","\\,","paste0","\\(",
                                 "nm","\\[","g","\\]","\\)","\\]"," "), collapse = "|"),
                         "",colnames(x))
      colnames(x) <- paste0(nm[g],"_D",colnames(x))
      
      x1 = cbind(asd[,1:3],x)
      colnames(x1) <- gsub(" ","", colnames(x1))
      
      for (v in 4:(ncol(x1))) {
        
        kwt[v,2] <- kruskal.test(x1[,2], x1[,v])$p.value #Nivel 1
        kwt[v,3] <- kruskal.test(x1[,3], x1[,v])$p.value #Nivel 2
        kwt[v,1] = colnames(x1)[v]
        kwt[v,4] = nlcor(x1[,2], x1[,v], plt = F)
        kwt[v,5] = nlcor(x1[,3], x1[,v], plt = F)
        
        kwt[v,6] = gh$modelo
      }
      kwt1 = rbind(kwt1, subset(kwt, !is.na(kwt$Marker)))
      rm(kwt)
    }
    
  } else {
    asd = d[,c("model","nivel1","nivel2", nm[g])]
    colnames(asd) <- gsub(" ","",colnames(asd))
    aj = unique(asd[,nm[g]])
    if (length(aj) !=1 ) {
      
      kwt[g,1] <- nm[g]
      kwt[g,6] <- gh$modelo
      kwt[g,2] <- kruskal.test(asd[,"nivel1"], asd[,paste0(nm[g])])$p.value #Nivel 1
      kwt[g,3] <- kruskal.test(asd[,"nivel2"], asd[,paste0(nm[g])])$p.value #Nivel 2
      kwt[g,4] = nlcor(asd[,"nivel1"], asd[,nm[g]], plt = F)$adjusted.p.value #No utilizar
      kwt[g,5] = nlcor(asd[,"nivel2"], asd[,nm[g]], plt = F)$adjusted.p.value #No utilizar
      
      kwt1 = rbind(kwt1, subset(kwt, !is.na(kwt$Marker)))
      rm(kwt)
    }
    
  }
}


