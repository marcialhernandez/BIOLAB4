source("http://www.bioconductor.org/biocLite.R")
#biocLite("Biobase")
#biocLite("GEOquery")
#biocLite("genefilter")
library(Biobase)
library(GEOquery)
library(genefilter)

#options('download.file.method.GEOquery'='auto')

#datos <- getGEO(filename='GDS4974.soft.gz')

print (Meta(datos)$description[1])

print (Meta(datos)$feature_count)

print (Meta(datos)$sample_type)

#eset <- GDS2eSet(datos, do.log2=TRUE,GPL=NULL,AnnotGPL=FALSE)

# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 33297 features, 18 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM992693 GSM992694 ... GSM992710 (18 total)
# varLabels: sample stress description
# varMetadata: labelDescription
# featureData
# featureNames: 7896736 7896738 ... 7896730 (33297 total)
# fvarLabels: ID GB_LIST ... category (12 total)
# fvarMetadata: Column labelDescription
# experimentData: use 'experimentData(object)'
# pubMedIds: 23195993 


f2 <- function(x) {
  IQR(x) > 0.75
}

f3 <- function(x){
  is.na(x)
}

ffun2 <- filterfun(f3)

#sprdGns <- genefilter(exprs(eset), ffun2)
#print (sum(sprdGns))

#esetPreprocessed <-eset[sprdGns,]

#----------------------------
#Matriz de Expresion
#----------------------------

# mat_exp<-matrix(0, nrow = 33297, ncol = 18)
# for(a in 1:33297)
#   {
#   for(b in 1:18)
#     {
#     mat_exp[a,b]<-exprs(eset)[a,b]
#     }
# }

#print(mat_exp[1:10,1:5])

#----------------------------
#Probes
#----------------------------
# [1] "ID"              "GB_LIST"         "SPOT_ID"         "seqname"        
# [5] "RANGE_GB"        "RANGE_STRAND"    "RANGE_START"     "RANGE_STOP"     
# [9] "total_probes"    "gene_assignment" "mrna_assignment" "category" 


#$dimLabels
#[1] "featureNames"   "featureColumns"
featureDataEset<-featureData(eset)
print(methods(class=class(featureDataEset)))
print(methods(class=class(featureDataEset)[1]))
atributos<-attributes(featureDataEset)
#print (featureDataEset$"ID"[1])
#print (pData(featureDataEset)[2])
print(varMetadata(featureDataEset[1]))
print(featureDataEset[1]$GB_LIST)
print(c(varLabels(featureDataEset)[1],varLabels(featureDataEset)[4]))

#print(pData(featureDataEset)[1])
#print(sampleNames(featureDataEset)) ->NOMBRES (ID)

# atributosDatos<-featureData(eset)
# pID<-featureData(eset)$"ID"[1:100] #ok
# pGB_LIST<-featureData(eset)$"GB_LIST"[1:100] #ok nita split usando ','
# pSPOT_ID<-featureData(eset)$"SPOT_ID"[1:100]
# pseqname<-featureData(eset)$"seqname"[1:100]
# pRANGE_GB<-featureData(eset)$"RANGE_GB"[1:100]
# pRANGE_STRAND<-featureData(eset)$"RANGE_STRAND"[1:100]
# pRANGE_START<-featureData(eset)$"RANGE_START"[1:100]
# pRANGE_STOP<-featureData(eset)$"RANGE_STOP"[1:100]
# ptotal_probes<-featureData(eset)$"total_probes"[1:100]
# pgene_assignment<-featureData(eset)$"gene_assignment"[1:100]
#pmrna_assignment<-featureData(eset)$"mrna_assignment"
# pcategory<-featureData(eset)$"category"[1:100] 

# nombres<-c()
# print (featureData(eset)$"RANGE_GB"[1])
# print (typeof(featureData(eset)$"seqname"[1]))
# for (i in 1:100){
#   ultimoNombre<-tail(strsplit(featureData(eset)$"GB_LIST"[i], ","), n=1)
#   if (is.null(ultimoNombre) | is.na(ultimoNombre) | nchar(ultimoNombre)==0){
#     nombres<-c(nombres,featureData(eset)$"ID"[i])
#   }
#   
#   else{
#     nombres<-c(nombres,ultimoNombre)
#   }
# }
# 
# print (nombres)