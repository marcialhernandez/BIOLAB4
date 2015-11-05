source("http://www.bioconductor.org/biocLite.R")
#biocLite("Biobase")
#biocLite("GEOquery")
#biocLite("genefilter")
library(Biobase)
library(GEOquery)
library(genefilter)

#options('download.file.method.GEOquery'='auto')

#Carga de datos a variable global
#datos <- getGEO(filename='GDS2447.soft.gz')

#Impresion de caracteristicas de los datos por salida estandar
print (Meta(datos)$description[1])
print (Meta(datos)$feature_count)
print (Meta(datos)$sample_type)

#Conversion de datos a objeto tipo ExpressionSet
#eset <- GDS2eSet(datos, do.log2=TRUE,GPL=NULL,AnnotGPL=FALSE)

#Declaracion de funciones de filtrado

f2 <- function(x) {
  IQR(x) > 0.75
}

f3 <- function(x){
  is.na(x)
}

#Asignacion de funcion de filtrado
ffun2 <- filterfun(f3)

#Aplicacion de filtrado a objeto tipo 'ExpressionSet'
#sprdGns <- genefilter(exprs(eset), ffun2)
#print (sum(sprdGns))

#esetPreprocessed <-eset[sprdGns,]

#----------------------------
#Matriz de Expresion
#----------------------------

mat_exp<-matrix(0, nrow = 33202, ncol = 15)
for(a in 1:33202)
  {
  for(b in 1:15)
    {
    mat_exp[a,b]<-exprs(eset)[a,b]
    }
}

featureDataEset<-featureData(eset)

#[Validacion] Imprime una porcion de la matriz de expresion
print(mat_exp[1:10,1:5])

#Imprime los metodos disponibles para el objeto featureDataEset
#print(methods(class=class(featureDataEset)))

#Menciona los atributos disponibles de cada gen
#dataVarMetadata<-varMetadata(featureDataEset)

#Nombre de cada Gen. Siempre esta disponible
dataSampleNames<-sampleNames(featureDataEset) 

#Descripcion Gen si esta disponible
dataGene_Name<-featureData(eset)$"Gene_Name"[1:33202]

#Nombre estandar si esta disponible
dataGene_Symbol<-featureData(eset)$"Gene_Symbol"[1:33202]

#Descripcion Gen si esta disponible
dataGene_Product<-featureData(eset)$"Gene_Product"[1:33202]

#informacion disonible de cada muestra
#sample 
#disease.state 
#description
muestras<-sampleNames(phenoData(eset))[1]
