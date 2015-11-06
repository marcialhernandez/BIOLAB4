source("http://www.bioconductor.org/biocLite.R")
#biocLite("Biobase")
#biocLite("GEOquery")
#biocLite("genefilter")
#install.packages( pkgs= "gplots")
print ("Precargando Librerias")
library(Biobase)
library(GEOquery)
library(genefilter)
library(gplots)

#options('download.file.method.GEOquery'='auto')

print ("Carga de datos Soft a variable global")
datos <- getGEO(filename='GDS2447.soft.gz')

#Impresion de caracteristicas de los datos por salida estandar
print (Meta(datos)$description[1])
print (Meta(datos)$feature_count)
print (Meta(datos)$sample_type)

cantidadGenes<-strtoi(Meta(datos)$feature_count)
cantidadMuestras<-strtoi(Meta(datos)$sample_count)

print ("Conversion de datos a objeto tipo ExpressionSet")
eset <- GDS2eSet(datos, do.log2=TRUE,GPL=NULL,AnnotGPL=FALSE)

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

mat_exp<-matrix(0, nrow = cantidadGenes, ncol = cantidadMuestras)
for(a in 1:cantidadGenes)
  {
  for(b in 1:cantidadMuestras)
    {
    mat_exp[a,b]<-exprs(eset)[a,b]
    }
}

featureDataEset<-featureData(eset)

print ("[Validacion] Imprime una porcion de la matriz de expresion")
print(mat_exp[1:10,1:5])

#Imprime los metodos disponibles para el objeto featureDataEset
#print(methods(class=class(featureDataEset)))

#Menciona los atributos disponibles de cada gen
#dataVarMetadata<-varMetadata(featureDataEset)

#Nombre de cada Gen. Siempre esta disponible
dataSampleNames<-sampleNames(featureDataEset) 

#Descripcion Gen si esta disponible
dataGene_Name<-featureData(eset)$"Gene_Name"[1:cantidadGenes]

#Nombre estandar si esta disponible
dataGene_Symbol<-featureData(eset)$"Gene_Symbol"[1:cantidadGenes]

#Descripcion Gen si esta disponible
dataGene_Product<-featureData(eset)$"Gene_Product"[1:cantidadGenes]

print ("Generando listado de genes y guardando informacion en genInfo.txt")
genIdOrEstandar<-paste(dataSampleNames,dataGene_Symbol, sep=":")
genInfo<-paste(genIdOrEstandar,dataGene_Name, sep=":")
write.table(genInfo, "genInfo.txt", sep="\t") 

#informacion disonible de cada muestra
#print (phenoData(eset))
#sample 
#disease.state 
#description
#Analisis de los posibles valores que pueden tomar
#phenoData(eset)$"disease.state"[1:15]
#disease.state:control|nicotine dependence
####### Los primeros 6 son nicotino dependientes   #######
####### Los demas (9) no son nicotino dependientes #######

#Informacion de los datos
#experimentData(eset)

muestras<-sampleNames(phenoData(eset))

muestras.clases<-0
muestras.cantidadControl<-0
muestras.cantidadDependientes<-0

for(i in 1:cantidadMuestras){
  if (phenoData(eset)$"disease.state"[i]=="control"){
    muestras.clases[i]=0
    muestras.cantidadControl<-muestras.cantidadControl+1
  }
  else{
    muestras.clases[i]=1
    muestras.cantidadDependientes<-muestras.cantidadDependientes+1
  }
}

print ("Verificando nivel de significancia")
pvalue<-0
c0 <- which(muestras.clases==0)
c1 <- which(muestras.clases==1)

mat_exp2<-matrix(0, nrow = cantidadGenes, ncol = cantidadMuestras)
for(i in 1:cantidadGenes){
  mat_exp2[i,]=c(mat_exp[i,c0],mat_exp[i,c1])
  pvalue[i]<-t.test(mat_exp[i,c0], mat_exp[i,c1])$p.value
}

muestras.clases2<-c(muestras.clases[c0],muestras.clases[c1])

#Ordenar matriz
mat_exp2<-mat_exp2[order(pvalue),]
probes<-genIdOrEstandar[order(pvalue)]
pvalue<-pvalue[order(pvalue)]

muestras.asignaColor <- function(valorClase) { 
  if (valorClase=="0"){
    return ("#FF0000") 
    }
  else {
    return ("#0000FF")
  }
}

muestras.mapeadoColor<-unlist(lapply(muestras.clases2, muestras.asignaColor))

print ("Creando EuclidianHeatMap15.pdf")
pdf("EuclidianHeatMap15.pdf", paper="a4", width=8, height=8)
muestras.EuclidianHeatMap15<-heatmap.2(mat_exp2[1:15,],col=greenred(50),xlab='Muestras',ylab='',
          labRow=probes[1:15],ColSideColors=muestras.mapeadoColor, scale="row", key=TRUE,
          symkey=FALSE, density.info="none", trace="none", cexRow=0.5,
          distfun = function(x) dist(x,method = 'euclidean'),labCol=muestras,
          hclustfun = function(x) hclust(x,method = 'complete'))
dev.off()

print ("Creando EuclidianHeatMap.pdf")
pdf("EuclidianHeatMap.pdf", paper="a4", width=8, height=8)
muestras.EuclidianClustergrama<-heatmap.2(mat_exp2[1:100,],col=greenred(50),xlab='Muestras',ylab='',
                                            labRow=probes[1:100],ColSideColors=muestras.mapeadoColor, scale="row", key=TRUE,
                                            symkey=FALSE, density.info="none", trace="none", cexRow=0.5,
                                            distfun = function(x) dist(x,method = 'euclidean'),labCol=muestras,
                                            hclustfun = function(x) hclust(x,method = 'complete'))
dev.off()

print ("Creando PearsonHeatMap15.pdf")
pdf("PearsonHeatMap15.pdf", paper="a4", width=8, height=8)
heatmap.2(t(mat_exp2[1:15,]), trace="none", density="none", 
          labRow=probes[1:15], scale="row", key=TRUE,
          symkey=FALSE, cexRow=0.5,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2))
dev.off()

print ("Creando PearsonHeatMap.pdf")
pdf("PearsonHeatMap.pdf", paper="a4", width=8, height=8)
heatmap.2(t(mat_exp2[1:100,]), trace="none", density="none", 
          labRow=probes[1:100], scale="row", key=TRUE,
          symkey=FALSE, cexRow=0.5,
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2))
dev.off()

print ("Creando ClusterizacionKMeans.pdf")
pdf("ClusterizacionKMeans.pdf", paper="a4", width=8, height=8)
mat_exp3 <- t(mat_exp2[1:1000,])[,nrow(mat_exp2[1:1000,]):1]
cl <- kmeans(mat_exp3, 2)
grafico<-plot(mat_exp3, col = cl$cluster, type='n')
text(mat_exp3, labels=muestras.clases2, col=cl$cluster)
points(cl$centers, col = 1:2, pch = 16, cex = 2)
title(main="Aplicacion de K-means")
dev.off()