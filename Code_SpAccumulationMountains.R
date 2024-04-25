######################################################################################################


#### El MARS debería tener un componente espacial de acuerdo a como lo venia estudiando. Leer cómo hacerlo

#Clean workspace
rm(list=ls())


#####################Load required packages
library(rgdal)
library(raster)
library(sp)
library(sf)
library(spatialEco)
library(phytools)
library(tibble)
library(picante)
library(BAT)
library(geiger)
library(FD)
library(matrixStats)
library(dplyr)
library(vegan)
library(ggplot2)
library(gridExtra)
library(lavaan)
library(quantreg)
require(mgcv) #GAM models
require(MuMIn) #for AICc
library(ggplot2)
library(sesem)
library(fields)

altura=raster("~/OneDrive - City University of New York/Documentos/Estadistica R SIG/Capas SIG/Rasters/Altitude_World.tif")

altura=crop(altura, c(-120, -20, -23, 23))
altura=aggregate(altura, 100, "mean")

America=st_read("~/OneDrive - City University of New York/Documentos/Estadistica R SIG/Capas SIG/Shapes/Political_World_Borders/America_Border", 'America_Border')
America=America[1]
America=as_Spatial(America)

altura2=mask(altura,  America)

Montanas=st_read('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/Documentos/IAvH Humboldt/Montañas/Códigos y Plots/Mapa/Shapes reducidos', "union")
Montanas=Montanas[1]
Montanas=as_Spatial(Montanas)

altura3=mask(altura2,  Montanas)

####Code from Daniel Rabosky Dan Rabosky
#https://github.com/macroevolution/fisse/blob/master/run_fisse/traitDependent_functions.R

DR_statistic <- function(x, return.mean = FALSE){
  
  rootnode <- length(x$tip.label) + 1
  
  sprates <- numeric(length(x$tip.label))
  
  for (i in 1:length(sprates)){
    node <- i
    index <- 1
    qx <- 0
    while (node != rootnode){
      el <- x$edge.length[x$edge[,2] == node]
      node <- x$edge[,1][x$edge[,2] == node]
      
      qx <- qx + el* (1 / 2^(index-1))
      
      index <- index + 1
      
    }
    
    sprates[i] <- 1/qx
    
  }
  
  
  if (return.mean) {
    
    return(mean(sprates)) 
    
  }else{
    
    names(sprates) <- x$tip.label
    return(sprates)
    
  }
  
  
}


##########################################################
##########################################################

setwd("~/OneDrive - City University of New York/Documentos/Estadistica R SIG/Capas SIG/Shapes/Mapas BirdLife 2020/Shapefiles BirdLife 2020")

tabla=read.delim("~/OneDrive - City University of New York/CUNY PhD/Thesis/Cap 1 Macro Montanas/Tablas, filogenias y bases de datos/Final_Birdfile_V5_2021.txt", header=T)

Na_Min=which(!is.na(tabla$Max) & is.na(tabla$Min))
tabla[Na_Min,'Min']<-0
tabla[Na_Min,'MidPoint']<-(tabla$Min[Na_Min]+tabla$Max[Na_Min])/2

Na_Min_Average=which(!is.na(tabla$Max_Average) & is.na(tabla$Min_Average))
tabla[Na_Min_Average,'Min_Average']<-0
tabla[Na_Min,'MidPoint_Average']<-(tabla$Min_Average[Na_Min_Average]+tabla$Max_Average[Na_Min_Average])/2

tabla[which(!is.na(tabla$Min) & is.na(tabla$Max)),c("Max", "Max_Average")]<-5500


Grupos=c("McGuire_Names_Trochilidae", "Harvey_etal_T400F_AOS_HowardMoore", "Barker_Names_Emberizoidea") #Poner las columnas de familia que hagan falta

Nombre_Grupo=c("Trochilidae", "Suboscines", "Oscines")

tree_Trochilidae <- read.tree('~/OneDrive - City University of New York/CUNY PhD/Thesis/Cap 1 Macro Montanas/Tablas, filogenias y bases de datos/Phylogenies/Aves/Trochilidae_PhyHenry.tre')
is.ultrametric(tree_Trochilidae)


tree_Suboscines <- read.tree("~/OneDrive - City University of New York/CUNY PhD/Thesis/Cap 1 Macro Montanas/Tablas, filogenias y bases de datos/Phylogenies/Aves/mgharvey-tyranni-f73aa7f/T400F_AOS_HowardMoore.tre")
is.ultrametric(tree_Suboscines)
tree_Suboscines <- force.ultrametric(tree_Suboscines)

tree_Oscines <- readNexus("~/OneDrive - City University of New York/CUNY PhD/Thesis/Cap 1 Macro Montanas/Tablas, filogenias y bases de datos/Phylogenies/Aves/Emberizoidea_Barker_etal_2015/Emberizoid_trees/sptree_MCC.tre")
is.ultrametric(tree_Oscines)

lista_arboles=list(tree_Trochilidae, tree_Suboscines, tree_Oscines)


#ClaDS data

ClaDS_Trochilidae=read.table('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/CUNY PhD/Thesis/Cap 1 Macro Montanas/ClaDS/ClaDS_Julia_Results/Humms_ClaDS/Humms_ClaDS_Pertips_Values.txt', header=T)
ClaDS_Oscines=read.table('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/CUNY PhD/Thesis/Cap 1 Macro Montanas/ClaDS/ClaDS_Julia_Results/Emberizoidea_ClaDS/Embe_ClaDS_Pertips_Values.txt', header=T)
ClaDS_Suboscines=read.table('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/CUNY PhD/Thesis/Cap 1 Macro Montanas/ClaDS/ClaDS_Julia_Results/Tyranni_ClaDS/Tyranni_ClaDS_Pertips_Values.txt', header=T)

list_ClaDS=list(ClaDS_Trochilidae, ClaDS_Suboscines, ClaDS_Oscines)

Finales=vector(mode = "list", length =3)
Mountains=vector(mode = "list", length =3)
Comunidades=vector(mode = "list", length =3)
Mapas=vector(mode = "list", length =3)


for(i in 1:3)
{
  setwd('~/OneDrive - City University of New York/Documentos/Estadistica R SIG/Capas SIG/Shapes/Mapas BirdLife 2020/Shapefiles BirdLife 2020')
  
  taxa=tabla[which(!is.na(tabla[,which(Grupos[i]==colnames(tabla))])),]
  
  taxa=taxa[which(!is.na(taxa$Birdlife_shapes_2021V1)),]
  
  community=as.data.frame(altura)
  
  for (j in 1:dim(taxa)[1])
    
  {
    print(j)
    
    X=taxa[j, 'Birdlife_shapes_2021V1']
    
    mapa=st_read("~/OneDrive - City University of New York/Documentos/Estadistica R SIG/Capas SIG/Shapes/Mapas BirdLife 2020/Shapefiles BirdLife 2020", strsplit(as.character(X), split="[.]")[[1]][1])
    
    r <- mapa[mapa$SEASONA == 1 | mapa$SEASONA == 2,]
    
    r=as_Spatial(r)
    
    mapa=rasterize(r,altura,field=1)
    
    mapa[is.na(mapa)]<-0
    
    #map <- raster(mapas_Sp[i])
    
    
    
    points=rasterToPoints(mapa,spatial=F)
    points=as.data.frame(points)
    points$layer=NULL
    coordinates(points)=c("x","y")
    alturas=extract(altura,points)
    altMax=taxa[which(X==taxa[,'Birdlife_shapes_2021V1']),'Max']
    altMin=taxa[which(X==taxa[,'Birdlife_shapes_2021V1']),'Min']
    mapa[which(alturas<altMin)]<-0
    mapa[which(alturas>altMax)]<-0
    
    dataRasterMap <- as.data.frame(mapa)
    colnames(dataRasterMap)<-taxa[j,which(Grupos[i]==colnames(taxa))]
    #dataRasterMap <- dataRasterMap[-comm2,]
    
    community <- cbind(community,dataRasterMap)
    
  }
  
  Elev=as.data.frame(altura)
  community[,1]<-NULL
  ColSUMA=colSums(community)
  community=community[,which(ColSUMA>0)]  ##################################### Pensa de que forma afecta las estimaciones de beta o divergencia morphologica
  community_2=community
  #community[community == 0] <- NA 
  community <- community %>% mutate_all(~na_if(., 0))

  ######
  #Species richness
  ######	
  
  Riq=rowSums(community[,1:dim(community)[2]], na.rm=T)
  Riqueza=altura
  Riqueza[]<-Riq
  Riqueza[which(is.na(Elev[,1]))] <- NA
  
  
  ######
  #Diversificiation rates
  ######	
  
  #Equal split
  DR=as.data.frame(DR_statistic(lista_arboles[[i]]))
  colnames(DR)=c("DR")
  DR2=as.data.frame(log(DR$DR))
  rownames(DR2)=rownames(DR)
  colnames(DR2)=c("DR")
  DR=DR2
  
  #ClaDS data
  ClaDS=list_ClaDS[[i]]
  rownames(ClaDS)=ClaDS$tiplabel

  Div.community_ES = community
  Div.community_ClaDS = community
  
  for (k in 1:dim(Div.community_ES)[2])
    
  {
    print(k)
    
    #ES
    sp.match = match(colnames(Div.community_ES[k]),rownames(DR))
    Div.community_ES[which(Div.community_ES[,k]>0),k] <- DR[sp.match,1]
    
    #ClaDS
    sp.match = match(colnames(Div.community_ClaDS[k]),rownames(ClaDS))
    Div.community_ClaDS[which(Div.community_ClaDS[,k]>0),k] <- ClaDS[sp.match,2]
  }
  
  Div.community_ES=rowMeans(Div.community_ES, na.rm=T)
  Div.community_ClaDS=rowMeans(Div.community_ClaDS, na.rm=T)
  
  Diversificacion_ES=altura
  Diversificacion_ClaDS=altura
  
  Diversificacion_ES[]<-Div.community_ES
  Diversificacion_ClaDS[]<-Div.community_ClaDS
  
  Diversificacion_ES[which(is.na(Elev[,1]))] <- NA
  Diversificacion_ClaDS[which(is.na(Elev[,1]))] <- NA
  
  
  ######
  #Community ages
  ######
  
  #extrac ages, function from Liam Revell
  n<-length(lista_arboles[[i]]$tip.label)
  Ages<-data.frame(setNames(lista_arboles[[i]]$edge.length[sapply(1:n,function(x,y) which(y==x),y=lista_arboles[[i]]$edge[,2])],lista_arboles[[i]]$tip.label))
  colnames(Ages)=c('Age')
  
  community.ages = community
  
  for (k in 1:dim(community.ages)[2])
    
  {
    print(k)
    sp.match = match(colnames(community.ages[k]),rownames(Ages))
    community.ages[which(community.ages[,k]>0),k] <- Ages[sp.match,1]
  }
  
  #community.ages=do.call(pmax, c(community.ages, list(na.rm=TRUE)))
  community.ages=rowMedians(as.matrix(community.ages), na.rm=T)
  communityAges=altura
  communityAges[]<-NA
  communityAges[]<-community.ages
  communityAges[which(is.na(Elev[,1]))] <- NA
  
  
  #communityAges=altura
  #communityAges[]<-NA
  #for (k in 1:dim(community.ages)[1])
  #
  #{
  #	communityAges[k]<-quantile(community.ages[k,], probs=0.75, na.rm=T)
  #}
  #communityAges[which(is.na(Elev[,1]))] <- NA
  
  
  
  
  ######
  # Beta-Diversity
  ######	
  
  Beta_Diversity=altura
  Beta_Diversity[]<-NA
  
  Cells=which(values(Riqueza)>4)
  
  for (h in 1:length(Cells))
    
  {
    
    print(paste(((round(h/length(Cells), 2))*100), "%", sep=" "))
    
    vecinos=adjacent(altura, Cells[h], directions=8)
    
    Vec=c()
    
    for (p in 1:nrow(vecinos))
      
    {
      AA=vecinos[p,]
      
      Vec=c(Vec, betadiver(community_2[AA,], method="e")[1])
      
    }
    
    Beta_Diversity[Cells[h]]<-mean(Vec)
    
    
  }
  
  
  
  ######
  # Phylo_beta-Diversity
  ######	
  
  Phylo_Beta_Diversity=altura
  Phylo_Beta_Diversity[]<-NA
  
  
  TabNombres=t(community_2)
  TabNombres=TabNombres[,1]
  
  TabNombres=treedata(lista_arboles[[i]], TabNombres)
  
  newPhy=TabNombres$phy
  
  
  Cells=which(values(Riqueza)>4)
  
  for (h in 1:length(Cells))
    
  {
    
    print(paste(((round(h/length(Cells), 2))*100), "%", sep=" "))
    
    vecinos=adjacent(altura, Cells[h], directions=8)
    
    Vec=c()
    
    for (p in 1:nrow(vecinos))
      
    {
      AA=vecinos[p,]
      
      if(is.na(values(Riqueza)[AA[1]]) | is.na(values(Riqueza)[AA[2]])) {next}
      if(values(Riqueza)[AA[1]] < 4 | values(Riqueza)[AA[2]] < 4) {next}
      
      if(values(Riqueza)[AA[1]] > 3 | values(Riqueza)[AA[2]] > 3)
      {
        
        Vec=c(Vec, comdist(community_2[AA,], as.dist(cophenetic(newPhy)))[1])
        #Vec=c(Vec, phylosor(community_2[AA,], newPhy)[1])
        
      }
      
    }
    
    Phylo_Beta_Diversity[Cells[h]]<-mean(Vec)
    
    
  }
  
  
  ######
  #Beta-Diversification
  ######
  
  
  dis=t(outer(DR[,1], DR[,1], "-"))
  
  colnames(dis)=rownames(DR)
  
  rownames(dis)=rownames(DR)
  
  dis=abs(dis)
  
  Cells=which(values(Riqueza)>4)
  
  BetaDiversification=altura
  BetaDiversification[]<-NA
  
  for (h in 1:length(Cells))
    
  {
    
    print(paste(((round(h/length(Cells), 2))*100), "%", sep=" "))
    
    vecinos=adjacent(altura, Cells[h], directions=8)
    
    Vec=c()
    
    for (p in 1:nrow(vecinos))
      
    {
      AA=vecinos[p,]
      
      Vec=c(Vec, comdist(community_2[AA,], dis)[1])
      #phylosor(community3, lista_arboles[[i]])
      
    }
    
    
      BetaDiversification[Cells[h]]<-max(Vec)
    
  }
  
  
  ######
  # Morphological disparity
  ######
  
  
  rasgos=tabla[match(colnames(community), tabla[,Grupos[i]]), c(Grupos[i], "Beak.Length_Culmen", "Beak.Width", "Beak.Depth", "Tarsus.Length", "Wing.Length", "Hand.Wing.Index", "Tail.Length", "Mass")]
  
  rownames(rasgos)=rasgos[,1]
  rasgos=rasgos[,-1]
  rasgos=na.omit(rasgos)
  rasgos=scale(rasgos)
  
  com=community
  
  differences<-setdiff(colnames(com), rownames(rasgos))
  
  com=com[,!(names(com) %in% differences)]
  
  rasgos <- rasgos[order(row.names(rasgos)), ]
  
  distFunc <- dist(rasgos)
  distFunc <- as.matrix(distFunc)
  
  com[is.na(com)]<-0
  com=com[,sort(colnames(com))]
  
  Cells=which(values(Riqueza)>3)
  
  mpd.result<-as.data.frame(mpd(com[Cells,], distFunc))
  colnames(mpd.result)=c("mpd")
  
  MorphoDisp=altura
  MorphoDisp[]<-NA
  
  MorphoDisp[Cells]<-log(mpd.result$mpd)
  
  
  #################################################
  
  
  Mapas[[i]][[1]]<-Riqueza
  Mapas[[i]][[2]]<-Diversificacion_ES
  Mapas[[i]][[3]]<-Diversificacion_ClaDS
  Mapas[[i]][[4]]<-Beta_Diversity
  Mapas[[i]][[5]]<-Phylo_Beta_Diversity
  Mapas[[i]][[6]]<-BetaDiversification
  Mapas[[i]][[7]]<-communityAges
  Mapas[[i]][[8]]<-MorphoDisp

  
  Final=cbind(coordinates(altura),as.data.frame(altura), as.data.frame(altura3), as.data.frame(Riqueza), as.data.frame(Diversificacion_ES), as.data.frame(Diversificacion_ClaDS), as.data.frame(Beta_Diversity), as.data.frame(Phylo_Beta_Diversity), as.data.frame(BetaDiversification), as.data.frame(communityAges), as.data.frame(MorphoDisp))
  
  
  colnames(Final)=c('x', 'y', 'Elevation', "Mountains", 'Richness', 'Diversification_ES', 'Diversification_ClaDS', 'Beta_Diversity', 'Phylo_Beta_Diversity', 'BetaDiversification', 'Time', 'MorpoDisp')
  
  #Final=Final[-which(is.na(Final$Elevation)),]
  Final=Final[which(Final$Richness>3),]
  Final2=Final[which(!is.na(Final$Mountains)),]
  
  #Final[,11]<-residuals(lm(Final[which(!is.na(Final$BetaDiversification)),'BetaDiversification']~Final[which(!is.na(Final$BetaDiversification)), 'Diversification'])) ##Esta columna se desactiva porque en este momento no necesito los reisduales de la betadiversificación
  
  #colnames(Final)=c('x', 'y', 'Elevation','Richness', 'Diversification', 'Beta_Diversity', 'BetaDiversification', 'Time', 'MorpoDisp', 'ResidualsBD')
  
  
  #Res_BetaDiversification=altura
  #Res_BetaDiversification[]<-NA
  #Res_BetaDiversification[as.numeric(rownames(Final))]<-Final$ResidualsBD
  

  setwd("/Users/elkintenorio/Downloads/Results/")
  
  Finales[[i]]<-Final
  Mountains[[i]]<-Final2
  write.table(Final, paste(Nombre_Grupo[i], 'RasterValues.txt', sep="_"))
  library(PerformanceAnalytics)
  pdf(paste(Nombre_Grupo[i], 'Correlations.pdf',sep='_'))
  chart.Correlation(Mountains[[i]][,c(5,6,9,11,12)], histogram=FALSE, pch=19)
  chart.Correlation(Mountains[[i]][,5:12], histogram=FALSE, pch=19)
  dev.off()
  
  Comunidades[[i]]<-community
  write.table(community, paste(Nombre_Grupo[i], 'Community.txt', sep="_"))
  
  writeRaster(Riqueza, paste(Nombre_Grupo[i], 'Riqueza.asc', sep='_'), format="ascii")
  writeRaster(Diversificacion_ES, paste(Nombre_Grupo[i], 'Diversificacion_ES.asc', sep='_'), format="ascii")
  writeRaster(Diversificacion_ClaDS, paste(Nombre_Grupo[i], 'Diversificacion_ClaDS.asc', sep='_'), format="ascii")
  writeRaster(Beta_Diversity, paste(Nombre_Grupo[i], 'Beta_Diversity.asc', sep='_'), format="ascii")
  writeRaster(Phylo_Beta_Diversity, paste(Nombre_Grupo[i], 'PhyloBeta_Diversity.asc', sep='_'), format="ascii")
  writeRaster(BetaDiversification, paste(Nombre_Grupo[i], 'BetaDiversification.asc', sep='_'), format="ascii")
  #writeRaster(Res_BetaDiversification, paste(Nombre_Grupo[i], 'Res_BetaDiversification.asc', sep='_'), format="ascii")
  writeRaster(communityAges, paste(Nombre_Grupo[i], 'communityAges.asc', sep='_'), format="ascii")
  writeRaster(MorphoDisp, paste(Nombre_Grupo[i], 'MorphoDisp.asc', sep='_'), format="ascii")
  
  
  pdf(paste(Nombre_Grupo[i], "1", 'matrix_plot.pdf',sep='_'))
  par(mfrow=c(2,4), oma=c(0,0,1,0))
  plot(Riqueza, main='Species richness')
  plot(America, add=T, lwd=0.1)
  plot(Diversificacion_ES, main='Diversification rates (ES)')
  plot(America, add=T, lwd=0.1)
  plot(Diversificacion_ClaDS, main='Diversification rates (ClaDS)')
  plot(America, add=T, lwd=0.1)
  plot(BetaDiversification, main='Beta Diversification')
  plot(America, add=T, lwd=0.1)
  plot(Beta_Diversity, main='Beta Diversity')
  plot(America, add=T, lwd=0.1)
  plot(Phylo_Beta_Diversity, main='PhyloBeta Diversity')
  plot(America, add=T, lwd=0.1)
  plot(MorphoDisp, main='Morphological disparity')
  plot(America, add=T, lwd=0.1)
  plot(communityAges, main='community Ages')
  plot(America, add=T, lwd=0.1)
  dev.off()
  
  pdf(paste(Nombre_Grupo[i], "2", 'matrix_plot.pdf',sep='_'))
  plot(Riqueza, main='Species richness')
  plot(America, add=T, lwd=0.1)
  plot(Diversificacion_ES, main='Diversification rates (ES)')
  plot(America, add=T, lwd=0.1)
  plot(Diversificacion_ClaDS, main='Diversification rates (ClaDS)')
  plot(America, add=T, lwd=0.1)
  plot(BetaDiversification, main='Beta Diversification')
  plot(America, add=T, lwd=0.1)
  plot(Beta_Diversity, main='Beta Diversity')
  plot(America, add=T, lwd=0.1)
  plot(Phylo_Beta_Diversity, main='PhyloBeta Diversity')
  plot(America, add=T, lwd=0.1)
  plot(MorphoDisp, main='Morphological disparity')
  plot(America, add=T, lwd=0.1)
  plot(communityAges, main='community Ages')
  plot(America, add=T, lwd=0.1)
  dev.off()
  
}


#########################
### Plots with GGplot ###
#########################


#########################################
#    Species richness vs. Elevation     #
#########################################

L_plot_Tro_DivRates <- ggplot(Mountains[[1]], aes(x=Elevation, y=Richness)) + geom_point(size = 1.5, alpha = 0.7) + geom_smooth(stat= "smooth", method = "gam", color='blue', size = 0.3, formula = y ~ s(x)) + theme_classic() + theme(axis.title = element_text(size = 22), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 15)) + xlab("Elevation (m)") + ylab("Number of species") + coord_cartesian(xlim =c(0, 4500))
L_plot_Sub_DivRates <- ggplot(Mountains[[2]], aes(x=Elevation, y=Richness)) + geom_point(size = 1.5, alpha = 0.7) + geom_smooth(stat= "smooth", method = "gam", color='green', size = 0.3, formula = y ~ s(x)) + theme_classic() + theme(axis.title = element_text(size = 22), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 15)) + xlab("Elevation (m)") + ylab("Number of species") + coord_cartesian(xlim =c(0, 4500))
L_plot_Osc_DivRates <- ggplot(Mountains[[3]], aes(x=Elevation, y=Richness)) + geom_point(size = 1.5, alpha = 0.7) + geom_smooth(stat= "smooth", method = "gam", color='red', size = 0.3, formula = y ~ s(x)) + theme_classic() + theme(axis.title = element_text(size = 22), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 15))+ xlab("Elevation (m)") + ylab("Number of species") + coord_cartesian(xlim =c(0, 4500))

##########################
#    Multiplot Mapas     #
##########################

pdf("Multiplot_Mapas.pdf", height=566/72, width=800/72)

library(RColorBrewer)
cols <- brewer.pal(9, "BuGn")
pal <- colorRampPalette(cols)

blue.col <- colorRampPalette(c("lightblue", "royalblue"))

par(mfrow=c(4,3), mar=c(1,1,1,1), oma=c(1,1,0,0.5))

plantilla=crop(Mapas[[1]][[1]], c(-117, -28, -23, 23))
plantilla[]<-NA

plot(plantilla, xaxt='n', ann=FALSE)
plot(Mapas[[1]][[2]], col=pal(20), axes = FALSE, add=T)
plot(Montanas, add=T, lwd=0.5)
plot(America, add=T, lwd=0.5)

plot(plantilla, xaxt='n', yaxt='n', ann=FALSE)
plot(Mapas[[2]][[2]], col=pal(20), add=T)
plot(Montanas, add=T, lwd=0.5)
plot(America, add=T, lwd=0.5)

plot(plantilla, xaxt='n', yaxt='n', ann=FALSE)
plot(Mapas[[3]][[2]], col=pal(20), add=T)
plot(Montanas, add=T, lwd=0.5)
plot(America, add=T, lwd=0.5)

plot(plantilla, xaxt='n', xaxt='n', ann=FALSE)
plot(Mapas[[1]][[7]], col=topo.colors(100), add=T) #heat.colors(20)
plot(Montanas, add=T, lwd=0.5)
plot(America, add=T, lwd=0.5)

plot(plantilla, xaxt='n', yaxt='n', ann=FALSE)
plot(Mapas[[2]][[7]], col=topo.colors(100), add=T)
plot(Montanas, add=T, lwd=0.5)
plot(America, add=T, lwd=0.5)

plot(plantilla, xaxt='n', yaxt='n', ann=FALSE)
plot(Mapas[[3]][[7]], col=topo.colors(100), add=T)
plot(Montanas, add=T, lwd=0.5)
plot(America, add=T, lwd=0.5)

plot(plantilla, xaxt='n', ann=FALSE)
plot(Mapas[[1]][[8]], col=rainbow(100), add=T)
plot(Montanas, add=T, lwd=0.5)
plot(America, add=T, lwd=0.5)

plot(plantilla, xaxt='n', yaxt='n', ann=FALSE)
plot(Mapas[[2]][[8]], col=rainbow(100), add=T)
plot(Montanas, add=T, lwd=0.5)
plot(America, add=T, lwd=0.5)

plot(plantilla, xaxt='n', yaxt='n', ann=FALSE)
plot(Mapas[[3]][[8]], col=rainbow(100), add=T)
plot(Montanas, add=T, lwd=0.5)
plot(America, add=T, lwd=0.5)

plot(plantilla)
plot(Mapas[[1]][[5]], col=terrain.colors(10), add=T)
plot(Montanas, add=T, lwd=0.5)
plot(America, add=T, lwd=0.5)

plot(plantilla, yaxt='n', ann=FALSE)
plot(Mapas[[2]][[5]], col=terrain.colors(10), add=T)
plot(Montanas, add=T, lwd=0.5)
plot(America, add=T, lwd=0.5)

plot(plantilla, yaxt='n', ann=FALSE)
plot(Mapas[[3]][[5]], col=terrain.colors(10), add=T)
plot(Montanas, add=T, lwd=0.5)
plot(America, add=T, lwd=0.5)

dev.off()


################################
#   Multiplot Figure 4 and 5   #
################################


pdf("Figure4.pdf", height=305/72, width=600/72)
#div rates
L_plot_Tro_DivRates <- ggplot(Mountains[[1]], aes(x=Elevation, y=Diversification_ES)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='blue', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 0)),  plot.margin=unit(c(0.3,0,-0.5,0.1), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("Diversification rates") + coord_cartesian(xlim =c(0, 4500), ylim = c(-1.8, -1.1))
L_plot_Sub_DivRates<- ggplot(Mountains[[2]], aes(x=Elevation, y=Diversification_ES)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='green', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)),  plot.margin=unit(c(0.3,0,-0.5,-0.8), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("") + coord_cartesian(xlim =c(0, 4500), ylim = c(-1.8, -1.1))
L_plot_Osc_DivRates <- ggplot(Mountains[[3]], aes(x=Elevation, y=Diversification_ES)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='red', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)),  plot.margin=unit(c(0.3,0.2,-0.5,-0.8), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("") + coord_cartesian(xlim =c(0, 4500), ylim = c(-1.8, -1.1))
#Community ages
L_plot_Tro_Time <- ggplot(Mountains[[1]], aes(x=Elevation, y=Time)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='blue', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 14.5, b = 0, l = -2.5)),  plot.margin=unit(c(0,0,-0.5,0.2), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("Average community age") + coord_cartesian(xlim =c(0, 4500), ylim = c(1.0,6.0))
L_plot_Sub_Time <- ggplot(Mountains[[2]], aes(x=Elevation, y=Time)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='green', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 23.5, b = 0, l = 0)), plot.margin=unit(c(0,0,-0.5,-0.8), "cm"),  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("") + coord_cartesian(xlim =c(0, 4500), ylim = c(1.0,6.0))
L_plot_Osc_Time <- ggplot(Mountains[[3]], aes(x=Elevation, y=Time)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='red', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12),  legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 23.5, b = 0, l = 0)),  plot.margin=unit(c(0,0.2,-0.5,-0.8), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("") + coord_cartesian(xlim =c(0, 4500), ylim = c(1.0,6.0))
grid.arrange(L_plot_Tro_DivRates, L_plot_Sub_DivRates, L_plot_Osc_DivRates, L_plot_Tro_Time, L_plot_Sub_Time, L_plot_Osc_Time, ncol = 3, nrow = 2)
dev.off()

pdf("Figure4_Supple.pdf", height=305/72, width=600/72)
#div rates
L_plot_Tro_DivRates <- ggplot(Mountains[[1]], aes(x=Elevation, y=Diversification_ClaDS)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='blue', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 0)),  plot.margin=unit(c(0.3,0,-0.5,0.1), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("Diversification rates") + coord_cartesian(xlim =c(0, 4500), ylim = c(0.10, 0.353))
L_plot_Sub_DivRates<- ggplot(Mountains[[2]], aes(x=Elevation, y=Diversification_ClaDS)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='green', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)),  plot.margin=unit(c(0.3,0,-0.5,-0.8), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("") + coord_cartesian(xlim =c(0, 4500), ylim = c(0.10, 0.353))
L_plot_Osc_DivRates <- ggplot(Mountains[[3]], aes(x=Elevation, y=Diversification_ClaDS)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='red', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)),  plot.margin=unit(c(0.3,0.2,-0.5,-0.8), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("") + coord_cartesian(xlim =c(0, 4500), ylim = c(0.10, 0.353))
#Community ages
L_plot_Tro_Time <- ggplot(Mountains[[1]], aes(x=Elevation, y=Time)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='blue', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 14.5, b = 0, l = -2.5)),  plot.margin=unit(c(0,0,-0.5,0.2), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("Average community age") + coord_cartesian(xlim =c(0, 4500), ylim = c(1.0,6.0))
L_plot_Sub_Time <- ggplot(Mountains[[2]], aes(x=Elevation, y=Time)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='green', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 23.5, b = 0, l = 0)), plot.margin=unit(c(0,0,-0.5,-0.8), "cm"),  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("") + coord_cartesian(xlim =c(0, 4500), ylim = c(1.0,6.0))
L_plot_Osc_Time <- ggplot(Mountains[[3]], aes(x=Elevation, y=Time)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='red', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12),  legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 23.5, b = 0, l = 0)),  plot.margin=unit(c(0,0.2,-0.5,-0.8), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("") + coord_cartesian(xlim =c(0, 4500), ylim = c(1.0,6.0))
grid.arrange(L_plot_Tro_DivRates, L_plot_Sub_DivRates, L_plot_Osc_DivRates, L_plot_Tro_Time, L_plot_Sub_Time, L_plot_Osc_Time, ncol = 3, nrow = 2)
dev.off()


pdf("Figure5.pdf", height=305/72, width=600/72)
#Morpho
L_plot_Tro_Morph <- ggplot(Mountains[[1]], aes(x=Elevation, y=MorpoDisp)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='blue', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 4, b = 0, l = 0)),  plot.margin=unit(c(0.3,0,-0.5,0.1), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("Morphological disparity") + coord_cartesian(xlim =c(0, 4500), ylim = c(0.5,2))
L_plot_Sub_Morph <- ggplot(Mountains[[2]], aes(x=Elevation, y=MorpoDisp)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='green', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)),  plot.margin=unit(c(0.3,0,-0.5,-0.8), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("") + coord_cartesian(xlim =c(0, 4500), ylim = c(0.5, 2))
L_plot_Osc_Morph <- ggplot(Mountains[[3]], aes(x=Elevation, y=MorpoDisp)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='red', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)),  plot.margin=unit(c(0.3,0.2,-0.5,-0.8), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("") + coord_cartesian(xlim =c(0, 4500), ylim = c(0.5, 2))

#Beta-Div
L_plot_Tro_Beta <- ggplot(Mountains[[1]], aes(x=Elevation, y=Beta_Diversity)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='blue', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 3.0, b = 0, l =10)),  plot.margin=unit(c(0,0,-0.5,-0.2), "cm"),  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("Phylobeta diversity") + coord_cartesian(xlim =c(0, 4500), ylim = c(0,0.8))
L_plot_Sub_Beta <- ggplot(Mountains[[2]], aes(x=Elevation, y=Beta_Diversity)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='green', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 13.5, b = 0, l = 10)),  plot.margin=unit(c(0,0,-0.5,-1.2), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("Elevation (m)") + ylab("")  + coord_cartesian(xlim =c(0, 4500), ylim = c(0,0.8))
L_plot_Osc_Beta <- ggplot(Mountains[[3]], aes(x=Elevation, y=Beta_Diversity)) + geom_point(size = 0.5, alpha = 0.15) + geom_smooth(stat= "smooth", method = "gam", color='red', size = 0.3, formula = y ~ s(x)) + theme_light() + theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10), legend.title = element_text(size = 15), axis.text = element_text(size = 8), axis.title.y = element_text(margin = margin(t = 0, r = 13.5, b = 0, l = 10)),  plot.margin=unit(c(0,0.2,-0.5,-1.2), "cm"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) + xlab("") + ylab("") + coord_cartesian(xlim =c(0, 4500), ylim =  c(0,0.8))
grid.arrange(L_plot_Tro_Morph, L_plot_Sub_Morph, L_plot_Osc_Morph, L_plot_Tro_Beta, L_plot_Sub_Beta, L_plot_Osc_Beta, ncol = 3, nrow = 2)
dev.off()


#######################################
### Scale variables and set to zero ###
#######################################

tabla1=as.data.frame(scale(data.frame(Mountains[[1]]$Richness, Mountains[[1]]$Elevation, Mountains[[1]]$Diversification_ES, Mountains[[1]]$Diversification_ClaDS, Mountains[[1]]$Time, Mountains[[1]]$MorpoDisp, Mountains[[1]]$Beta_Diversity)))
colnames(tabla1)=c('Richness', 'Elevation_Scale', 'Diversification_ES', 'Diversification_ClaDS', 'Time', 'MorpoDisp', 'Beta_Diversity')
tabla1=data.frame(Richness_Scale=tabla1$Richness+abs(min(tabla1$Richness, na.rm=T)), Elevation_Scale=tabla1$Elevation+abs(min(tabla1$Elevation, na.rm=T)), Diversification_ES_Scale=tabla1$Diversification_ES+abs(min(tabla1$Diversification_ES, na.rm=T)), Diversification_ClaDS_Scale=tabla1$Diversification_ClaDS+abs(min(tabla1$Diversification_ClaDS, na.rm=T)), Time_Scale=tabla1$Time+abs(min(tabla1$Time, na.rm=T)), MorpoDisp_Scale=tabla1$MorpoDisp+abs(min(tabla1$MorpoDisp, na.rm=T)), Beta_Scale=tabla1$Beta_Diversity+abs(min(tabla1$Beta_Diversity, na.rm=T)))
tabla1=data.frame(tabla1, Richness_ScaleLog=log(tabla1$Richness_Scale+0.01)+abs(min(log(tabla1$Richness_Scale+0.01), na.rm=T)), Elevation_ScaleLog=log(tabla1$Elevation_Scale+0.01)+abs(min(log(tabla1$Elevation_Scale+0.01), na.rm=T)), Diversification_SE_ScaleLog=log(tabla1$Diversification_ES_Scale+0.01)+abs(min(log(tabla1$Diversification_ES_Scale+0.01), na.rm=T)), Diversification_ClaDS_ScaleLog=log(tabla1$Diversification_ClaDS_Scale+0.01)+abs(min(log(tabla1$Diversification_ClaDS_Scale+0.01), na.rm=T)), Time_ScaleLog=log(tabla1$Time_Scale+0.01)+abs(min(log(tabla1$Time_Scale+0.01), na.rm=T)), MorpoDisp_ScaleLog=log(tabla1$MorpoDisp_Scale+0.01)+abs(min(log(tabla1$MorpoDisp_Scale+0.01), na.rm=T)), Beta_ScaleLog=log(tabla1$Beta_Scale+0.01)+abs(min(log(tabla1$Beta_Scale+0.01), na.rm=T)))
tabla1=data.frame(Mountains[[1]], tabla1)
tabla1=na.omit(tabla1)

tabla2=as.data.frame(scale(data.frame(Mountains[[2]]$Richness, Mountains[[2]]$Elevation, Mountains[[2]]$Diversification_ES, Mountains[[2]]$Diversification_ClaDS, Mountains[[2]]$Time, Mountains[[2]]$MorpoDisp, Mountains[[2]]$Beta_Diversity)))
colnames(tabla2)=c('Richness', 'Elevation_Scale', 'Diversification_ES', 'Diversification_ClaDS', 'Time', 'MorpoDisp', 'Beta_Diversity')
tabla2=data.frame(Richness_Scale=tabla2$Richness+abs(min(tabla2$Richness, na.rm=T)), Elevation_Scale=tabla2$Elevation+abs(min(tabla2$Elevation, na.rm=T)), Diversification_ES_Scale=tabla2$Diversification_ES+abs(min(tabla2$Diversification_ES, na.rm=T)), Diversification_ClaDS_Scale=tabla2$Diversification_ClaDS+abs(min(tabla2$Diversification_ClaDS, na.rm=T)), Time_Scale=tabla2$Time+abs(min(tabla2$Time, na.rm=T)), MorpoDisp_Scale=tabla2$MorpoDisp+abs(min(tabla2$MorpoDisp, na.rm=T)), Beta_Scale=tabla2$Beta_Diversity+abs(min(tabla2$Beta_Diversity, na.rm=T)))
tabla2=data.frame(tabla2, Richness_ScaleLog=log(tabla2$Richness_Scale+0.01)+abs(min(log(tabla2$Richness_Scale+0.01), na.rm=T)), Elevation_ScaleLog=log(tabla2$Elevation_Scale+0.01)+abs(min(log(tabla2$Elevation_Scale+0.01), na.rm=T)), Diversification_SE_ScaleLog=log(tabla2$Diversification_ES_Scale+0.01)+abs(min(log(tabla2$Diversification_ES_Scale+0.01), na.rm=T)), Diversification_ClaDS_ScaleLog=log(tabla2$Diversification_ClaDS_Scale+0.01)+abs(min(log(tabla2$Diversification_ClaDS_Scale+0.01), na.rm=T)), Time_ScaleLog=log(tabla2$Time_Scale+0.01)+abs(min(log(tabla2$Time_Scale+0.01), na.rm=T)), MorpoDisp_ScaleLog=log(tabla2$MorpoDisp_Scale+0.01)+abs(min(log(tabla2$MorpoDisp_Scale+0.01), na.rm=T)), Beta_ScaleLog=log(tabla2$Beta_Scale+0.01)+abs(min(log(tabla2$Beta_Scale+0.01), na.rm=T)))
tabla2=data.frame(Mountains[[2]], tabla2)
tabla2=na.omit(tabla2)

tabla3=as.data.frame(scale(data.frame(Mountains[[3]]$Richness, Mountains[[3]]$Elevation, Mountains[[3]]$Diversification_ES, Mountains[[3]]$Diversification_ClaDS, Mountains[[3]]$Time, Mountains[[3]]$MorpoDisp, Mountains[[3]]$Beta_Diversity)))
colnames(tabla3)=c('Richness', 'Elevation_Scale', 'Diversification_ES', 'Diversification_ClaDS', 'Time', 'MorpoDisp', 'Beta_Diversity')
tabla3=data.frame(Richness_Scale=tabla3$Richness+abs(min(tabla3$Richness, na.rm=T)), Elevation_Scale=tabla3$Elevation+abs(min(tabla3$Elevation, na.rm=T)), Diversification_ES_Scale=tabla3$Diversification_ES+abs(min(tabla3$Diversification_ES, na.rm=T)), Diversification_ClaDS_Scale=tabla3$Diversification_ClaDS+abs(min(tabla3$Diversification_ClaDS, na.rm=T)), Time_Scale=tabla3$Time+abs(min(tabla3$Time, na.rm=T)), MorpoDisp_Scale=tabla3$MorpoDisp+abs(min(tabla3$MorpoDisp, na.rm=T)), Beta_Scale=tabla3$Beta_Diversity+abs(min(tabla3$Beta_Diversity, na.rm=T)))
tabla3=data.frame(tabla3, Richness_ScaleLog=log(tabla3$Richness_Scale+0.01)+abs(min(log(tabla3$Richness_Scale+0.01), na.rm=T)), Elevation_ScaleLog=log(tabla3$Elevation_Scale+0.01)+abs(min(log(tabla3$Elevation_Scale+0.01), na.rm=T)), Diversification_SE_ScaleLog=log(tabla3$Diversification_ES_Scale+0.01)+abs(min(log(tabla3$Diversification_ES_Scale+0.01), na.rm=T)), Diversification_ClaDS_ScaleLog=log(tabla3$Diversification_ClaDS_Scale+0.01)+abs(min(log(tabla3$Diversification_ClaDS_Scale+0.01), na.rm=T)), Time_ScaleLog=log(tabla3$Time_Scale+0.01)+abs(min(log(tabla3$Time_Scale+0.01), na.rm=T)), MorpoDisp_ScaleLog=log(tabla3$MorpoDisp_Scale+0.01)+abs(min(log(tabla3$MorpoDisp_Scale+0.01), na.rm=T)), Beta_ScaleLog=log(tabla3$Beta_Scale+0.01)+abs(min(log(tabla3$Beta_Scale+0.01), na.rm=T)))
tabla3=data.frame(Mountains[[3]], tabla3)
tabla3=na.omit(tabla3)

tablas<-list(tabla1, tabla2, tabla3)


##################
# Spatial regression with Multivariate Adaptive Regression Splines (MARS)
##################

library(earth)

pdf('MARS_Plot.pdf')
sink("MARSs.txt")

for (i in 1:3)
  
{

list_mars<-list()

list_mars[[1]] <- earth(Diversification_ES_Scale ~ Elevation_Scale, data = tablas[[i]], nk = 9, degree = 1) ### degree = 2 is suggested if interactions want to be considered, but no more
list_mars[[2]] <- earth(Diversification_ClaDS_Scale ~ Elevation_Scale, data = tablas[[i]], nk = 9, degree = 1) ### degree = 2 is suggested if interactions want to be considered, but no more
list_mars[[3]] <- earth(Time_Scale ~ Elevation_Scale, data = tablas[[i]], nk = 9, degree = 1) ### degree = 2 is suggested if interactions want to be considered, but no more
list_mars[[4]] <- earth(MorpoDisp_Scale ~ Elevation_Scale, data = tablas[[i]], nk = 9, degree = 1) ### degree = 2 is suggested if interactions want to be considered, but no more
list_mars[[5]] <- earth(Beta_Scale ~ Elevation_Scale, data = tablas[[i]], nk = 9, degree = 1) ### degree = 2 is suggested if interactions want to be considered, but no more

#mars2 <- earth(Div ~ Elevation, data = tabla1, degree = 1, glm=list(family=Gamma))

print(summary(list_mars[[1]]))
print(summary(list_mars[[2]]))
print(summary(list_mars[[3]]))
print(summary(list_mars[[4]]))
print(summary(list_mars[[5]]))

Plot1=ggplot(tablas[[i]], aes(Elevation_Scale, Diversification_ES_Scale)) +
  geom_point(size = 1, alpha = .2) +
  geom_line(aes(y = as.vector(list_mars[[1]]$fitted.values)), size = 1, color = "blue") +
  ggtitle("")

Plot2=ggplot(tablas[[i]], aes(Elevation_Scale, Diversification_ClaDS_Scale)) +
  geom_point(size = 1, alpha = .2) +
  geom_line(aes(y = as.vector(list_mars[[2]]$fitted.values)), size = 1, color = "blue") +
  ggtitle("")


Plot3=ggplot(tablas[[i]], aes(Elevation_Scale, Time_Scale)) +
  geom_point(size = 1, alpha = .2) +
  geom_line(aes(y = as.vector(list_mars[[3]]$fitted.values)), size = 1, color = "blue") +
  ggtitle("")

Plot4=ggplot(tablas[[i]], aes(Elevation_Scale, MorpoDisp_Scale)) +
  geom_point(size = 1, alpha = .2) +
  geom_line(aes(y = as.vector(list_mars[[4]]$fitted.values)), size = 1, color = "blue") +
  ggtitle("")

Plot5=ggplot(tablas[[i]], aes(Elevation_Scale, Beta_Scale)) +
  geom_point(size = 1, alpha = .2) +
  geom_line(aes(y = as.vector(list_mars[[5]]$fitted.values)), size = 1, color = "blue") +
  ggtitle("")

grid.arrange(Plot1, Plot2, Plot3, Plot4, Plot5, ncol = 2, nrow = 3)


}

dev.off()
sink()


#######################################
# Structural Equation Modeling (SEM) ##
#######################################

library(lavaan)
library(semPlot)
library(knitr)
library(sesem)
library(mgcv)
library(fields)
library(raster)


################
# With Equal Split estimations
################

Complete <- '
         
         # Direct models
         Richness_Scale ~ e1 * Beta_Scale
         Richness_Scale ~ e2 * MorpoDisp_Scale

         Richness_Scale ~ e3 * Diversification_ES_Scale
         Richness_Scale ~ e4 * Time_Scale
         
         # Mediation models
         Beta_Scale ~ a * Diversification_ES_Scale
         MorpoDisp_Scale ~ b * Diversification_ES_Scale
         Beta_Scale ~ c * Time_Scale
         MorpoDisp_Scale ~ d * Time_Scale

         # Indirect effects
         ae1 := a * e1
         be2 := b * e2
         ce1 := c * e1
         de2 := d * e2
         
         # Total effects
         totalae1 := e3 + (ae1)
         totalbe2 := e3 + (be2)
         totalce1 := e4 + (ce1)
         totalde2 := e4 + (de2)'

#################################
##  Total
#################################

A=match(rownames(tablas[[1]]), rownames(tablas[[2]]))
B=match(rownames(tablas[[1]]), rownames(tablas[[3]]))

Richness_Scale=tablas[[1]]$Richness_Scale+
Diversification_ES
Elevation_Scale
Diversification_ES_Scale
Diversification_ClaDS_Scale
Time_Scale
MorpoDisp_Scale

Total_Tabla=data.frame()





sink("SEM_Total_WithES.txt")

set.seed(1234)

fComplete <- sem(Complete, data = Total_Tabla, se = "bootstrap", bootstrap = 1000)
summary(fComplete, standardized = TRUE, rsquare=T)
fitmeasures(fComplete)[c("chisq", "pvalue", "aic")]
parameterestimates(fComplete, boot.ci.type = "bca.simple", standardized = TRUE)

sink()


#################################
##  Trochilidae
#################################

sink("SEM_Trochilidae_WithES.txt")

set.seed(1234)

fComplete <- sem(Complete, data = tablas[[1]], se = "bootstrap", bootstrap = 1000)
summary(fComplete, standardized = TRUE, rsquare=T)
fitmeasures(fComplete)[c("chisq", "pvalue", "aic")]
parameterestimates(fComplete, boot.ci.type = "bca.simple", standardized = TRUE)

print('Spatial SEM results')

#Generate a distance matrix for use in the spatial structural equation modeldistancematrix <- round(rdist.earth(x1 = tablas[[1]][,c('x', 'y')])[lower.tri(matrix(nrow = nrow(tablas[[1]]), ncol = nrow(tablas[[1]])), diag=T)])temp_bins<-make.bin(distancematrix, type="n.bins",n.bin=3, p.dist=25)binsize<-temp_bins[1][[1]]binname<-temp_bins[2][[1]]covariances<-make.covar(tablas[[1]], distancematrix, binsize, binname)#Fit the modeloutput<-sesem::runModels(Complete, covariances)modelsummary(output)
plotmodelfit(output,rmsea_err=FALSE)plotpath(output,pch=11)
par(mar=c(2.5,2.5,1.5,1))gam.path(output, plot.points=F,se.plot=T)bin.rsquare(output, bin="binflat")
sink()

#################################
##  Suboscines
#################################

sink("SEM_Suboscines_WithES.txt")

set.seed(1234)

fComplete <- sem(Complete, data = tablas[[2]], se = "bootstrap", bootstrap = 1000)
summary(fComplete, standardized = TRUE, rsquare=T)
fitmeasures(fComplete)[c("chisq", "pvalue", "aic")]
parameterestimates(fComplete, boot.ci.type = "bca.simple", standardized = TRUE)

print('Spatial SEM results')

#Generate a distance matrix for use in the spatial structural equation model
distancematrix <- round(rdist.earth(x1 = tablas[[2]][,c('x', 'y')])[lower.tri(matrix(nrow = nrow(tablas[[2]]), ncol = nrow(tablas[[2]])), diag=T)])
temp_bins<-make.bin(distancematrix, type="n.bins",n.bin=3, p.dist=25)
binsize<-temp_bins[1][[1]]
binname<-temp_bins[2][[1]]
covariances<-make.covar(tablas[[2]], distancematrix, binsize, binname)

#Fit the model
output<-sesem::runModels(Complete, covariances)
modelsummary(output)
plotmodelfit(output,rmsea_err=FALSE)
plotpath(output,pch=11)
par(mar=c(2.5,2.5,1.5,1))
gam.path(output, plot.points=F,se.plot=T)
bin.rsquare(output, bin="binflat")


sink()


#################################
##  Oscines
#################################

sink("SEM_Oscines_WithES.txt")

set.seed(1234)

fComplete <- sem(Complete, data = tablas[[3]], se = "bootstrap", bootstrap = 1000)
summary(fComplete, standardized = TRUE, rsquare=T)
fitmeasures(fComplete)[c("chisq", "pvalue", "aic")]
parameterestimates(fComplete, boot.ci.type = "bca.simple", standardized = TRUE)

print('Spatial SEM results')

#Generate a distance matrix for use in the spatial structural equation model
distancematrix <- round(rdist.earth(x1 = tablas[[3]][,c('x', 'y')])[lower.tri(matrix(nrow = nrow(tablas[[3]]), ncol = nrow(tablas[[3]])), diag=T)])
temp_bins<-make.bin(distancematrix, type="n.bins",n.bin=3, p.dist=25)
binsize<-temp_bins[1][[1]]
binname<-temp_bins[2][[1]]
covariances<-make.covar(tablas[[3]], distancematrix, binsize, binname)

#Fit the model
output<-sesem::runModels(Complete, covariances)
modelsummary(output)
plotmodelfit(output,rmsea_err=FALSE)
plotpath(output,pch=11)
par(mar=c(2.5,2.5,1.5,1))
gam.path(output, plot.points=F,se.plot=T)
bin.rsquare(output, bin="binflat")

sink()


################
# With ClaDS estimations
################

Complete <- '
         
         # Direct models
         Richness_Scale ~ e1 * Beta_Scale
         Richness_Scale ~ e2 * MorpoDisp_Scale

         Richness_Scale ~ e3 * Diversification_ClaDS_Scale
         Richness_Scale ~ e4 * Time_Scale
         
         # Mediation models
         Beta_Scale ~ a * Diversification_ClaDS_Scale
         MorpoDisp_Scale ~ b * Diversification_ClaDS_Scale
         Beta_Scale ~ c * Time_Scale
         MorpoDisp_Scale ~ d * Time_Scale

         # Indirect effects
         ae1 := a * e1
         be2 := b * e2
         ce1 := c * e1
         de2 := d * e2
         
         # Total effects
         totalae1 := e3 + (ae1)
         totalbe2 := e3 + (be2)
         totalce1 := e4 + (ce1)
         totalde2 := e4 + (de2)'

#################################
##  Total
#################################

sink("SEM_Total_WithClaDS.txt")

set.seed(1234)

fComplete <- sem(Complete, data = Total_Tabla, se = "bootstrap", bootstrap = 1000)
summary(fComplete, standardized = TRUE, rsquare=T)
fitmeasures(fComplete)[c("chisq", "pvalue", "aic")]
parameterestimates(fComplete, boot.ci.type = "bca.simple", standardized = TRUE)

sink()

#################################
##  Trochilidae
#################################

sink("SEM_Trochilidae_WithClaDS.txt")

set.seed(1234)

fComplete <- sem(Complete, data = tablas[[1]], se = "bootstrap", bootstrap = 1000)
summary(fComplete, standardized = TRUE, rsquare=T)
fitmeasures(fComplete)[c("chisq", "pvalue", "aic")]
parameterestimates(fComplete, boot.ci.type = "bca.simple", standardized = TRUE)

print('Spatial SEM results')

#Generate a distance matrix for use in the spatial structural equation model
distancematrix <- round(rdist.earth(x1 = tablas[[1]][,c('x', 'y')])[lower.tri(matrix(nrow = nrow(tablas[[1]]), ncol = nrow(tablas[[1]])), diag=T)])
temp_bins<-make.bin(distancematrix, type="n.bins",n.bin=3, p.dist=25)
binsize<-temp_bins[1][[1]]
binname<-temp_bins[2][[1]]
covariances<-make.covar(tablas[[1]], distancematrix, binsize, binname)

#Fit the model
output<-sesem::runModels(Complete, covariances)
modelsummary(output)
plotmodelfit(output,rmsea_err=FALSE)
plotpath(output,pch=11)
par(mar=c(2.5,2.5,1.5,1))
gam.path(output, plot.points=F,se.plot=T)
bin.rsquare(output, bin="binflat")

sink()


#################################
##  Suboscines
#################################

sink("SEM_Suboscines_WithClaDS.txt")

set.seed(1234)

fComplete <- sem(Complete, data = tablas[[2]], se = "bootstrap", bootstrap = 1000)
summary(fComplete, standardized = TRUE, rsquare=T)
fitmeasures(fComplete)[c("chisq", "pvalue", "aic")]
parameterestimates(fComplete, boot.ci.type = "bca.simple", standardized = TRUE)

print('Spatial SEM results')

#Generate a distance matrix for use in the spatial structural equation model
distancematrix <- round(rdist.earth(x1 = tablas[[2]][,c('x', 'y')])[lower.tri(matrix(nrow = nrow(tablas[[2]]), ncol = nrow(tablas[[2]])), diag=T)])
temp_bins<-make.bin(distancematrix, type="n.bins",n.bin=3, p.dist=25)
binsize<-temp_bins[1][[1]]
binname<-temp_bins[2][[1]]
covariances<-make.covar(tablas[[2]], distancematrix, binsize, binname)

#Fit the model
output<-sesem::runModels(Complete, covariances)
modelsummary(output)
plotmodelfit(output,rmsea_err=FALSE)
plotpath(output,pch=11)
par(mar=c(2.5,2.5,1.5,1))
gam.path(output, plot.points=F,se.plot=T)
bin.rsquare(output, bin="binflat")


sink()



#################################
##  Oscines
#################################

sink("SEM_Oscines_WithClaDS.txt")

set.seed(1234)

fComplete <- sem(Complete, data = tablas[[3]], se = "bootstrap", bootstrap = 1000)
summary(fComplete, standardized = TRUE, rsquare=T)
fitmeasures(fComplete)[c("chisq", "pvalue", "aic")]
parameterestimates(fComplete, boot.ci.type = "bca.simple", standardized = TRUE)

print('Spatial SEM results')

#Generate a distance matrix for use in the spatial structural equation model
distancematrix <- round(rdist.earth(x1 = tablas[[3]][,c('x', 'y')])[lower.tri(matrix(nrow = nrow(tablas[[3]]), ncol = nrow(tablas[[3]])), diag=T)])
temp_bins<-make.bin(distancematrix, type="n.bins",n.bin=3, p.dist=25)
binsize<-temp_bins[1][[1]]
binname<-temp_bins[2][[1]]
covariances<-make.covar(tablas[[3]], distancematrix, binsize, binname)

#Fit the model
output<-sesem::runModels(Complete, covariances)
modelsummary(output)
plotmodelfit(output,rmsea_err=FALSE)
plotpath(output,pch=11)
par(mar=c(2.5,2.5,1.5,1))
gam.path(output, plot.points=F,se.plot=T)
bin.rsquare(output, bin="binflat")

sink()


############
# SEMs using piecewiseSEM
############

##Code adapated from Garcia-Andrade et al. 2021. Evolutionary and environmental drivers of species richness in poeciliid fishes across the Americas. Global Ecology and Biogeography.

# Load packages
library(ncf)
library(spdep)
library(spatialreg) # version 1.1-3
library(piecewiseSEM)

# Set your working directory
dir_work <- ("~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/CUNY PhD/Thesis/Cap 1 Macro Montanas/Results/piecewiseSEM") 
setwd(dir_work)

### Model selection for SAR neighbourhood distances and weighting schemes ###
# This code fits a series of spatial autoregressive models,
# using different neighbour distances and weights scheme to find the best fit for each path 

# Matrix weighting schemes to prove
scheme <- c('W', 'C', 'S')

Clade_Folder <- c('Humms', 'Subosci', 'Osci')

Div_Folder <- c('DR', 'ClaDS')

for (j in 1:3) ##change among clades

{
    CF=Clade_Folder[j]

    for (z in 1:2) ##change among diversification estimations
    
    {
        # Make a matrix of spatial coordinates (X and Y coordinates)
        sp <- SpatialPoints(data.frame(x=tablas[[j]]$x, y=tablas[[j]]$y))
        crs(sp) <- '+proj=longlat +ellps=WGS84 +no_defs'
        coords <- coordinates(sp)

        # find min and maximum nearest neighbour
        k1 = knn2nb(knearneigh(coords, k=1))
        dist <- unlist(nbdists(k1, coords))
        max1 <- max(dist)
        min1 <- min(dist)

        # create a series of neighbour matrices based on different distances (distances are in km)
        d1 <- dnearneigh(sp, longlat = T, d1=0, d2=min1)
        d2 <- dnearneigh(sp, longlat = T, d1=0, d2=max1)

        # Create a data frame with path number and equations that will be tested in the pSEM
        path_equation <- data.frame(path_number = paste0('path',1:3))
        
        if(z==1)
        {path_equation$equation <-c('Richness_Scale ~ Diversification_ES_Scale + Time_Scale + Beta_Scale + MorpoDisp_Scale', 'Beta_Scale ~ Diversification_ES_Scale + Time_Scale', 'MorpoDisp_Scale ~ Diversification_ES_Scale + Time_Scale')}
        if(z==2)
        {path_equation$equation <-c('Richness_Scale ~ Diversification_ClaDS_Scale + Time_Scale + Beta_Scale + MorpoDisp_Scale', 'Beta_Scale ~ Diversification_ClaDS_Scale + Time_Scale', 'MorpoDisp_Scale ~ Diversification_ClaDS_Scale + Time_Scale')}

        ### Create models for each path and weighting scheme ###
        # We tested three schemes ('W', 'C', 'S'), two distance (d1, d2) and two models (OLS and SAR)
        # To find the best combination that minimize the autocorrelation and increase the fit for each path equation


        for (i in 1:length(scheme))

        {
            # Create a directory by scheme
            setwd(paste(dir_work, Clade_Folder[j], Div_Folder[z], sep='/'))
            dir_scheme <- scheme[i]
            if(!dir.exists(dir_scheme)) dir.create(dir_scheme)
            setwd(dir_scheme)
            # Create weighting scheme distance matrix
            spatial_weights_d1 <- nb2listw(d1, zero.policy=TRUE, style=scheme[i])
            spatial_weights_d2 <- nb2listw(d2, zero.policy=TRUE, style=scheme[i])
            
            for (p in 1:nrow(path_equation))
            {
                schemes <- rep(c(paste(scheme[i])), times=3)
                paths <- rep(c(paste(path_equation[p,1])), times=3)
                models <- c("lm_mod", "error_d1", "error_d2")
                results_path <- data.frame(scheme=schemes, path=paths, model=models)
                # OLS model 
                lm_mod <- lm(paste(path_equation[p,2]),
                        data = tablas[[j]])
                lm_mod_s <- summary(lm_mod)
                R2_lm <- lm_mod_s[["adj.r.squared"]]
                # SAR error models 
                # Minimum distance
                error_d1 <- spatialreg::errorsarlm(paste(path_equation[p,2]),
                data = tablas[[j]],listw = spatial_weights_d1,tol=1e-12,zero.policy=T)
                error_d1_s <-summary(error_d1, Nagelkerke=TRUE)
                R2_error_d1 <- error_d1_s$NK
                # Maximum distance
                error_d2 <- spatialreg::errorsarlm(paste(path_equation[p,2]),
                data = tablas[[j]],listw = spatial_weights_d2, tol=1e-12,zero.policy=T)
                error_d2_s <-summary(error_d2, Nagelkerke=TRUE)
                R2_error_d2 <- error_d2_s$NK
                # Save R2 (pseudo R2 for errorsar) and AIC
                results_path$R2_models <- c(R2_lm, R2_error_d1, R2_error_d2)
                results_path$AIC <- AIC(lm_mod, error_d1, error_d2)
                write.csv(results_path, file=paste(path_equation[p,1], ".csv"),
                row.names=F)
                # Make correlograms of residual autocorrelation
                cor.ols1.res<-correlog(tablas[[j]]$x, tablas[[j]]$y, z=residuals(lm_mod), na.rm=TRUE, increment=1, resamp=1)
                cor.sar1.res<-correlog(tablas[[j]]$x, tablas[[j]]$y, z=residuals(error_d1),na.rm=TRUE, increment=1, resamp=1)
                cor.sar2.res<-correlog(tablas[[j]]$x, tablas[[j]]$y, z=residuals(error_d2), na.rm=TRUE, increment=1, resamp=1)
                # Save correlograms in a jpeg file
                jpeg(filename = paste(path_equation[p,1],".jpg"),
                width = 215, height = 279, units = "mm",res = 600)
                par(mfrow = c(3, 1))
                plot(cor.ols1.res, xlab = "Distance (Km)", ylab = "Moran's I", ylim=c(-1,1), type = "l",
                lwd= 2, main = paste(models[1]), cex.main=2, cex.lab=1.8, cex.axis=1.5)
                abline(h=0, lty=5)
                plot(cor.sar1.res, xlab = "Distance (Km)", ylab = "Moran's I", ylim=c(-1,1), type = "l",
                lwd= 2, main =paste(models[2]), cex.main=2, cex.lab=1.8, cex.axis=1.5)
                abline(h=0, lty=5)
                plot(cor.sar2.res, xlab = "Distance (Km)", ylab = "Moran's I", ylim=c(-1,1), type = "l",
                lwd= 2, main = paste(models[3]), cex.main=2, cex.lab=1.8, cex.axis=1.5)
                abline(h=0, lty=5)
                dev.off()
                # Save path results
                save.image(file=paste(path_equation[p,1], ".RData"))
            }
        }
    


        ### Summarize the results in a data frame easily accessable ###
        results <- data.frame()
        for (i in 1:length(scheme)){
            setwd(paste(dir_work, Clade_Folder[j], Div_Folder[z], sep='/'))
        # Extract results from each scheme directory
            dir_scheme <- scheme[i]
            setwd(dir_scheme)
            scheme_results_list <- list.files(path = ".", pattern= '.csv$')
        # Join all results
            scheme_results <- lapply(scheme_results_list, read.csv)
            scheme_results2 <- do.call("rbind", scheme_results)
            results <- rbind(results, scheme_results2)
                        }
        # Back to your working directory
        setwd(paste(dir_work, Clade_Folder[j], Div_Folder[z], sep='/'))

        # Save the model selection results in a csv 
        write.csv(results, file="results_model_selection.csv", row.names=F)

        # Check results and select the better fitted scheme and model
        # by the lower AIC and the higher pseudo-R for each path

        ### Piecewise Structural Equation modelling evaluation ###
        # After the model selection we choose the "W" spatial weighting matrix with d2 distance (maximum distance) for all paths 

        # Create spatial weighting matrix

        results[,3]<-gsub('error_d1', 'd1', results$model)
        results[,3]<-gsub('error_d2', 'd2', results$model)

        rp1=subset(results, path=='path1')
        rp1=rp1[which(rp1$AIC.AIC==min(rp1$AIC.AIC)),]
        rp2=subset(results, path=='path2')
        rp2=rp2[which(rp2$AIC.AIC==min(rp2$AIC.AIC)),]
        rp3=subset(results, path=='path3')
        rp3=rp3[which(rp3$AIC.AIC==min(rp3$AIC.AIC)),]

        list_d=list(d1, d2)

        if (rp1$model=='d1') {spatial_Path1_weights <- nb2listw(list_d[[1]], zero.policy=TRUE, style=rp1$scheme)}
        if (rp1$model=='d2') {spatial_Path1_weights <- nb2listw(list_d[[2]], zero.policy=TRUE, style=rp1$scheme)}

        if (rp2$model=='d1') {spatial_Path2_weights <- nb2listw(list_d[[1]], zero.policy=TRUE, style=rp2$scheme)}
        if (rp2$model=='d2') {spatial_Path2_weights <- nb2listw(list_d[[2]], zero.policy=TRUE, style=rp2$scheme)}

        if (rp3$model=='d1') {spatial_Path3_weights <- nb2listw(list_d[[1]], zero.policy=TRUE, style=rp3$scheme)}
        if (rp3$model=='d2') {spatial_Path3_weights <- nb2listw(list_d[[1]], zero.policy=TRUE, style=rp3$scheme)}


        # Testing the theoretical pSEM model
        
        
        
        if(z==1) {sem_sar_model_teorico <- piecewiseSEM::psem(
        # 1 # Equation 1: species richness as response
        spatialreg::errorsarlm(Richness_Scale ~ Time_Scale + Diversification_ES_Scale + MorpoDisp_Scale + Beta_Scale, 
                                data = tablas[[j]],
                                listw = spatial_Path1_weights,
                                tol=1e-12,zero.policy=T),
                                
        spatialreg::errorsarlm(Beta_Scale ~ Time_Scale + Diversification_ES_Scale, 
                                data = tablas[[j]],
                                listw = spatial_Path2_weights,
                                tol=1e-12,zero.policy=T),
                                
        spatialreg::errorsarlm(MorpoDisp_Scale ~ Time_Scale + Diversification_ES_Scale, 
                                data = tablas[[j]],
                                listw = spatial_Path3_weights,
                                tol=1e-12,zero.policy=T), 
        data=tablas[[j]])}

        if(z==2) {sem_sar_model_teorico <- piecewiseSEM::psem(
        # 1 # Equation 1: species richness as response
        spatialreg::errorsarlm(Richness_Scale ~ Time_Scale + Diversification_ClaDS_Scale + MorpoDisp_Scale + Beta_Scale, 
                                data = tablas[[j]],
                                listw = spatial_Path1_weights,
                                tol=1e-12,zero.policy=T),
                                
        spatialreg::errorsarlm(Beta_Scale ~ Time_Scale + Diversification_ClaDS_Scale, 
                                data = tablas[[j]],
                                listw = spatial_Path2_weights,
                                tol=1e-12,zero.policy=T),
                                
        spatialreg::errorsarlm(MorpoDisp_Scale ~ Time_Scale + Diversification_ClaDS_Scale, 
                                data = tablas[[j]],
                                listw = spatial_Path3_weights,
                                tol=1e-12,zero.policy=T), 
        data=tablas[[j]])}


        # Run summary model 
        summary_sar_model_teorico <- summary(sem_sar_model_teorico)

        # Save RData with the pSEM model and summary model
        save(sem_sar_model_teorico, summary_sar_model_teorico, file="psem_sar_teorico_mollprj_50km_ahull.RData")

        # Extract path coefficients
        coefs_sar_model_teorico <- coefs(sem_sar_model_teorico)

        # Save coefficients
        write.csv(coefs_sar_model_teorico, file = "coefs_psem_sar_modelo_teorico_mollprj_50km.csv")

    }

}



########################################
# Correlation Lavaan vs piecewiseSEM   #
########################################

library(readxl)

table_correlation<-read_excel('~/Library/CloudStorage/OneDrive-CityUniversityofNewYork/CUNY PhD/Thesis/Cap 1 Macro Montanas/Results/Lavaan_piecewiseSEM_Correlation.xlsx',sheet=1)

par(mfrow=c(1,2))
plot(table_correlation$Lavaan_DR~table_correlation$piecewiseSEM_DR, xlab='piecewiseSEM', ylab='Lavaan', main='DR statistic', pch=19)

plot(table_correlation$Lavaan_ClaDS~table_correlation$piecewiseSEM_ClaDS, xlab='piecewiseSEM', ylab='Lavaan', main='ClaDS', pch=19)

#######################################################################
# Spatial models based on generalized least squares regression (GLSs) #
#######################################################################

library(nlme)
library(ape)
library(MuMIn)

autocor <- gls( Diversification ~ Elevation, data = tablas[[1]] )
semivario <- Variogram(autocor, form = ~x + y, resType = "normalized")
plot(semivario, smooth = TRUE)

exponential.autocor <- gls( Diversification ~ Elevation, correlation = corExp(form = ~x + y, nugget=T), data = tablas[[1]] )
gaussian.autocor <- gls( Diversification ~ Elevation, correlation = corGaus(form = ~x + y, nugget=T), data = tablas[[1]] )
spherical.autocor <- gls( Diversification ~ Elevation , correlation = corSpher(form = ~x + y, nugget=T), data = tablas[[1]]  )
linear.autocor <- gls( Diversification ~ Elevation, correlation = corLin(form = ~x + y, nugget=T), data = tablas[[1]]  )
ratio.autocor <- gls( Diversification ~ Elevation, correlation = corRatio(form = ~x + y, nugget=T), data = tablas[[1]]  )

model.sel(exponential.autocor, gaussian.autocor, spherical.autocor, linear.autocor, ratio.autocor)

plot(fitted(spherical.autocor), residuals(spherical.autocor))
abline(h=0,lty=3)

semivario <- Variogram(spherical.autocor, form = ~x + y, resType = "normalized")
plot(semivario, smooth = TRUE)

summary(spherical.autocor)


##################################################
### Generalized Lineal Model regression (GLMs) ###
##################################################

library(spdep)
library(ncf)
library(MASS)

AA=tablas[[1]][,1:2]

model1=glm((tablas[[1]]$Diversification-min(tablas[[1]]$Diversification))+1~tablas[[1]]$Elevation, family=gaussian)
AA[,3]<-resid(model1)
correlog1.1 <- correlog(AA[,1], AA[,2], AA[,3], na.rm=T, increment=1, resamp=100)

model2=glm((tablas[[1]]$Diversification-min(tablas[[1]]$Diversification))+1~tablas[[1]]$Elevation, family=Gamma(link = "log"))
AA[,4]<-resid(model2)
correlog2 <- correlog(AA[,1], AA[,2], AA[,4], na.rm=T, increment=1, resamp=100)

model7=glm((tablas[[1]]$Diversification-min(tablas[[1]]$Diversification))+1~tablas[[1]]$Elevation, family = gaussian(link = "log"))
AA[,9]<-residuals(model7)
correlog7 <- correlog(AA[,1], AA[,2], AA[,9], na.rm=T, increment=1, resamp=100)

model8=glm((tablas[[1]]$Diversification-min(tablas[[1]]$Diversification))+1~tablas[[1]]$Elevation, family = gaussian(link = "inverse"))
AA[,10]<-residuals(model8)
correlog8 <- correlog(AA[,1], AA[,2], AA[,10], na.rm=T, increment=1, resamp=100)



#################################################
###           GLMs with glmmfields            ###
#################################################

#m_spatial <- glmmfields(Div ~ Elevation,
#                        data = tabla1, family = Gamma(link = "log"),
#                        lat = "y", lon = "x", nknots = 12, iter = 500, chains = 1,
#                        prior_intercept = student_t(3, 0, 10), 
#                        prior_beta = student_t(3, 0, 3),
#                        prior_sigma = half_t(3, 0, 3),
#                        prior_gp_theta = half_t(3, 0, 10),
#                        prior_gp_sigma = half_t(3, 0, 3),
#                        seed = 123 # passed to rstan::sampling()
#)


#plot(m_spatial, type = "spatial-residual", link = TRUE) +
#  geom_point(size = 3)

##################################################
### Spatial autoregressive models (SAR) models ###
##################################################

library(spatialreg)
library(spdep)

sink("SARs.txt")

for (i in 1:3)
  
{
  
  # Define connectivity matrix (0/1)
  nbdist<-dnearneigh(x=as.matrix(tablas[[i]][,1:2]), d1=0, d2=1) # 
  
  # Compute the Euclidean distance between neighbouring sites
  neigh.dist<-nbdists(nbdist,
                      as.matrix(tablas[[i]][,1:2]), longlat=F)
  
  # Compute the inverse distance weigthed matrix
  inverse<-lapply(neigh.dist, function(x) (1/(x^2)))
  
  # Coding style W = row standardised
  nlw <- nb2listw(neighbours=nbdist, 
                  glist=inverse,
                  style="W", 
                  zero.policy=T) # Use zero policy to avoid error for any empty neighbour sets
  
  DivRate.sar <- spatialreg::errorsarlm(Diversification ~ Elevation,
                                        data = tablas[[i]],
                                        listw = nlw, 
                                        zero.policy=T) # Use zero policy to avoid error for any empty neighbour sets
  
  print(summary(DivRate.sar, Nagelkerke=TRUE))
  
  
  Time.sar <- spatialreg::errorsarlm(Time ~ Elevation,
                                     data = tablas[[i]],
                                     listw = nlw, 
                                     zero.policy=T) # Use zero policy to avoid error for any empty neighbour sets
  
  print(summary(Time.sar, Nagelkerke=TRUE))
  
  Morpho.sar <- spatialreg::errorsarlm(MorpoDisp ~ Elevation,
                                       data = tablas[[i]],
                                       listw = nlw, 
                                       zero.policy=T) # Use zero policy to avoid error for any empty neighbour sets
  
  print(summary(Morpho.sar, Nagelkerke=TRUE))
  
  Beta.sar <- spatialreg::errorsarlm(Beta_Diversity ~ Elevation,
                                     data = tablas[[i]],
                                     listw = nlw, 
                                     zero.policy=T) # Use zero policy to avoid error for any empty neighbour sets
  
  print(summary(Beta.sar, Nagelkerke=TRUE))
}

sink()



###############################################################
# Spatial regression with GAMs (generalized additive models) ##
###############################################################

library(mgcv)

Familias=c("###### Trochilidae ######", "###### Suboscines #######", "###### OSCINES ######")

sink("GAMSs.txt")

for (i in 1:3)
  
{
  
  print(Familias[i])
  
  elevation_gam <- gam(Diversification ~ s(Elevation, bs='cr'), data = tablas[[i]])
  print(summary(elevation_gam))
  elevation_gam <- gam(Time ~ s(Elevation, bs='cr'), data = tablas[[i]])
  print(summary(elevation_gam))
  elevation_gam <- gam(MorpoDisp ~ s(Elevation, bs='cr'), data = tablas[[i]])
  print(summary(elevation_gam))
  elevation_gam <- gam(Beta_Diversity ~ s(Elevation, bs='cr'), data = tablas[[i]])
  print(summary(elevation_gam))
  
}


sink()

#gam.check(elevation_gam)
#plot(elevation_gam)

anova()

##################
# Autocovariate Spatial regression
##################

#Fit the model with environmental variables

env_glm <- glm(Time ~ Elevation, tabla2, family = gaussian)

#RAC model (autocovariate derived from residuals of model with environmental predictors)

xy_residuals <- cbind(tabla2[,1:2], resid(env_glm))

rast=RiquezaTotal

rast[cellFromXY(rast, xy_residuals)] <- xy_residuals[,3]

plot(rast)

#Calculate residuals autocovariate
#Focal operations: ngb is neighbourhood size, set to 3 by 3 cells; fun is function, #here the mean value within the defined neighbourhood

focal_rac_rast <- raster::focal(rast, w=matrix(1/9, nc=3, nr=3), fun = 'mean', na.rm = TRUE)
plot(focal_rac_rast)

#Extract the values of the focal operation from “focal_rac_rast” rasterfile using the #co‐ordinates stored in “xy”
focal_rac_vect <- extract(focal_rac_rast, coordinates(tabla2[,1:2]))

tabla2<- cbind(tabla2, focal_rac_vect)

rac_glm <- glm(Time ~ Elevation +focal_rac_vect, tabla2, family = gaussian)
summary(rac_glm)

#calculate McFadden's R-squared for model
with(summary(rac_glm), 1 - deviance/null.deviance)


