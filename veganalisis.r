#R
#2022
#veg analysis peat communities

std <- function(x) sd(x)/sqrt(length(x))

library(vegan)
library(reshape)
library(ggplot2)

library(ade4)
library(adespatial)
library(vegan)
library(gclus)
library(cluster)
library(FD)
library(agricolae)
library(RColorBrewer)

library(labdsv)


dirmac<-"/Users/juanbenavides/Dropbox/papers/peatland comunities/stats3"
setwd(dirmac)

source("boxplerk.R")
source("boxplert.R")
source("drawmap.R")
source("coldiss.R")
source("hcoplot.R")




veg.all<-read.csv("vegetation3.csv",header=TRUE)
env<-read.csv("cations2.csv",header=TRUE)


veg.all$logcover<-log(veg.all$COVER+1.1)
############
#Melt and cast operations for all species and plot type

veg.all.1<-melt(veg.all,id.vars=c(2,3),measure.vars=4)


veg.all.cast<-data.frame(cast(veg.all.1,PLOT~SPECIES,sum))

#nmds of all samples and species including rare species (abundance =1)
all.nmds<-metaMDS(veg.all.cast[,-1],k=2,trymax=100)

#morethanone.nmds<-metaMDS(veg.all.morethan1.dist,k=2,trymax=100)

plot(all.nmds)

######
#removing species that are only once
veg.all.counts<-data.frame(cast(veg.all.1,SPECIES~PLOT,sum))
veg.all.counts.1<-veg.all.counts[,2:dim(veg.all.counts)[2]]
veg.all.counts2<-data.frame(ifelse(veg.all.counts.1>0,1,FALSE))
count<-data.frame(species=veg.all.counts[,1],frequency=rowSums(veg.all.counts2))
morethanone<-subset(count,frequency>1)
veg.morethanone<-merge(morethanone,veg.all,by.x='species',by.y='SPECIES',all.y=FALSE)
  ############
#Melt and cast operations for more than one species and all plot types
morethanone.1<-melt(veg.morethanone,id.vars=c(1,2,4),measure.vars=5)
morethanone.2<-cast(morethanone.1,PLOT~species,sum)
###########################################################################
#######    ##################################
#######NMDS##################################
#######    ##################################
#### NMDS of ground layer species including singletons and all comm. types
#nmds of all samples and species including rare species (abundance =1)

morethanone.nmds<-metaMDS(morethanone.2[c(2:188)],k=2,trymax=100, autotransform=FALSE)


#species list for centroids

#from veganalisis

spec2<-read.csv("ind.spe.table3.csv",header=TRUE)
list1<-data.frame(Species=colnames(morethanone.2[,-1]))
spec3<-merge(spec2,list1,by="Species")
spec4<-spec3[,10]


#plot nmds and spec

plot(morethanone.nmds,type="n")
points(morethanone.nmds, display =  "species", pch=20,cex=0.4)
text(morethanone.nmds,  display ="species", spec4,cex=0.4,col="blue")

plot(morethanone.nmds,type="n")
points(morethanone.nmds, display =  "sites", pch=20,cex=0.4)
text(morethanone.nmds,  display ="sites", label= morethanone.2[,1], cex=0.4,col="blue")


env.vectors<-vectorfit(morethanone.nmds$points, env[,c(6:8)], permutations=999)
plot(env.vectors,p.max = 0.05, col = "red")

########
#heat map dissimilarities
spe.ds2 <- vegdist(morethanone.2[,-1], method = "bray", binary = TRUE)
coldiss(spe.ds2, byrank = FALSE, diag = TRUE)
#coldiss is a script from numerical ecology book coldiss.R http://www.numericalecology.com/numecolR/

spe.KM.cascade <- cascadeKM(
morethanone.2[,-1],
inf.gr = 2,
sup.gr = 10,
iter = 100,
criterion = "calinski")

plot(spe.KM.cascade)
summary(spe.KM.cascade, sortg = TRUE)

veg.bray <- vegdist(morethanone.2[,-1], method = "bray", binary = TRUE)

veg.ward<-hclust(d = veg.bray, method = "ward.D")
#order labels with plot names

rownames(labels_ord)<- morethanone.2[,1]
veg.ward.g <- cutree(veg.ward, k = 5)


plotnames<-data.frame(id=1:121,plot=morethanone.2[,1])

plotnames<-plotnames[veg.ward$order,]

par(cex=0.4)
plot(veg.ward,labels=plotnames$plot,xlab="", ylab="", main="", sub="", axes=FALSE)
rect.hclust(veg.ward, k = 5)

par(mfrow=c(1,2),cex=0.4)
hcoplot(veg.ward, veg.bray, lab = morethanone.2[,1], k =5)

labels_ord<-data.frame(morethanone.nmds$points)

rownames(labels_ord)<- morethanone.2[,1]
veg.ward.g <- cutree(veg.ward, k = 5)


#env variables among groups
groups<-data.frame(plot=rownames(labels_ord),clusters=veg.ward.g)
groups$names<-ifelse(groups$clusters==1,"Sphagnum",NA)
groups$names<-ifelse(groups$clusters==2,"True mosses",groups$names)
groups$names<-ifelse(groups$clusters==3,"Sedges",groups$names)
groups$names<-ifelse(groups$clusters==4,"Grasses",groups$names)
groups$names<-ifelse(groups$clusters==5,"Cushion",groups$names)
groups$names<-as.factor(groups$names)

#rownames(labels_ord)<- morethanone.2[,1]
veg.ward.g <- cutree(veg.ward, k = 5)
drawmap(labels_ord,clusters=veg.ward.g)

#####################################################################################
#main gradients not using all peat chemistry or carbon
env.vectors<-vectorfit(morethanone.nmds$points, env[,6:8], permutations=999)

data.env<-data.frame(variables= rownames(env.vectors$arrows), arrowsX=env.vectors$arrows[,1],arrowsY=env.vectors$arrows[,2],X0=c(0,0,0),Y0=c(0,0,0))

arrows(data.env$X0,data.env$Y0,data.env$arrowsX,data.env$arrowsY)
text(data.env$arrowsX,data.env$arrowsY-0.1,labels=data.env$variables)

#text(labels_ord,labels=rownames(labels_ord))



spe.ordered <- reorder.hclust(veg.ward, veg.bray)

spe.order.syn<-vegemite(morethanone.2[,-1] , spe.ordered,scale="Hill")

dend <- as.dendrogram(veg.ward)


heatmap(
# Heat map of the doubly ordered community table, with dendrogram
t(morethanone.2[,-1][(spe.order.syn$species)]),
    Rowv = NA, Colv = dend, col = c("white", brewer.pal(5, "Greens")),
    scale = "none", margin = c(4, 4), labRow = spe.order.syn$table[-c(1:3),1],
    ylab = NA, xlab = "Sites")


#indicator species analysis


peat_ind<-indval(morethanone.2[,-1],veg.ward.g,numitr=999)
sink("peat_ind.txt")
summary(peat_ind)
sink()


############################################################
############################################################
#load environmental data with veg groups

env1<-merge(env,groups,by="plot")

veg.exp<-merge(veg.all,groups,by.x="PLOT",by.y="plot")

species<-read.csv("species cat.csv",header=TRUE)
veg.all2<-merge(veg.exp,species,by="SPECIES")

write.csv(veg.all2,"vegall2.csv")



##################################################################
##################################################################
#permanova analysis

permanova<-adonis2(morethanone.2[,-1]~names, data=env1)

#load paired test for permanova adonis2 
#https://rdrr.io/github/GuillemSalazar/EcolUtils/man/adonis.pair.html

braydist<-vegdist(morethanone.2[,-1], method="bray") 

dist.mat<-vegdist(morethanone.2[,-1], method="bray")
  Factor<-as.factor(env1$names)
  nper<-9999
  comb.fact<-combn(levels(Factor),2)
  corr.method<-"fdr"
  pv<-NULL
  R2<-NULL
  SS<-NULL
  MeanSqs<-NULL
  F.Model<-NULL

  for (i in 1:dim(comb.fact)[2]){
    model.temp<-adonis2(as.dist(as.matrix(dist.mat)[Factor==comb.fact[1,i] | Factor==comb.fact[2,i],Factor==comb.fact[1,i] | Factor==comb.fact[2,i]])~Factor[Factor==comb.fact[1,i] | Factor==comb.fact[2,i]],permutations=nper)
    pv<-c(pv,model.temp[[5]][1])
    R2<-c(R2,model.temp[[3]][1])
    SS<-c(SS,model.temp[[2]][1])
    MeanSqs<-c(MeanSqs,model.temp[[2]][1]/model.temp[[1]][1])
    F.Model<-c(F.Model,model.temp[[4]][1])
    print(i)
    print(model.temp[[3]][1])
    }

  pv.corr<-p.adjust(pv,method=corr.method)
  permanova.table<-data.frame(combination=paste(comb.fact[1,],comb.fact[2,],sep=" <-> "),SumsOfSqs=SS,MeanSqs=MeanSqs,F.Model=F.Model,R2=R2,P.value=pv,P.value.corrected=pv.corr)
write.csv(permanova.table,"permanova.csv")




#write.csv(veg.exp,"vegetation3.1.csv")


#######
#water content field samples at field capacity

env1$field.cap1<-(env1$Peat1.w.weight-env1$dry.weight.1)/env1$Peat1.w.weight*100
env1$field.cap2<-(env1$Peat2.w.weight-env1$dry.weight.2)/env1$Peat2.w.weight*100

env1$Carbon1.g.m2<- env1$Carbon1.g.g*env1$BD1*1E5
env1$Carbon2.g.m2<- env1$Carbon2.g.g*env1$BD2*1E5


env1$names<-factor(env1$names,levels=c("Cushion", "Sedges","Grasses", "True mosses","Sphagnum"))

             

########
#boxplots for water chemistry
k<-5
par(mfrow=c(3,2),mar=c(2,5,2,1))
# Use boxplert() or boxplerk() to plot results with post-hoc tests
with(env1, {
boxplert(K, names, xlab = "", ylab = expression("K"^"+" * "(mg l"^-1 *")"), bcol = (1:k) + 1, p.adj = "holm",ylim=c(0,10) )
boxplerk( (Mg), names, xlab = "", ylab = expression("Na"^"+" * "(mg l"^-1 *")"), bcol = (1:k) + 1, p.adj = "holm",ylim=c(0,10) )
boxplert( (Ca), names, xlab = "", ylab = expression("Ca"^"++" *" (mg l"^-1 *")"), bcol = (1:k) + 1, p.adj = "holm",ylim=c(0,40) )
boxplerk( (Conductivity), names, xlab = "", ylab = expression("Conductivity ("*mu* "S cm"^-1*")"), bcol = (1:k) + 1, p.adj = "holm",ylim=c(0,400) )
boxplert(pH, names, xlab = "", ylab = "pH", bcol = (1:k) + 1, p.adj = "holm" ,ylim=c(2.5,8) )
boxplert( Elevation,  names, xlab = "", ylab = "Elevation (m)", bcol = (1:k) + 1, p.adj = "holm",ylim=c(2500,4800) )
})

########
#boxplots for peat carbon
par(mfrow=c(3,2),mar=c(2,5,2,1))
with(env1, {
boxplert( Carbon1.g.g, names, xlab = "", ylab = expression("Carbon  0-10 cm (gC g"^-1*")"), bcol = (1:k) + 1, p.adj = "holm",ylim=c(0,0.5))
boxplerk( Carbon2.g.g, names, xlab = "", ylab = expression("Carbon  10-20 cm (gC g"^-1*")"), bcol = (1:k) + 1, p.adj = "holm",ylim=c(0,0.5))
boxplerk((BD1), names, xlab = "", ylab = expression(" Bulk density 0-10 cm (g cm"^-3*")"), bcol = (1:k) + 1, p.adj = "holm" , ylim=c(0,0.2))
boxplert((BD2),  names, xlab = "", ylab = expression(" Bulk density 10-20 cm (g cm"^-3*")"), bcol = (1:k) + 1, p.adj = "holm", ylim=c(0,0.2))
boxplert(Carbon1.g.m2, names, xlab = "", ylab = expression("Carbon content 0-10 cm (gC m"^-2*")"), bcol = (1:k) + 1, p.adj = "holm" , ylim=c(500,5000))
boxplert( (Carbon2.g.m2), names, xlab = "", ylab = expression("Carbon content 10-20 cm (gC m"^-2*")"), bcol = (1:k) + 1, p.adj = "holm", ylim=c(500,5000))
})


########
#boxplots for peat water and chemistry
par(mfrow=c(3,2),mar=c(2,5,2,1))
with(env1, {
boxplerk(field.cap1, names, xlab = "", ylab = "Water content (0-10cm) (%)", bcol = (1:k) + 1, p.adj = "holm" , ylim=c(80,100))
boxplert( field.cap2, names, xlab = "", ylab = "Water content (10-20cm) (%)", bcol = (1:k) + 1, p.adj = "holm", ylim=c(80,100))
boxplert( Ca.g.kg, names, xlab = "", ylab = expression("Peat Ca (mg g"^-1*")"), bcol = (1:k) + 1, p.adj = "holm", ylim=c(0,45))
boxplert( K.g.kg, names, xlab = "", ylab = expression("Peat K (mg g"^-1*")"), bcol = (1:k) + 1, p.adj = "holm", ylim=c(0,8))
boxplert( (Mg.g.kg), names, xlab = "", ylab = expression(" Peat Mg (mg g"^-1*")"), bcol = (1:k) + 1, p.adj = "holm" , ylim=c(0,3))
boxplert( (Na.g.kg),  names, xlab = "", ylab = expression("Peat Na (mg g"^-1*")"), bcol = (1:k) + 1, p.adj = "holm", ylim=c(0,3))
})


######
#anovas for all env variables

#########
#water chemistry

W.Ca<-aov((Ca)~names,data=env1) ; summary(W.Ca)
W.K<-aov(K~names,data=env1)   ; summary(W.K)
W.Mg<-kruskal.test(Mg~names,data=env1) ; (W.Mg)
W.Na<-aov((Na)~names,data=env1) ; summary(W.Na)
W.pH<-aov(pH~names,data=env1) ; summary(W.pH)
W.Conductivity<-kruskal.test((Conductivity)~names,data=env1)   ; (W.Conductivity)

#########
#peat carbon

P.C1gg<-aov(Carbon1.g.g~names,data=env1) ; summary(P.C1gg)
P.C2gg<-kruskal.test(Carbon2.g.g~names,data=env1) ; (P.C2gg)
P.BD1<-kruskal.test((BD1)~names,data=env1) ; (P.BD1)
P.BD2<-aov(BD2~names,data=env1) ; summary(P.BD2)
P.CC1<-aov(Carbon1.g.m2~names,data=env1) ; summary(P.CC1)
P.CC2<-aov(Carbon2.g.m2~names,data=env1) ; summary(P.CC2)

#########
#peat chemistry and water

P.FC1<-aov(field.cap1~names,data=env1) ; summary(P.FC1)
P.FC2<-aov(field.cap1~names,data=env1) ; summary(P.FC2)
P.Ca<-aov(Ca.g.kg~names,data=env1) ; summary(P.Ca)
P.K<-aov(K.g.kg~names,data=env1) ; summary(P.K)
P.Mg<-aov(Mg.g.kg~names,data=env1) ; summary(P.Mg)
P.Na<-aov(Na.g.kg~names,data=env1) ; summary(P.Na)

########################################################################
########################################################################
#disturbance

dist<-read.csv("disturbance2.csv",header=TRUE)

dist.DF<-data.frame(plot=dist[,1],Towndist1=scale(dist$Town.distance.km)*-1,
                  Housedist1=scale(dist$House.distance)*-1,
                  Cattledistance1=scale(dist$Cattle.distance)*-1,
                  Catleinten1=scale(dist$Cattle.intensity.indiv),
                  Ag.dist1=scale(dist$Ag.field.distance)*-1,
                  Ag.area1=scale(dist$Ag.field.area.ha),
                  Ag.runoff1=scale(dist$Ag.runoff.ha),
                  Townsize1=scale(log(dist$Town.size)),
                  Towndist=(dist$Town.distance.km),
                  Housedist=(dist$House.distance),
                  Cattledistance=(dist$Cattle.distance),
                  Catleinten=(dist$Cattle.intensity.indiv),
                  Ag.dist=(dist$Ag.field.distance),
                  Ag.area=(dist$Ag.field.area.ha),
                  Ag.runoff=(dist$Ag.runoff.ha),
                  Townsize=((dist$Town.size)))


dist.DF$distindex<-rowSums(dist.DF[,2:9])
dist.DF$distindex2<-dist.DF[,2]+dist.DF[,3]/3+dist.DF[,9]/3


env2<-merge(env1,dist.DF,by="plot")



########
#boxplots for disturbance
par(mfrow=c(3,2),mar=c(2,5,2,1))
with(env2, {
boxplert(distindex, names, xlab = "", ylab = "distindex", bcol = (1:k) + 1, p.adj = "holm")
boxplert( distindex2, names, xlab = "", ylab = "distindex2", bcol = (1:k) + 1, p.adj = "holm")
boxplert(Ag.dist, names, xlab = "", ylab = "Ag dist", bcol = (1:k) + 1, p.adj = "holm" )
boxplert(Ag.area,  names, xlab = "", ylab = "Agarea", bcol = (1:k) + 1, p.adj = "holm")
boxplert(Towndist, names, xlab = "", ylab = "Town dist", bcol = (1:k) + 1, p.adj = "holm" )
boxplerk(log(Townsize), names, xlab = "", ylab = "Town size" , bcol = (1:k) + 1, p.adj = "holm")
})


modeldist<-lm(Carbon1.g.g~ distindex2 * names, data=env2)
modeldist2<-lm(Carbon1.g.g~ distindex2 + names, data=env2)
summary(modeldist)
summary(modeldist2)
anova(modeldist)

modeldist3<-lm(distindex2~   names, data=env2)
anova(modeldist3)

plot(env2$distindex2,env2$Carbon1.g.m2,col=as.factor(env2$names))
text(env2$distindex2,env2$Carbon1.g.m2,label=env2[,1])


dist.plot<-ggplot(data=env2,aes(x=distindex2,y=Carbon1.g.m2))+theme_bw()

dist.plot.1<-dist.plot+geom_point(aes(colour=names,shape=names))+
            stat_smooth(method = "lm",se = FALSE,colour="black")+
            scale_x_continuous("Disturbance index")+
            scale_y_continuous(expression("Organic Soil Carbon  0-10 cm (gC m"^-2*")"))+
            scale_colour_discrete("")+
            scale_shape_discrete("")+
            theme(legend.position=c(0.85,0.82),panel.grid.minor=element_blank  (),
                     panel.grid.major=element_blank(), panel.background=element_blank(),
                     axis.text.x =element_text(size = 12,hjust=0.5,angle=0,colour='black'),
                     axis.text.y =element_text(size = 12,angle=0,colour='black'),
                     axis.title.x = element_text(size = 12,colour='black'),
                     axis.title.y =element_text(size = 12, angle = 90,colour='black'),
                     axis.line.x = element_line(color="black", size = 0.5),
               axis.line.y = element_line(color="black", size = 0.5 ),
               legend.background=element_rect(fill = alpha("white", 0)),
                legend.key=element_rect(fill = alpha("white", .5)))

dist.plot.1


specb_elev.plot<-ggplot(data=env2,aes(x=Elevation,y=Bryo.species))+theme_bw()

specb_elev.plot.1<-specb_elev.plot+geom_point(aes(colour=names,shape=names))+
            scale_x_continuous("Elevation (m)")+
            scale_y_continuous("Bryophyte species")+
            scale_colour_discrete("")+
            scale_shape_discrete("")+
            theme(legend.position=c(0.16,0.7),panel.grid.minor=element_blank  (),
                     panel.grid.major=element_blank(), panel.background=element_blank(),
                     legend.key.width=(unit(0.5, 'cm')),
                     legend.title = element_text(size=8), #change legend title font size
                    legend.text = element_text(size=8), #change legend text font size
                     axis.text.x =element_blank(),
                     axis.text.y =element_text(size = 12,angle=0,colour='black'),
                     axis.title.x = element_blank(),
                     axis.title.y =element_text(size = 12, angle = 90,colour='black'),
                     axis.line.x = element_line(color="black", size = 0.5),
               axis.line.y = element_line(color="black", size = 0.5 ))

specb_elev.plot.1


specv_elev.plot<-ggplot(data=env2,aes(x=Elevation,y=Vasc.species))+theme_bw()

specv_elev.plot.1<-specv_elev.plot+geom_point(aes(colour=names,shape=names))+
            scale_x_continuous("Elevation (m)")+
            scale_y_continuous("Vascular species")+
            scale_colour_discrete("")+
            theme(legend.position="none",panel.grid.minor=element_blank  (),
                     panel.grid.major=element_blank(), panel.background=element_blank(),
                     axis.text.x =element_blank(),
                     axis.text.y =element_text(size = 12,angle=0,colour='black'),
                     axis.title.x = element_blank(),
                     axis.title.y =element_text(size = 12, angle = 90,colour='black'),
                     axis.line.x = element_line(color="black", size = 0.5),
               axis.line.y = element_line(color="black", size = 0.5 ))

specv_elev.plot.1

env2$Elevation2<-env2$Elevation^2

quadraticModel <- lm(Species ~ Elevation + Elevation2, data=env2)

#view model summary
summary(quadraticModel)

spec_elev.plot<-ggplot(data=env2,aes(x=Elevation,y=Species))+theme_bw()

spec_elev.plot.1<-spec_elev.plot+geom_point(aes(colour=names,shape=names))+
            stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1,se=FALSE,colour="black") +            
            scale_x_continuous("Elevation (m)")+
            scale_y_continuous("Plant species")+
            scale_colour_discrete("")+
            theme(legend.position="none",panel.grid.minor=element_blank  (),
                     panel.grid.major=element_blank(), panel.background=element_blank(),
                     axis.text.x =element_text(size = 12,hjust=0.5,angle=0,colour='black'),
                     axis.text.y =element_text(size = 12,angle=0,colour='black'),
                     axis.title.x = element_text(size = 12,colour='black'),
                     axis.title.y =element_text(size = 12, angle = 90,colour='black'),
                     axis.line.x = element_line(color="black", size = 0.5),
               axis.line.y = element_line(color="black", size = 0.5 ))

spec_elev.plot.1

grid.arrange(nrow=3,specv_elev.plot.1,specb_elev.plot.1,spec_elev.plot.1)



##################################################################
##################################################################
#export tables
#env table

table1<-aggregate(cbind(K,Na,Ca,Mg,pH,Conductivity,Elevation ,
Carbon1.g.g, Carbon2.g.g, BD1,BD2,Carbon1.g.m2, Carbon2.g.m2,
Mg.g.kg,K.g.kg,Ca.g.kg,Na.g.kg,wc1,wc2,Hum.holl.average, hum.holl.max, WT.depth, 
Sphagnum,Species,Bryo.species, Vasc.species,Towndist, 
Housedist, Cattledistance, Catleinten, Ag.dist,Ag.area,
Ag.runoff, Townsize)  ~ names,data=env2, FUN=mean)

table1.1<-round(table1[,-1],2)
table1.2<-cbind(table1[,1],table1.1)

tablestd<-aggregate(cbind(K,Na,Ca,Mg,pH,Conductivity,Elevation ,
Carbon1.g.g, Carbon2.g.g, BD1,BD2,Carbon1.g.m2, Carbon2.g.m2,
Mg.g.kg,K.g.kg,Ca.g.kg,Na.g.kg,wc1,wc2,Hum.holl.average, hum.holl.max, WT.depth, 
Sphagnum,Species,Bryo.species, Vasc.species,Towndist, 
Housedist, Cattledistance, Catleinten, Ag.dist,Ag.area,
Ag.runoff, Townsize)  ~ names,data=env2, FUN=std)

tablestd.1<-round(tablestd[,-1],3)
tablestd.2<-cbind(tablestd[,1],tablestd.1)

write.csv(table1.2,"table12.csv")
write.csv(tablestd.2,"tablestd2.csv")

##################################################################
#export tables
#ind species

ind.spe2<-read.csv("ind_sp_out.csv",header=TRUE)
all.spe2<-read.csv("all.spe.names.csv",header=TRUE)

ind.spe.table2<-merge(all.spe2,ind.spe2,by="Species",all=TRUE)

write.csv(ind.spe.table2,"ind.spe.table2.csv")

