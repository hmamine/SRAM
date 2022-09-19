#require packages for the analysis the analysis
pkg=c("plyr","ggplot2", "data.table", "reshape2","PMCMRplus","venn","ape","UpSetR")


lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])


list.color<-c("#4f7b95","#ffa500","#155832")

theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.title.x=element_blank(),
	axis.title.y=element_blank()
	)
theme_change <- theme(
	legend.position="none",
	plot.background = element_blank(),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x = element_blank(),	
	)

Man1=read.table("manresa_network_sub.txt", sep="\t", row.names=1, header=1)
Man1=data.table(Man1, keep.rownames=T, key="rn")
Man<-unique(c(Man1$from, Man1$to))
Mar1=read.table("marquis_network_sub.txt", sep="\t", row.names=1, header=1)
Mar1=data.table(Mar1, keep.rownames=T, key="rn")
Mar<-unique(c(Mar1$from, Mar1$to))
Swe1=read.table("sweetann_network_sub.txt", sep="\t", row.names=1, header=1)
Swe1=data.table(Swe1, keep.rownames=T, key="rn")
Swe<-unique(c(Swe1$from, Swe1$to))


man.full=read.table("manresa_metrics_full.txt", sep="\t", row.names=1, header=1)
man.full=data.table(man.full, key="rn")
Man.taxa=man.full[man.full$rn %in% Man]
Man.taxa$Cultivar<-"Manresa"
man.full$Sub<-ifelse(man.full$rn %in% Man, man.full$Phylum,"NotSub")



mar.full=read.table("marquis_metrics_full.txt", sep="\t", row.names=1, header=1)
mar.full=data.table(mar.full, key="rn")
Mar.taxa=mar.full[mar.full$rn %in% Mar]
Mar.taxa$Cultivar<-"Marquis"
mar.full$Sub<-ifelse(mar.full$rn %in% Mar, mar.full$Phylum,"NotSub")

swe.full=read.table("sweetann_metrics_full.txt", sep="\t", row.names=1, header=1)
swe.full=data.table(swe.full, key="rn")
Swe.taxa=swe.full[swe.full$rn %in% Swe]
Swe.taxa$Cultivar<-"SweetAnn"
swe.full$Sub<-ifelse(swe.full$rn %in% Swe, swe.full$Phylum,"NotSub")


All.taxa=rbind(Man.taxa,Mar.taxa,Swe.taxa)


Man2=read.table("manresa_metrics_sub.txt", sep="\t", row.names=1, header=1)
Man2=data.table(Man2, key="rn")
Man2$Col<-"manresa"

Mar2=read.table("marquis_metrics_sub.txt", sep="\t", row.names=1, header=1)
Mar2=data.table(Mar2,key="rn")
Mar2$Col<-"marquis"

Swe2=read.table("sweetann_metrics_sub.txt", sep="\t", row.names=1, header=1)
Swe2=data.table(Swe2, key="rn")
Swe2$Col<-"sweetann"

ManPos=length(Man1[Man1$color == "darkgreen"]$color)
ManNeg=length(Man1[Man1$color != "darkgreen"]$color)
#ManTot=ManPos+ManNeg
MarPos=length(Mar1[Mar1$color == "darkgreen"]$color)
MarNeg=length(Mar1[Mar1$color != "darkgreen"]$color)
#MarTot=MarPos+MarNeg

SwePos=length(Swe1[Swe1$color == "darkgreen"]$color)
SweNeg=length(Swe1[Swe1$color != "darkgreen"]$color)
#SweTot=SwePos+SweNeg

Pos=rbind(ManPos,MarPos,SwePos)
Neg=rbind(ManNeg,MarNeg,SweNeg)
#Tot=rbind(ManTot,MarTot,SweTot)
tab1=rbind(Pos,Neg)
colnames(tab1)<-"value"
tab1<-data.frame(tab1)
tab1$cultivar<-c("Manresa","Marquis","SweetAnn","Manresa","Marquis","SweetAnn")
tab1$type<-c("Pos","Pos","Pos","Neg","Neg","Neg")

	pp1=ggplot(tab1,aes( x=cultivar,y=value,fill=type))
	p1=pp1+geom_bar(position="stack", stat="identity")+
	scale_fill_manual(values=c("darkred","darkgreen"))+
	theme_bw()+theme_change

tab2=rbind(Man2,Mar2,Swe2)


Man3=read.table("manresa_sub_relative.txt", sep="\t", row.names=1, header=1)
Man3=data.table(Man3, key=c("Var1","Var2"))
data.Man3=dcast(setDT(Man3), Var1 ~ Var2 , value.var= "value")
rownames(data.Man3)<-data.Man3$Var1
data.Man3=data.Man3[,-1]
dataMan3=data.table(t(data.Man3),keep.rownames=TRUE)

ManPos<-c(unique(Man3[Man3$color=="darkgreen",]$Var1))
ManNeg<-c(unique(Man3[Man3$color!="darkgreen",]$Var1))

dataMan3$Pos.mean<-dataMan3[,.(rowSums(.SD,na.rm=TRUE)),.SDcols = ManPos]
dataMan3$Neg.mean<-dataMan3[,.(rowSums(.SD,na.rm=TRUE)),.SDcols = ManNeg]
dataMan3=dataMan3[,c(1,29,30)]
meltMan3=melt(dataMan3)
meltMan3$cultivar<-"Manresa"


Mar3=read.table("marquis_sub_relative.txt", sep="\t", row.names=1, header=1)
Mar3=data.table(Mar3, key=c("Var1","Var2"))
data.Mar3=dcast(setDT(Mar3), Var1 ~ Var2 , value.var= "value")
rownames(data.Mar3)<-data.Mar3$Var1
data.Mar3=data.Mar3[,-1]
dataMar3=data.table(t(data.Mar3),keep.rownames=TRUE)

MarPos<-c(unique(Mar3[Mar3$color=="darkgreen",]$Var1))
MarNeg<-c(unique(Mar3[Mar3$color!="darkgreen",]$Var1))

dataMar3$Pos.mean<-dataMar3[,.(rowSums(.SD,na.rm=TRUE)),.SDcols = MarPos]
dataMar3$Neg.mean<-dataMar3[,.(rowSums(.SD,na.rm=TRUE)),.SDcols = MarNeg]
dataMar3=dataMar3[,c(1,27,28)]
meltMar3=melt(dataMar3)
meltMar3$cultivar<-"Marquis"


Swe3=read.table("sweetann_sub_relative.txt", sep="\t", row.names=1, header=1)
Swe3=data.table(Swe3, key=c("Var1","Var2"))
data.Swe3=dcast(setDT(Swe3), Var1 ~ Var2 , value.var= "value")
rownames(data.Swe3)<-data.Swe3$Var1
data.Swe3=data.Swe3[,-1]
dataSwe3=data.table(t(data.Swe3),keep.rownames=TRUE)

SwePos<-c(unique(Swe3[Swe3$color=="darkgreen",]$Var1))
SweNeg<-c(unique(Swe3[Swe3$color!="darkgreen",]$Var1))

dataSwe3$Pos.mean<-dataSwe3[,.(rowSums(.SD,na.rm=TRUE)),.SDcols = SwePos]
dataSwe3$Neg.mean<-dataSwe3[,.(rowSums(.SD,na.rm=TRUE)),.SDcols = SweNeg]
dataSwe3=dataSwe3[,c(1,7,8)]
meltSwe3=melt(dataSwe3)
meltSwe3$cultivar<-"SweetAnn"


meltAll=rbind(meltMan3,meltMar3,meltSwe3)
meltAll$concat=paste0(meltAll$variable,".",meltAll$cultivar)
Order1<-c("Pos.mean.Manresa","Neg.mean.Manresa","Pos.mean.Marquis","Neg.mean.Marquis",
"Pos.mean.SweetAnn","Neg.mean.SweetAnn")


	Man.count=read.table("Species_table_manresa.txt", sep="\t", row.names=1, header=1)
	i=colnames(Man.count)
	Man.count$Mean=rowMeans(Man.count)
	Man.count=data.table(Man.count, keep.rownames=TRUE, key="rn")
	Man.count$mean<-Man.count[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = i]
	for (col in i) set (Man.count, j=col, value = +(Man.count[[col]]>0))
	Man.count$freq<-Man.count[,.(rowSums(.SD,na.rm=TRUE)),.SDcols = i]
	Man.count$freq<-Man.count$freq/length(i)
	Man.count=Man.count[,c("rn","mean","freq")]
	man.full=merge(man.full,Man.count)
	man.full=man.full[order(freq,decreasing=F),]
	man.full$lab=ifelse(man.full$rn %in% ManNeg ,man.full$rn, "")
	
	

	Mar.count=read.table("Species_table_marquis.txt", sep="\t", row.names=1, header=1)
	i=colnames(Mar.count)
	Mar.count$Mean=rowMeans(Mar.count)
	Mar.count=data.table(Mar.count, keep.rownames=TRUE, key="rn")
	Mar.count$mean<-Mar.count[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = i]
	for (col in i) set (Mar.count, j=col, value = +(Mar.count[[col]]>0))
	Mar.count$freq<-Mar.count[,.(rowSums(.SD,na.rm=TRUE)),.SDcols = i]
	Mar.count$freq<-Mar.count$freq/length(i)
	Mar.count=Mar.count[,c("rn","mean","freq")]
	mar.full=merge(mar.full,Mar.count)
	mar.full=mar.full[order(freq,decreasing=F),]
	mar.full$lab=ifelse(mar.full$rn %in% MarNeg ,mar.full$rn, "")
	
	
	
	Swe.count=read.table("Species_table_sweetann.txt", sep="\t", row.names=1, header=1)
	i=colnames(Swe.count)
	Swe.count$Mean=rowMeans(Swe.count)
	Swe.count=data.table(Swe.count, keep.rownames=TRUE, key="rn")
	Swe.count$mean<-Swe.count[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = i]
	for (col in i) set (Swe.count, j=col, value = +(Swe.count[[col]]>0))
	Swe.count$freq<-Swe.count[,.(rowSums(.SD,na.rm=TRUE)),.SDcols = i]
	Swe.count$freq<-Swe.count$freq/length(i)
	Swe.count=Swe.count[,c("rn","mean","freq")]
	swe.full=merge(swe.full,Swe.count)
	swe.full=swe.full[order(freq,decreasing=F),]
	swe.full$lab=ifelse(swe.full$rn %in% SweNeg, swe.full$rn, "")
	
	
	Uniq.Neg<-unique(c(ManNeg,MarNeg,SweNeg))
	Uniq.Neg<-c("BASV5371","BASV3385","BASV2392","FASV349","BASV497","FASV355","FASV410",
	"FASV2088","BASV76","BASV478","BASV267","BASV99","BASV31","FASV35")

		
	Man.count=read.table("Species_table_manresa.txt", sep="\t", row.names=1, header=1)
	Man.count=data.table(Man.count, keep.rownames=TRUE, key="rn")
	Man.count=Man.count[Man.count$rn %in% Uniq.Neg]
	Man.count=melt(Man.count)
	Man.count$Cultivar="Manresa"
		
	Mar.count=read.table("Species_table_marquis.txt", sep="\t", row.names=1, header=1)
	Mar.count=data.table(Mar.count, keep.rownames=TRUE, key="rn")
	Mar.count=Mar.count[Mar.count$rn %in% Uniq.Neg]
	Mar.count=melt(Mar.count)
	Mar.count$Cultivar="Marquis"
	
	Swe.count=read.table("Species_table_sweetann.txt", sep="\t", row.names=1, header=1)
	Swe.count=data.table(Swe.count, keep.rownames=TRUE, key="rn")
	Swe.count=Swe.count[Swe.count$rn %in% Uniq.Neg]
	Swe.count=melt(Swe.count)
	Swe.count$Cultivar="SweetAnn"
	
	Count=rbind(Man.count,Mar.count,Swe.count)
	
	DT.taxa1=read.table("taxa_Bac.txt", sep="\t", row.names=1, header=1)
	DT.taxa2=read.table("taxa_Fun.txt", sep="\t", row.names=1, header=1)
	DT.taxa=rbind(DT.taxa1,DT.taxa2)
	DT.taxa=data.table(DT.taxa,keep.rownames=TRUE, key="rn")
	DT.taxa=DT.taxa[DT.taxa$rn %in% Uniq.Neg]
	
	
	DT=merge(Count, DT.taxa, by.x="rn", by.y="rn")
	DT$mean <- with(DT, ave(value, rn))
	DT=data.table(DT, key="rn")
	DT=DT[order(mean,decreasing=F),]
	
	pp=ggplot(DT, aes(x=rn, y=log(value)))
	pp$data$rn <- ordered(pp$data$rn, levels= rev(unique(DT$rn)))
	p=pp+geom_boxplot(aes(color=Cultivar))+coord_flip()+
	scale_color_manual(values=list.color)+
	theme_bw()+theme_new
	
	
	gridExtra::grid.arrange(p, nrow=2, ncol=2)

#	for (i in Uniq.Neg) 
#	{
#		print ( i ) 
#		print(kruskal.test(data=p$data[p$data$rn == i ,], value ~ Cultivar))
#		print ("###########################################################")
#	}	



	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	



