#require packages for the analysis the analysis
pkg=c("ggplot2", "data.table", "reshape2","PMCMRplus")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

print(" These R scripts reproduce the plottings indicated in Fig 2A and 2B")

list.color<-c("#4f7b95","#ffa500","#155832")

theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	legend.position="none",
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x = element_text(angle=90, vjust=1, size=9, color="black")
	)

tab=read.table("qPCR_Mp.txt", sep="\t", row.names=1, header=1)
tab=data.table(tab, keep.rownames=T, key="rn")
DT<-as.data.table( melt( tab ))

pp5<-ggplot(data=DT, aes(x=Concat, y=value, colour=Cv))
p5=pp5+geom_boxplot()+
	geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=Cv),binwidth = 0.02,show.legend=F)+
	scale_color_manual(values=list.color)+ theme_bw()+	
	scale_fill_manual(values=list.color) + theme_new 



tab1=read.table("data_survival.txt", sep="\t", row.names=1, header=1)
p4<-ggplot(tab1, aes(x=Plot.ID, y=survival, colour=Plot.ID))+ geom_boxplot()+
	geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=Plot.ID),binwidth = 1,show.legend=F)+
	scale_color_manual(values=list.color)+ theme_bw()+	
	scale_fill_manual(values=list.color) + theme_new 
	
	

		
#	pdf("Supp_fig_qPCR_Mp_Roots.pdf",useDingbats=FALSE)
	gridExtra::grid.arrange(p4,p5, nrow=2, ncol=2)
#	dev.off()

kruskal.test(data=p4$data, survival ~ Plot.ID)
kwAllPairsConoverTest(survival ~ as.factor(Plot.ID),data=p4$data,  p.adjust.method="BH" )

kruskal.test(data=p5$data, value ~ Concat)
kwAllPairsConoverTest(value ~ as.factor(Concat) , data=p5$data,  p.adjust.method="BH" )









