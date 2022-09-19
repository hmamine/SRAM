#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table","agricolae", "metagenomeSeq", "ape", "vegan", "dplyr", "seqtime")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])
print(" R script to compute BC distance to centroids indicated in Fig 1C and 1D");
color_palette1<-c("#4f7b95","#ffa500","#000000","#155832")
color_palette2<-c("#b02b76","#ffa500","#000000","#155832")
theme_new <- theme (
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	axis.ticks.x=element_blank(),
#	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
#	legend.position="none",
	axis.text.x = element_text(angle=90, vjust=1, size=9, color="black"),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	)
dist=( "bray" )
Ord1=c("Bulk.none.none","Bulk.Mac.none",
	"Rhizo.none.Manresa","Rhizo.Mac.Manresa",
	"Rhizo.none.Marquis","Rhizo.Mac.Marquis",
	"Rhizo.none.SweetAnne","Rhizo.Mac.SweetAnne",
	"Root.none.Manresa","Root.Mac.Manresa",
	"Root.none.Marquis","Root.Mac.Marquis",
	"Root.none.SweetAnne","Root.Mac.SweetAnne")
	
Ord2=c("Bulk.none.none","Bulk.Vert.none",
	"Rhizo.none.Festival","Rhizo.Vert.Festival",
	"Rhizo.none.Marquis","Rhizo.Vert.Marquis",
	"Rhizo.none.SweetAnne","Rhizo.Vert.SweetAnne",
	"Root.none.Festival","Root.Vert.Festival",
	"Root.none.Marquis","Root.Vert.Marquis",
	"Root.none.SweetAnne","Root.Vert.SweetAnne")

###upload bacterial data
	tmp<-get(load( "/home/mahassani/Biodata/NewHaven/SRAM/phyloseq_nochim_silva.RData") )
	OTU.B=tmp@otu_table
	TAXA.B=tmp@tax_table
	SD.B=tmp@sam_data
	SD.B$Concat1=paste0(SD.B$SampleType,".",SD.B$Pathogen,".",SD.B$Cultivar)
	TREE=tmp@phy_tree
	physeq.B=phyloseq(OTU.B, TAXA.B, SD.B,TREE) 

#	physeq.B=subset_samples( physeq.B, Field != "Nursery" )
	physeq.BA=subset_taxa( physeq.B, Kingdom %in% c("Archaea","Bacteria") )
	physeq.BAM=subset_taxa( physeq.BA, Genus != "Mitochondria")
	TAXA.BAM=tax_table(physeq.BAM)
	
#normalization of count reads using CSS 
	otumat=t(as( otu_table( physeq.BAM ), "matrix" ) )
	mp=newMRexperiment( ( otumat ) )
	Norm.BA=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), TAXA.BAM, SD.B,TREE )

	Fie1.BA=subset_samples( Norm.BA, Field == "Field1" )	
	
#computing distance to centroid
	dist.B.Fie1=distance( Fie1.BA, dist )	
	SD.B.Fie1=data.table(as(sample_data(Fie1.BA), "data.frame"),keep.rownames=T, key="rn")	
	disp.B.Fie1 <- with( SD.B.Fie1, betadisper( dist.B.Fie1, Concat1, type="centroid" ) )	
#	print(disp)
#	print(anova(disp))
#	print(TukeyHSD(disp))
	DT.B.Fie1=data.table(as.data.frame(disp.B.Fie1$distances), keep.rownames=T, key="rn")
	setnames(DT.B.Fie1, "disp.B.Fie1$distances","value")
#	DT.sd=data.table( as( Sd, "data.frame" ),keep.rownames=T, key="rn" )
	DT.B.Fie1=DT.B.Fie1[SD.B.Fie1]
	pp1=ggplot(DT.B.Fie1, aes(x=Concat1, y=value))
	pp1$data$Concat1 <- ordered(pp1$data$Concat1, levels=Ord1 )
	p1=pp1+geom_boxplot(aes(color=Cultivar),width=0.75)+
	geom_point(size=.5, color="black")+theme_bw()+ theme_new +
	scale_color_manual(values=color_palette1)+
	scale_y_continuous( limits=c( .2,.65  ), breaks=c(.2,.3,.4,.5,.6,.7) )+
	ylab("BC distances to centroid")
	
kruskal.test(data=p1$data[p1$data$SampleType == "Bulk",], value ~ Pathogen)
kruskal.test(data=p1$data[p1$data$SampleType == "Rhizo" & p1$data$Cultivar == "Manresa",], value ~ Pathogen)
kruskal.test(data=p1$data[p1$data$SampleType == "Rhizo" & p1$data$Cultivar == "Marquis",], value ~ Pathogen)
kruskal.test(data=p1$data[p1$data$SampleType == "Rhizo" & p1$data$Cultivar == "SweetAnne",], value ~ Pathogen)
kruskal.test(data=p1$data[p1$data$SampleType == "Root" & p1$data$Cultivar == "Manresa",], value ~ Pathogen)
kruskal.test(data=p1$data[p1$data$SampleType == "Root" & p1$data$Cultivar == "Marquis",], value ~ Pathogen)
kruskal.test(data=p1$data[p1$data$SampleType == "Root" & p1$data$Cultivar == "SweetAnne",], value ~ Pathogen)
###upload fungal data
	tmp<-get(load("/home/mahassani/Biodata/NewHaven/SRAM/phyloseq_nochim_unite.RData"))
	OTU.F=tmp@otu_table
	TAXA.F=tmp@tax_table
	SD.F=tmp@sam_data
	SD.F$Concat1=paste0(SD.F$SampleType,".",SD.F$Pathogen,".",SD.F$Cultivar)
	physeq.F=phyloseq(OTU.F, TAXA.F, SD.F) 

	physeq.Fun=subset_taxa( physeq.F, Kingdom == "k__Fungi")
	TAXA.Fun=tax_table(physeq.Fun)
	
#normalization of count reads using CSS 
	otumat=t(as( otu_table( physeq.Fun), "matrix" ) )
	mp=newMRexperiment( ( otumat ) )
	Norm.Fun=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), TAXA.Fun, SD.F)

	Fie1.Fun=subset_samples( Norm.Fun, Field == "Field1" )
	
#computing distance to centroid
	dist.F.Fie1=distance( Fie1.Fun, dist )	
	SD.F.Fie1=data.table(as(sample_data(Fie1.Fun), "data.frame"),keep.rownames=T, key="rn")	
	disp.F.Fie1 <- with( SD.F.Fie1, betadisper( dist.F.Fie1, Concat1, type="centroid" ) )	
#	print(disp)
#	print(anova(disp))
#	print(TukeyHSD(disp))
	DT.F.Fie1=data.table(as.data.frame(disp.F.Fie1$distances), keep.rownames=T, key="rn")
	setnames(DT.F.Fie1, "disp.F.Fie1$distances","value")
#	DT.sd=data.table( as( Sd, "data.frame" ),keep.rownames=T, key="rn" )
	DT.F.Fie1=DT.F.Fie1[SD.F.Fie1]
	pp3=ggplot(DT.F.Fie1, aes(x=Concat1, y=value))
	pp3$data$Concat1 <- ordered(pp3$data$Concat1, levels=Ord1 )
	p3=pp3+geom_boxplot(aes(color=Cultivar),width=0.75)+
	geom_point(size=.5, color="black")+theme_bw()+ theme_new +
	scale_color_manual(values=color_palette1)+
	scale_y_continuous( limits=c( .2,.65  ), breaks=c(.2,.3,.4,.5,.6,.7) )+
	ylab("BC distances to centroid")

kruskal.test(data=p3$data[p3$data$SampleType == "Bulk",], value ~ Pathogen)
kruskal.test(data=p3$data[p3$data$SampleType == "Rhizo" & p3$data$Cultivar == "Manresa",], value ~ Pathogen)
kruskal.test(data=p3$data[p3$data$SampleType == "Rhizo" & p3$data$Cultivar == "Marquis",], value ~ Pathogen)
kruskal.test(data=p3$data[p3$data$SampleType == "Rhizo" & p3$data$Cultivar == "SweetAnne",], value ~ Pathogen)
kruskal.test(data=p3$data[p3$data$SampleType == "Root" & p3$data$Cultivar == "Manresa",], value ~ Pathogen)
kruskal.test(data=p3$data[p3$data$SampleType == "Root" & p3$data$Cultivar == "Marquis",], value ~ Pathogen)
kruskal.test(data=p3$data[p3$data$SampleType == "Root" & p3$data$Cultivar == "SweetAnne",], value ~ Pathogen)

#	pdf("figure_3.pdf",useDingbats=FALSE)
	gridExtra::grid.arrange(p1,p3, nrow=2, ncol=2)	
#	dev.off()

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	



