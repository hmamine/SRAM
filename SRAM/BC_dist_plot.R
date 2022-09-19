#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table","agricolae", "metagenomeSeq", "ape", "vegan", "dplyr", "seqtime","picante","PMCMRplus")
lapply(pkg, library, character.only = TRUE)


print("These scripts compute BC distance displayed in Fig 2E and 2F") ;

rm(list = ls()[grep("*.*", ls())])
set.seed(131)
color_palette1<-c("#4f7b95","#ffa500","#155832")
color_palette2<-c("#b02b76","#ffa500","#155832")
theme_new <- theme (
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	axis.ticks.x=element_blank(),
	axis.text.x = element_blank(),
#	axis.text.y = element_blank(),
	legend.position="none",
#	axis.text.x = element_text(angle=90, vjust=1, size=9, color="black"),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	)
dist=( "bray" )

###upload bacterial data
	tmp<-get(load( "/home/mahassani/Biodata/NewHaven/SRAM/phyloseq_nochim_silva.RData") )
	OTU.B=tmp@otu_table
	TAXA.B=tmp@tax_table
	SD.B=tmp@sam_data
	SD.B$Concat1=paste0(SD.B$SampleType,".",SD.B$Pathogen,".",SD.B$Cultivar)
	TREE=tmp@phy_tree
	physeq.B=phyloseq(OTU.B, TAXA.B, SD.B,TREE) 

	physeq.BA=subset_taxa( physeq.B, Kingdom %in% c("Archaea","Bacteria") )
	physeq.BAM=subset_taxa( physeq.BA, Genus != "Mitochondria")
	physeq.Mp=subset_samples( physeq.BAM, Field == "Field1" )

#trimming OTUs with low occurance
	otu.I=t(otu_table(physeq.Mp))
	taxa.I=tax_table(physeq.Mp)
	sd.I=sample_data(physeq.Mp)

	flt=filterTaxonMatrix(otu.I,minocc=3, keepSum = TRUE, return.filtered.indices = TRUE)
	otu.I.filtered=flt$mat
	taxa.I.filtered=taxa.I[setdiff(1:nrow(taxa.I),flt$filtered.indices),]
	dummyTaxonomy=c("k__dummy","p__","c__","o__","f__","g__","s__","sc__","ot_")
	TAXA.filt=rbind(taxa.I.filtered, dummyTaxonomy)
	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"
	rownames(otu.I.filtered)[nrow(otu.I.filtered)]="SUM"
	NewOTU=otu_table(otu.I.filtered, taxa_are_rows = TRUE)
	NewTaxa=tax_table(TAXA.filt)
	physeq.filtB=phyloseq(NewOTU, NewTaxa, sd.I)

	
#normalization of count reads using CSS 
	otumat=as( otu_table( physeq.filtB ), "matrix" ) 
	mp=newMRexperiment( ( otumat ) )
	physeq.norm=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), NewTaxa, sd.I)
	
	DT.sd=data.table(as(sample_data(physeq.norm), "data.frame"),keep.rownames=T, key="rn")
	dist.BC=distance(physeq.norm, dist )
	melt.dist.BC=reshape2::melt(as.matrix( dist.BC ))
	
	DT.dist = melt.dist.BC %>%
   	filter(as.character(Var1) != as.character(Var2)) %>%
   	mutate_if(is.factor,as.character)
 	
   	DT.dist=data.table(DT.dist ,  key=c("Var1","Var2"))
   	DT.merge=merge( DT.dist, DT.sd, by.x="Var2", by.y="rn" )
	DT.merge=merge( DT.merge, DT.sd, by.x="Var1", by.y="rn" )
	
	DT.rhizo=DT.merge[DT.merge$SampleType.x=="Rhizo" & DT.merge$SampleType.y=="Rhizo"  ]
	DT.rhizo=DT.rhizo[DT.rhizo$Cultivar.x=="Manresa"]
	p1=ggplot(DT.rhizo, aes(x=Cultivar.y, y=value))+
	geom_violin(aes(fill=Cultivar.y))+geom_boxplot(width=0.1,outlier.size =0.5)+
	scale_fill_manual(values=color_palette1)+
	theme_bw()+ theme_new + ylab("BC distances")
#	kruskal.test(data=p1$data, value ~ Cultivar.y)
#	kwAllPairsConoverTest(value ~ as.factor(Cultivar.y) , data=p1$data,  p.adjust.method="BH" )
#	aov.out<-aov(data=p1$data, value ~ Cultivar.y)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Cultivar.y"))
	
	DT.root=DT.merge[DT.merge$SampleType.x=="Root" & DT.merge$SampleType.y=="Root"  ]
	DT.root=DT.root[DT.root$Cultivar.x=="Manresa"]
	p2=ggplot(DT.root, aes(x=Cultivar.y, y=value))+
	geom_violin(aes(fill=Cultivar.y))+geom_boxplot(width=0.1,outlier.size =0.5)+
	scale_fill_manual(values=color_palette1)+
	theme_bw()+ theme_new + ylab("BC distances")
#	aov.out<-aov(data=p2$data, value ~ Cultivar.y)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Cultivar.y"))

	
	
###upload fungal data
	tmp<-get(load("/home/mahassani/Biodata/NewHaven/SRAM/phyloseq_nochim_unite.RData"))
	OTU.F=tmp@otu_table
	TAXA.F=tmp@tax_table
	SD.F=tmp@sam_data
	SD.F$Concat1=paste0(SD.F$SampleType,".",SD.F$Pathogen,".",SD.F$Cultivar)
	physeq.F=phyloseq(OTU.F, TAXA.F, SD.F) 
	physeq.Fun=subset_taxa( physeq.F, Kingdom == "k__Fungi")
	physeq.Mp=subset_samples( physeq.Fun, Field == "Field1" )
	
#trimming OTUs with low occurance
	otu.I=t(otu_table(physeq.Mp))
	taxa.I=tax_table(physeq.Mp)
	sd.I=sample_data(physeq.Mp)

	flt=filterTaxonMatrix(otu.I,minocc=3, keepSum = TRUE, return.filtered.indices = TRUE)
	otu.I.filtered=flt$mat
	taxa.I.filtered=taxa.I[setdiff(1:nrow(taxa.I),flt$filtered.indices),]
	dummyTaxonomy=c("k__dummy","p__","c__","o__","f__","g__","s__","sc__","ot_")
	TAXA.filt=rbind(taxa.I.filtered, dummyTaxonomy)
	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"
	rownames(otu.I.filtered)[nrow(otu.I.filtered)]="SUM"
	NewOTU=otu_table(otu.I.filtered, taxa_are_rows = TRUE)
	NewTaxa=tax_table(TAXA.filt)
	physeq.filtF=phyloseq(NewOTU, NewTaxa, sd.I)

#normalization of count reads using CSS 
	otumat=as( otu_table( physeq.filtF ), "matrix" ) 
	mp=newMRexperiment( ( otumat ) )
	physeq.norm=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), NewTaxa, sd.I)


	DT.sd=data.table(as(sample_data(physeq.norm), "data.frame"),keep.rownames=T, key="rn")
	dist.BC=distance(physeq.norm, dist )
	melt.dist.BC=reshape2::melt(as.matrix( dist.BC ))
	
	DT.dist = melt.dist.BC %>%
   	filter(as.character(Var1) != as.character(Var2)) %>%
   	mutate_if(is.factor,as.character)
 	
   	DT.dist=data.table(DT.dist ,  key=c("Var1","Var2"))
   	DT.merge=merge( DT.dist, DT.sd, by.x="Var2", by.y="rn" )
	DT.merge=merge( DT.merge, DT.sd, by.x="Var1", by.y="rn" )
	
	DT.rhizo=DT.merge[DT.merge$SampleType.x=="Rhizo" & DT.merge$SampleType.y=="Rhizo"  ]
	DT.rhizo=DT.rhizo[DT.rhizo$Cultivar.x=="Manresa"]
	p3=ggplot(DT.rhizo, aes(x=Cultivar.y, y=value))+
	geom_violin(aes(fill=Cultivar.y))+geom_boxplot(width=0.1,outlier.size =0.5)+
	scale_fill_manual(values=color_palette1)+
	theme_bw()+ theme_new + ylab("BC distances")
#	aov.out<-aov(data=p3$data, value ~ Cultivar.y)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Cultivar.y"))	

	DT.root=DT.merge[DT.merge$SampleType.x=="Root" & DT.merge$SampleType.y=="Root"  ]
	DT.root=DT.root[DT.root$Cultivar.x=="Manresa" ]
	p4=ggplot(DT.root, aes(x=Cultivar.y, y=value))+
	geom_violin(aes(fill=Cultivar.y))+geom_boxplot(width=0.1,outlier.size =0.5)+
	scale_fill_manual(values=color_palette1)+
	theme_bw()+ theme_new + ylab("BC distances")
#	aov.out<-aov(data=p4$data, value ~ Cultivar.y)
#	print(summary(aov.out))
#	print(TukeyHSD(aov.out))
#	print(HSD.test(aov.out,trt="Cultivar.y"))

#	pdf("BC_distances.pdf",useDingbats=FALSE)
	gridExtra::grid.arrange(p1,p2,p3,p4, nrow=2, ncol=2)
#	dev.off()
	
	
	
	
	
	
	

