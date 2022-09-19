#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table","agricolae", "metagenomeSeq", "ape", "vegan", "dplyr", "seqtime","picante","PMCMRplus")
lapply(pkg, library, character.only = TRUE)


print(" These scripts computed permanova analysis indicated in Table-1 ");

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

#computing permanova on BC distances
	physeq.Rhizo=subset_samples( physeq.norm, SampleType == "Rhizo" )
	DT.MM=data.table( as( sample_data( physeq.Rhizo ), "data.frame" ),keep.rownames=T, key="rn" )
	DT.MM$test<-ifelse(DT.MM$Cultivar %in% c("Manresa","Marquis"), "CC","SA")
	BC.MM=distance(physeq.Rhizo, dist )
	print("PERMANOVA for *** Bacterial Comm. in Rhizo. ***")
	print(with( DT.MM, adonis2 ( BC.MM ~ test ) ))
	print("###########################################################################")
	physeq.Root=subset_samples( physeq.norm, SampleType == "Root" )
	DT.MM=data.table( as( sample_data( physeq.Root ), "data.frame" ),keep.rownames=T, key="rn" )
	DT.MM$test<-ifelse(DT.MM$Cultivar %in% c("Manresa","Marquis"), "CC","SA")
	BC.MM=distance(physeq.Root, dist )
	print("##### PERMANOVA for *** Bacterial Comm. in Roots *** #####")
	print(with( DT.MM, adonis2 ( BC.MM ~ test ) ))
	print("###########################################################################")

	
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


#computing permanova on BC distances
	physeq.Rhizo=subset_samples( physeq.norm, SampleType == "Rhizo" )
	DT.MM=data.table( as( sample_data( physeq.Rhizo ), "data.frame" ),keep.rownames=T, key="rn" )
	DT.MM$test<-ifelse(DT.MM$Cultivar %in% c("Manresa","Marquis"), "CC","SA")
	BC.MM=distance(physeq.Rhizo, dist )
	print("##### PERMANOVA for ***** Fungal Comm. in Rhizosphere ***** #####")
	print(with( DT.MM, adonis2 ( BC.MM ~ test ) ))
	print("###########################################################################")
	physeq.Root=subset_samples( physeq.norm, SampleType == "Root" )
	DT.MM=data.table( as( sample_data( physeq.Root ), "data.frame" ),keep.rownames=T, key="rn" )
	DT.MM$test<-ifelse(DT.MM$Cultivar %in% c("Manresa","Marquis"), "CC","SA")
	BC.MM=distance(physeq.Root, dist )
	print("##### PERMANOVA for ***** Fungal Comm. in Roots ***** #####")
	print(with( DT.MM, adonis2 ( BC.MM ~ test ) ))
	print("###########################################################################")
	
	
	
	
	
	

