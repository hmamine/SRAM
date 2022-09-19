#***require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "phytools", "metagenomeSeq","PMCMR")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])


print("These scripts display figures depicted in Fig 2C and Fig 2D"); 

color_palette1<-c("#4f7b95","#ffa500","#155832")
color_palette2<-c("#b02b76","#ffa500","#155832")
shape_palette=c(24,17,22,15)

theme_new <- theme (
	panel.grid.major = element_blank(),
	legend.position="none",
	panel.grid.minor = element_blank(),
	)
meth1=("PCoA")
meth2=("CAP")
dist1=("bray")

###upload bacterial data
	tmp<-get(load("/home/mahassani/Biodata/NewHaven/SRAM/phyloseq_nochim_silva.RData"))
	OTU.B=tmp@otu_table
	TAXA.B=tmp@tax_table
	SD.B=tmp@sam_data
	SD.B$Concat1=paste0(SD.B$SampleType,".",SD.B$Pathogen)
	TREE=tmp@phy_tree
	physeq.B=phyloseq(OTU.B, TAXA.B, SD.B,TREE) 

	physeq.BA=subset_taxa( physeq.B, Kingdom %in% c("Archaea","Bacteria") )
	physeq.BAM=subset_taxa( physeq.BA, Genus != "Mitochondria")
	physeq.BAM=subset_samples( physeq.BAM, SampleType != "Bulk" )
	TAXA.BAM=tax_table(physeq.BAM)
	
#normalization of count reads using CSS 
	otumat=t(as( otu_table( physeq.BAM ), "matrix" ))
	mp=newMRexperiment( ( otumat ) )
	physeq_norm=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), TAXA.B, SD.B,TREE )

	Fie1_Bac=subset_samples( physeq_norm, Field == "Field1" )

#computing Bray-Curtis distances
	dist_BC_Fie1=distance( Fie1_Bac, dist1 )

#computing unconstrained PCoA
	pcoa_BC_Fie1=ordinate( Fie1_Bac, meth1 , dist_BC_Fie1 )
	pFie1<-plot_ordination( Fie1_Bac, pcoa_BC_Fie1, shape="Concat1", color="Cultivar")
	pFie1$layers<-pFie1$layers[-1]
	p1=pFie1+geom_point(size=3)+theme_bw()+theme_new+
	scale_colour_manual(values=color_palette1)+
	scale_shape_manual(values=shape_palette)

#computing constrained PCoA	
	Bac.pc.BC=ordinate( Fie1_Bac, meth2 , dist_BC_Fie1, ~Cultivar*Pathogen*SampleType )
	cpBC<-plot_ordination( Fie1_Bac, Bac.pc.BC, shape="Concat1", color="Cultivar")
	cpBC$layers<-cpBC$layers[-1]
	p2=cpBC+geom_point(size=3)+theme_bw()+theme_new+
	scale_colour_manual(values=color_palette1)+
	scale_shape_manual(values=shape_palette)
	
	print(anova.cca(Bac.pc.BC, by="term", perm.max=1000))

#	DT.L=data.table( as( sample_data( physeq_norm.L ), "data.frame" ),keep.rownames=T, key="rn" )
#	with( DT.L, adonis ( dist_BC.L ~ Geno ) )
	
#	print(p3)

###upload fungal data 
	tmp<-get(load("/home/mahassani/Biodata/NewHaven/SRAM/phyloseq_nochim_unite.RData"))
	OTU.f=tmp@otu_table
	TAXA.f=tmp@tax_table
	SD.f=tmp@sam_data
	SD.f$Concat1=paste0(SD.f$SampleType,".",SD.f$Pathogen)
	physeq.f=phyloseq(OTU.f, TAXA.f, SD.f) 
	physeq.Fun=subset_taxa( physeq.f, Kingdom == "k__Fungi")
	physeq.Fun=subset_samples( physeq.Fun, SampleType != "Bulk" )
	TAXA.Fun=tax_table(physeq.Fun)
#normalization of count reads using CSS 
	otumat=t(as( otu_table( physeq.Fun ), "matrix" ))
	mp=newMRexperiment( ( otumat ) )
	Norm_Fun=phyloseq( otu_table( MRcounts( cumNorm( mp, p=cumNormStat( mp ) ),
	norm=TRUE, log=TRUE ),taxa_are_rows=TRUE ), TAXA.Fun, SD.f)

	Fie1_Fun=subset_samples( Norm_Fun, Field == "Field1" )

#computing Bray-Curtis distances
	BC_Fie1_Fun=distance( Fie1_Fun, dist1 )
	
#computing unconstrained PCoA
	pcoa_Fie1_Fun=ordinate( Fie1_Fun, meth1 , BC_Fie1_Fun )
	pF1F<-plot_ordination( Fie1_Fun, pcoa_Fie1_Fun, shape="Concat1", color="Cultivar")
	pF1F$layers<-pF1F$layers[-1]
	p3=pF1F+geom_point(size=3)+theme_bw()+theme_new+
	scale_colour_manual(values=color_palette1)+
	scale_shape_manual(values=shape_palette)
	
#computing constrained PCoA	
	Fun.pc.BC=ordinate( Fie1_Fun, meth2 , BC_Fie1_Fun, ~Cultivar*Pathogen*SampleType )
	cpBC.Fun<-plot_ordination( Fie1_Fun, Fun.pc.BC, shape="Concat1", color="Cultivar")
	cpBC.Fun$layers<-cpBC.Fun$layers[-1]
	p4=cpBC.Fun+geom_point(size=3)+theme_bw()+theme_new+
	scale_colour_manual(values=color_palette1)+
	scale_shape_manual(values=shape_palette)
	
	print(anova.cca(Fun.pc.BC, by="term", perm.max=1000))
	
#	print(p3)	
	
#	pdf("figure_cpcoa.pdf",useDingbats=FALSE)
	gridExtra::grid.arrange(p2,p4, nrow=2, ncol=2)
#	dev.off()

	


