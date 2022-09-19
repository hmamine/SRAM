#require packages for the analysis the analysis
pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "ape", "metagenomeSeq","PMCMR","UpSetR", "venn","seqtime")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

print("The scripts reporduce the figure indicated in Supp fig Supp fig 2 A")

Palette.phylum <-c("#EF5656","#47B3DA","#F7A415","#88888a","#000000","#2BB065")
ColPhylum<-c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")

Shape.1=c(21,19)

theme_new1 <- theme(
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	legend.position="none",
	)

theme_new2 <- theme(
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	legend.position="none",
	)

Manhat<-function(mapping){
		p<- ggplot( DT2 , mapping )
		p+geom_jitter( size=2 )+theme_bw( )+
		geom_hline(yintercept = -log2(0.01), linetype="longdash")+
		scale_color_manual(values=Palette.phylum)+
		scale_shape_manual(values=Shape.1)+
		coord_flip()+theme_new2
	}

Volcan<-function(mapping){
		p<- ggplot( DT2 , mapping )
		p+geom_jitter( size=2 )+ theme_bw( )+
		scale_color_manual(values=Palette.phylum)+
		geom_hline(yintercept = -log2(0.01), linetype="longdash")+
		geom_vline(xintercept = c(-1,0,1), linetype="longdash")+
		theme_new1
	}

DT<-vector( "list" )

#log2 fold change
	folch=1

#Adjusted p-value
	alpha=0.01

# upload and prepare phyloseq objects***
	tmp<-get(load("/home/mahassani/Biodata/NewHaven/SRAM/phyloseq_nochim_silva.RData"))
	OTU=tmp@otu_table
	TAXA=tmp@tax_table
	SD=tmp@sam_data
	TREE=tmp@phy_tree
	physeq=phyloseq(OTU, TAXA, SD,TREE) 
	
	physeq.BA=subset_taxa( physeq, Kingdom %in% c("Archaea","Bacteria"))
	physeq.BA=subset_taxa( physeq.BA, Genus != "Mitochondria")
	physeqB1=subset_samples( physeq.BA, Field == "Field1" )
	
	physeqB2=subset_samples( physeqB1, SampleType == "Rhizo" )
	
	DT.taxa=data.table( as(tax_table(physeqB2),"matrix"), keep.rownames=T, key="rn" )
	
	LIST<-list ("Manresa","Marquis","SweetAnne")
	for( i in LIST ) 
	{
		print( i )
		I=c(i)
		physeq.i=subset_samples( physeqB2, Cultivar == i )

###computing differentially abundant fungal ASVs
		m=t(as( otu_table( physeq.i ), "matrix" ) + 1L)
		t=data.frame( as( tax_table( physeq.i ), "matrix" ) )
		T=AnnotatedDataFrame( t )
		s=as( sample_data( physeq.i ), "data.frame" )
		S=AnnotatedDataFrame( s )
	
		obj=newMRexperiment( m, phenoData=S, featureData=T ) 
	
		p=cumNormStatFast( obj )
		objTrim=cumNorm( obj, p=p )
		Treatment = pData( obj )$Pathogen

		settings = zigControl( maxit=30, verbose=TRUE )
		dsg1=model.matrix( ~0+Pathogen, data=s )

		res1=fitZig( obj=objTrim, mod=dsg1, control=settings, useCSSoffset = TRUE)	
		zigFit1=res1@fit
		finalMod1=res1@fit$design
		c.mat1 = makeContrasts ( PathogenMac	-	Pathogennone, levels = finalMod1)
		fit1 = contrasts.fit( zigFit1, c.mat1 )
		fit1 = eBayes( fit1 )

		DT_1=fit1$coefficients
		DT_1p=fit1$p.value
		DT.zig<-topTable(fit1, number="inf",  adjust.method="BH" )
		DT.zig$rn<-rownames(DT.zig)
		DT.zig<-data.table(DT.zig, key="rn")
		DT3=DT.taxa[ DT.zig, ]
	
		ColPhylum<-c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria")
		DT3$ColPhylum <- ifelse(!DT3$Phylum %in% ColPhylum, "Other", ifelse(DT3$Phylum == "Actinobacteriota", "Actinobacteria", ifelse(DT3$Phylum == "Bacteroidota","Bacteroidetes", ifelse(DT3$Phylum == "Firmicutes","Firmicutes", "Proteobacteria"))))

		DT3$FC <- ifelse ( DT3$logFC > 0, "Enr", "Dep")
#		DT3$Col <- ifelse (DT3$adj.P.Val < alpha & abs( DT3$logFC ) > folch, DT3$ColPhylum,"NS" )
		DT3$Col <- ifelse (DT3$adj.P.Val < alpha, DT3$ColPhylum,"NS" )
		
		DT[[i]]<-DT3
	}

	PhylumFun <-c("#EF5656","#88888a","#000000","#2BB065")
	DT2<-DT$Manresa
	p3<-Manhat( aes( y=-log2(adj.P.Val), x=ColPhylum, color = Col, shape=FC )) +
	scale_color_manual(values=PhylumFun) +
	scale_y_continuous(limits=c(0,50))
	p4<-Volcan( aes( y=-log2(adj.P.Val), x=logFC, color = Col )) +
	scale_color_manual(values=PhylumFun)
	
	
	PhylumFun <-c("#EF5656","#88888a","#000000","#2BB065")
	DT2<-DT$Marquis
	p5<-Manhat( aes( y=-log2(adj.P.Val), x=ColPhylum, color = Col, shape=FC )) +
	scale_color_manual(values=PhylumFun) +
	scale_y_continuous(limits=c(0,50))
	p6<-Volcan( aes( y=-log2(adj.P.Val), x=logFC, color = Col )) +
	scale_color_manual(values=PhylumFun)
	
	
	PhylumFun <-c("#EF5656","#88888a","#000000","#2BB065")
	DT2<-DT$SweetAnne
	p7<-Manhat( aes( y=-log2(adj.P.Val), x=ColPhylum, color = Col, shape=FC )) +
	scale_color_manual(values=PhylumFun) +
	scale_y_continuous(limits=c(0,50))
	p8<-Volcan( aes( y=-log2(adj.P.Val), x=logFC, color = Col )) +
	scale_color_manual(values=PhylumFun)
	
	

#	pdf("Supp_fig_Manhattan.pdf",useDingbats=FALSE)
#	gridExtra::grid.arrange(p3,p5,p7, ncol=2, nrow=2)
#	dev.off()



