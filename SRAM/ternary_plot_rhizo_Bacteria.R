pkg=c("plyr", "ggplot2", "phyloseq", "data.table", "reshape2", "ape", "phytools","ggtern","PMCMR","venn", "metagenomeSeq","UpSetR", "venn","seqtime")

lapply(pkg, library, character.only = TRUE)


print (" These R scripts display the ternary plot indicated in Fig 3A and corresponding venn diagram in 2D")

rm(list = ls()[grep("*.*", ls())])
set.seed(131)
Palette.phylum <-c("#EF5656","#47B3DA","#F7A415","#88888a","#56137c","#2BB065")
ColPhylum<-c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria")
Shape.sign<-c(21, 16, 15)
lines <-data.frame(x=c(0.5,0,0.5),y=c(0.5,0.5,0),z=c(0,0.5,0.5), xend=c(1,1,1)/3, yend=c(1,1,1)/3, zend=c(1,1,1)/3)

minZ=1
maxZ=20
RAN=c(minZ,maxZ)

minX=1e-07
maxX=0.80
LIM=c(minX,maxX)

SEQ=seq(0.05,0.8,0.15)


Palette.phylum <-c("#EF5656","#47B3DA","#F7A415","#000000","#2BB065")
ColPhylum<-c("Actinobacteriota","Bacteroidota","Firmicutes","Proteobacteria")

add_seg<- geom_segment(data=lines, aes(x,y,z, xend=xend, yend=yend, zend=zend),color="grey", size=1,linetype="dotdash")
theme_change <- theme(
	legend.position="none",
	plot.background = element_blank(),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank()
	)
DT<-vector( "list" )
DT.out<-list()
RN<-list()

#log2 fold change
	folch=2
#adjusted p-value
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
	DT.taxa$ColPhylum <- ifelse(!DT.taxa$Phylum %in% ColPhylum, "Other", ifelse(DT.taxa$Phylum == "Actinobacteriota", "Actinobacteria", ifelse(DT.taxa$Phylum == "Bacteroidota","Bacteroidetes", ifelse(DT.taxa$Phylum == "Firmicutes","Firmicutes", "Proteobacteria"))))
	
	
	
	LIST<-list ("Manresa","Marquis","SweetAnne")

	MnMr=c("Manresa", "Marquis")
	MnSw=c("Manresa", "SweetAnne")
	MrSw=c("Marquis", "SweetAnne")

	LIST2<-list(MnMr,MnSw,MrSw)

###computing differentially ASVs === Manresa Vs Marquis	

	print("##### D.A. ASVs - Manresa Vs Marquis #####")	
	physeq.i=subset_samples( physeqB2, Cultivar %in% c("Manresa", "Marquis"))
	
#trimming OTUs with low occurance
	otu.I=t(otu_table(physeq.i))
	taxa.I=tax_table(physeq.i)
	sd.I=sample_data(physeq.i)
	flt=filterTaxonMatrix(otu.I,minocc=3, keepSum = TRUE, return.filtered.indices = TRUE)
	otu.I.filtered=flt$mat
	taxa.I.filtered=taxa.I[setdiff(1:nrow(taxa.I),flt$filtered.indices),]
	dummyTaxonomy=c("k__","p__","c__","o__","f__","g__","s__")
	TAXA.filt=rbind(taxa.I.filtered, dummyTaxonomy)
	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"
	rownames(otu.I.filtered)[nrow(otu.I.filtered)]="SUM"
	NewOTU=otu_table(otu.I.filtered, taxa_are_rows = TRUE)
	NewTaxa=tax_table(TAXA.filt)

	physeq.i=phyloseq(NewOTU, NewTaxa, sd.I)

#computing DE ASVs using fitZig	
	m=as( otu_table( physeq.i ), "matrix" ) + 1L
	t=data.frame( as( tax_table( physeq.i ), "matrix" ) )
	T=AnnotatedDataFrame( t )
	s=as( sample_data( physeq.i ), "data.frame" )
	S=AnnotatedDataFrame( s )
	obj=newMRexperiment( m, phenoData=S, featureData=T ) 
	p=cumNormStatFast( obj )
	objTrim=cumNorm( obj, p=p )
	Treatment = pData( obj )$Cultivar
	settings = zigControl( maxit=30, verbose=TRUE )
	dsg1=model.matrix( ~0+Cultivar, data=s )
	res1=fitZig( obj=objTrim, mod=dsg1, control=settings, useCSSoffset = TRUE)	
	zigFit1=res1@fit
	finalMod1=res1@fit$design
	c.mat1 = makeContrasts ( CultivarManresa	-	CultivarMarquis, levels = finalMod1)
	fit1 = contrasts.fit( zigFit1, c.mat1 )
	fit1 = eBayes( fit1 )
	DT_1=fit1$coefficients
	DT_1p=fit1$p.value
	DT.zig<-topTable(fit1, number="inf",  adjust.method="BH" )
	DT.zig$rn<-rownames(DT.zig)
	DT.zig<-data.table(DT.zig, key="rn")

	DT1=DT.taxa[ DT.zig, ]

	DT1$Col <- ifelse (DT1$adj.P.Val < alpha, DT1$ColPhylum,"NS" )
	
	Man1=DT1[DT1$adj.P.Val < alpha & DT1$logFC > folch]$rn
	Mar1=DT1[DT1$adj.P.Val < alpha & DT1$logFC < -folch]$rn
	
	
	
###computing differentially ASVs === Manresa Vs SweetAnne
	print("##### D.A. ASVs - Manresa Vs SweetAnne #####")	
	physeq.i=subset_samples( physeqB2, Cultivar %in% c("Manresa", "SweetAnne"))
	
#trimming OTUs with low occurance
	otu.I=t(otu_table(physeq.i))
	taxa.I=tax_table(physeq.i)
	sd.I=sample_data(physeq.i)
	flt=filterTaxonMatrix(otu.I,minocc=3, keepSum = TRUE, return.filtered.indices = TRUE)
	otu.I.filtered=flt$mat
	taxa.I.filtered=taxa.I[setdiff(1:nrow(taxa.I),flt$filtered.indices),]
	dummyTaxonomy=c("k__","p__","c__","o__","f__","g__","s__")
	TAXA.filt=rbind(taxa.I.filtered, dummyTaxonomy)
	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"
	rownames(otu.I.filtered)[nrow(otu.I.filtered)]="SUM"
	NewOTU=otu_table(otu.I.filtered, taxa_are_rows = TRUE)
	NewTaxa=tax_table(TAXA.filt)

	physeq.i=phyloseq(NewOTU, NewTaxa, sd.I)

#computing DE ASVs using fitZig	
	m=as( otu_table( physeq.i ), "matrix" ) + 1L
	t=data.frame( as( tax_table( physeq.i ), "matrix" ) )
	T=AnnotatedDataFrame( t )
	s=as( sample_data( physeq.i ), "data.frame" )
	S=AnnotatedDataFrame( s )
	obj=newMRexperiment( m, phenoData=S, featureData=T ) 
	p=cumNormStatFast( obj )
	objTrim=cumNorm( obj, p=p )
	Treatment = pData( obj )$Cultivar
	settings = zigControl( maxit=30, verbose=TRUE )
	dsg1=model.matrix( ~0+Cultivar, data=s )
	res1=fitZig( obj=objTrim, mod=dsg1, control=settings, useCSSoffset = TRUE)	
	zigFit1=res1@fit
	finalMod1=res1@fit$design
	c.mat1 = makeContrasts ( CultivarManresa	-	CultivarSweetAnne, levels = finalMod1)
	fit1 = contrasts.fit( zigFit1, c.mat1 )
	fit1 = eBayes( fit1 )
	DT_1=fit1$coefficients
	DT_1p=fit1$p.value
	DT.zig<-topTable(fit1, number="inf",  adjust.method="BH" )
	DT.zig$rn<-rownames(DT.zig)
	DT.zig<-data.table(DT.zig, key="rn")

	DT2=DT.taxa[ DT.zig, ]
	DT2$Col <- ifelse (DT2$adj.P.Val < alpha, DT2$ColPhylum,"NS" )
	
	Man2=DT2[DT2$adj.P.Val < alpha & DT2$logFC > folch]$rn
	Sw2=DT2[DT2$adj.P.Val < alpha & DT2$logFC < -folch]$rn
	
	
###computing differentially ASVs === Marquis Vs SweetAnne
	print("##### D.A. ASVs - Marquis Vs SweetAnne #####")	
	physeq.i=subset_samples( physeqB2, Cultivar %in% c("Marquis", "SweetAnne"))
	
#trimming OTUs with low occurance
	otu.I=t(otu_table(physeq.i))
	taxa.I=tax_table(physeq.i)
	sd.I=sample_data(physeq.i)
	flt=filterTaxonMatrix(otu.I,minocc=3, keepSum = TRUE, return.filtered.indices = TRUE)
	otu.I.filtered=flt$mat
	taxa.I.filtered=taxa.I[setdiff(1:nrow(taxa.I),flt$filtered.indices),]
	dummyTaxonomy=c("k__","p__","c__","o__","f__","g__","s__")
	TAXA.filt=rbind(taxa.I.filtered, dummyTaxonomy)
	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"
	rownames(otu.I.filtered)[nrow(otu.I.filtered)]="SUM"
	NewOTU=otu_table(otu.I.filtered, taxa_are_rows = TRUE)
	NewTaxa=tax_table(TAXA.filt)

	physeq.i=phyloseq(NewOTU, NewTaxa, sd.I)

#computing DE ASVs using fitZig	
	m=as( otu_table( physeq.i ), "matrix" ) + 1L
	t=data.frame( as( tax_table( physeq.i ), "matrix" ) )
	T=AnnotatedDataFrame( t )
	s=as( sample_data( physeq.i ), "data.frame" )
	S=AnnotatedDataFrame( s )
	obj=newMRexperiment( m, phenoData=S, featureData=T ) 
	p=cumNormStatFast( obj )
	objTrim=cumNorm( obj, p=p )
	Treatment = pData( obj )$Cultivar
	settings = zigControl( maxit=30, verbose=TRUE )
	dsg1=model.matrix( ~0+Cultivar, data=s )
	res1=fitZig( obj=objTrim, mod=dsg1, control=settings, useCSSoffset = TRUE)	
	zigFit1=res1@fit
	finalMod1=res1@fit$design
	c.mat1 = makeContrasts ( CultivarMarquis	-	CultivarSweetAnne, levels = finalMod1)
	fit1 = contrasts.fit( zigFit1, c.mat1 )
	fit1 = eBayes( fit1 )
	DT_1=fit1$coefficients
	DT_1p=fit1$p.value
	DT.zig<-topTable(fit1, number="inf",  adjust.method="BH" )
	DT.zig$rn<-rownames(DT.zig)
	DT.zig<-data.table(DT.zig, key="rn")

	DT3=DT.taxa[ DT.zig, ]
	DT3$Col <- ifelse (DT3$adj.P.Val < alpha, DT3$ColPhylum,"NS" )
	
	Mar3=DT3[DT3$adj.P.Val < alpha & DT3$logFC > folch]$rn
	Sw3=DT3[DT3$adj.P.Val < alpha & DT3$logFC < -folch]$rn
	
	
	Man<-unique(c(Man1,Man2))
	Mar<-unique(c(Mar1,Mar3))
	Sw<-unique(c(Sw2,Sw3))
	
#	pdf("venn_rhizo.pdf", useDingbats=FALSE)
	venn(list(Man,Mar,Sw), zcolor="style", snames=c("Manresa","Marquis","SweetAnne"))
#	dev.off()
# transformt counts to rel. abundance
	physeq.ra = transform_sample_counts(physeqB2, function(x) x/sum(x))

# computing arithmetic mean of rel. abund. by genotype
	DT.ra=data.table(t(otu_table(physeq.ra)), keep.rownames=T)
	DT.melt=data.table(melt(t(otu_table(physeq.ra))), keep.rownames=F)
	DT_sd=data.table(data.frame(sample_data(physeq.ra)), keep.rownames=T, key="rn")
	DT.merge=merge(DT.melt, DT_sd, by.x="Var2", by.y="rn")

	list.Man<-c(unique(DT.merge[DT.merge$Cultivar=="Manresa",]$Var2))
	DT.ra$Man.mean<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = list.Man]
	
	list.Mar<-c(unique(DT.merge[DT.merge$Cultivar=="Marquis",]$Var2))
	DT.ra$Mar.mean<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = list.Mar]
	
	list.Sw<-c(unique(DT.merge[DT.merge$Cultivar=="SweetAnne",]$Var2))
	DT.ra$Sw.mean<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = list.Sw]
	
	list.all<-c(unique(DT.merge[DT.merge$Cultivar %in% c("Manresa","Marquis","SweetAnne")]$Var2))
	DT.ra$mean<-DT.ra[,.(rowMeans(.SD,na.rm=TRUE)),.SDcols = list.all]
	
	DT.mean=DT.ra[,c("rn","Man.mean","Mar.mean","Sw.mean","mean")]
	DT.mean=merge(DT.mean, DT.taxa, by.x="rn", by.y="rn") 
	DT.mean=DT.mean[DT.mean$mean!=0]
	DT.mean=DT.mean[DT.mean$mean>0.00001]
	
	pp1=ggtern(DT.mean, mapping=aes(Man.mean, Mar.mean, Sw.mean ,colour=ColPhylum))
	p1=pp1 + geom_point(alpha=.8)+ theme_bw () + 
	scale_color_manual(values=Palette.phylum)+ 
	scale_size_continuous(range = RAN,limits=LIM, breaks=SEQ) + 
	theme_change + add_seg + labs(x="Manresa", y="Marquis", z="SweetAnne")

#	pdf("ternary_rhizo.pdf", useDingbats=FALSE)
	print(p1)
#	dev.off()

	VENN=venn(list(Man,Mar,Sw), zcolor="style", snames=c("Manresa","Marquis","SweetAnne"))
	ATTR<-attr(VENN,"intersection")
	
	
	uniq<-unique(c(Man,Mar,Sw))
	DT.uniq<-DT.taxa[DT.taxa$rn %in% uniq]
#	write.table(DT.uniq, "TAXA_root_Bac.txt", sep="\t")

	






























