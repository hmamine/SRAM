#require packages for the analysis the analysis
pkg=c( "ggplot2", "phyloseq", "data.table", "reshape2", "seqtime" ,"SpiecEasi", "huge")
lapply(pkg, library, character.only = TRUE)


print (" These R scripts to plot Fig 4D showing within-module and between-module connectivity scores");

rm(list = ls()[grep("*.*", ls())])
set.seed(131)
Palette.phylum <-c("#EF5656","#47B3DA","#F7A415","#88888a","#2BB065")

theme_change <- theme(
	legend.position="none",
	plot.background = element_blank(),
	panel.grid.minor = element_blank(),
	panel.grid.major = element_blank()
	)

# upload and prepare phyloseq objects***
	##upload and prepare phyloseq objects bacteria data
#	tmp<-get(load("/home/mahassani/Biodata/SRAM/phyloseq_nochim_silva.RData"))
#	OTU=tmp@otu_table
#	TAXA=tmp@tax_table
#	SD=tmp@sam_data
#	TREE=tmp@phy_tree
#	physeq1=phyloseq(OTU, TAXA, SD,TREE) 
	
#	physeq.BA=subset_taxa( physeq1, Kingdom %in% c("Archaea","Bacteria"))
#	physeq.BA=subset_taxa( physeq.BA, Genus != "Mitochondria")
#	physeq.Bac=subset_samples( physeq.BA, Field == "Field1" )
#	physeq.BacRA=transform_sample_counts(physeq.Bac, function(x) x / sum(x) )
	
#	write.table(tax_table(physeq.BacRA),"taxa_Bac.txt", sep="\t")
#	write.table(t(otu_table(physeq.BacRA)),"otu_Bac.txt", sep="\t")
#	write.table(sample_data(physeq.BacRA),"sample_dataB.txt", sep="\t")

	mat=read.table( "otu_Bac.txt", sep="\t", row.names=1, header=T)
	mat=as.matrix(mat)
	OTU=otu_table(mat, taxa_are_rows=T) 
	tax=read.table("taxa_Bac.txt", sep="\t", row.names=1, header=1)
	tax=as.matrix(tax)
	TAXA=tax_table(tax)
	sd=read.table("sample_data_Bac.txt", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)

	physeq.Bac=phyloseq(OTU, TAXA, SD) 
	
##upload and prepare phyloseq objects fungi data
#	tmp<-get(load("/home/mahassani/Biodata/SRAM/phyloseq_nochim_unite.RData"))
#	OTU=tmp@otu_table
#	TAXA=tmp@tax_table
#	SD=tmp@sam_data
#	physeq2=phyloseq(OTU, TAXA, SD) 
	
#	physeq.Fun=subset_taxa( physeq2, Kingdom == "k__Fungi")
#	physeq.Fun=subset_samples( physeq.Fun, Field == "Field1" )
#	physeq.FunRA=transform_sample_counts(physeq.Fun, function(x) x / sum(x) )

#	write.table(tax_table(physeq.FunRA),"taxa_Fun.txt", sep="\t")
#	write.table(t(otu_table(physeq.FunRA)),"otu_Fun.txt", sep="\t")
#	write.table(sample_data(physeq.FunRA),"sample_dataF.txt", sep="\t")

	mat=read.table( "otu_Fun.txt", sep="\t", row.names=1, header=T)
	mat=as.matrix(mat)
	OTU=otu_table(mat, taxa_are_rows=T) 
	tax=read.table("taxa_Fun.txt", sep="\t", row.names=1, header=1)
	tax=as.matrix(tax)
	TAXA=tax_table(tax)
	sd=read.table("sample_data_Fun.txt", sep="\t", row.names=1, header=1)
	SD=sample_data(sd)
	physeq.Fun=phyloseq(OTU, TAXA, SD) 
	physeq.merge<-merge_phyloseq(physeq.Bac, physeq.Fun)
	physeq.merge<-merge_phyloseq(physeq.Bac, physeq.Fun)
	physeq.merge=tax_glom(physeq.merge, "Species")

#subset data to Manresa
	physeq.subset=subset_samples(physeq.merge, Cultivar == "Manresa")

##filter OTUs with low occurence across samples
	otu.cv.ts=otu_table(physeq.subset)
	taxa.cv.ts=tax_table(physeq.subset)
	sd.cv.ts=sample_data(physeq.subset)
	flt=filterTaxonMatrix(otu.cv.ts,minocc=6, keepSum = TRUE, return.filtered.indices = TRUE)
	otu.cv.ts.filtered=flt$mat
	taxa.cv.ts.filtered=taxa.cv.ts[setdiff(1:nrow(taxa.cv.ts),flt$filtered.indices),]
	dummyTaxonomy=c("k__","p__","c__","o__","f__","g__","s__")
	TAXA.filt=rbind(taxa.cv.ts.filtered, dummyTaxonomy)
	ColPhylum<-c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Basidiomycota", "Ascomycota")
	TAXA.filt[,"Phylum"] <- ifelse(! TAXA.filt[,"Phylum"] %in% ColPhylum, "Other", ifelse(TAXA.filt[,"Phylum"] == "Bacteroidota", "Bacteroidetes", ifelse(TAXA.filt[,"Phylum"] == "Firmicutes", "Firmicutes", ifelse (TAXA.filt[,"Phylum"] == "Actinobacteriota", "Actinobacteria", ifelse (TAXA.filt[,"Phylum"] == "Proteobacteria","Proteobacteria", ifelse (TAXA.filt[,"Phylum"] == "Ascomycota", "Ascomycota", "Basidiomycota" ))))))
	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"	
	rownames(otu.cv.ts.filtered)[nrow(otu.cv.ts.filtered)]="SUM"
	NewOTU=otu_table(otu.cv.ts.filtered, taxa_are_rows = TRUE)	
	NewTaxa=tax_table(TAXA.filt)
	physeq.filter=phyloseq(NewOTU, NewTaxa, sd.cv.ts)
#computing bacterial and fungal co-occurrences network 
print("##### computing newtwork using mb method #####")
	spiec.man <- spiec.easi(physeq.filter, method="mb", verbose=TRUE, nlambda=99)
	ig2.man <- adj2igraph(getRefit(spiec.man), vertex.attr=list(name=taxa_names(physeq.filter)))
#computing within- (Zi) and among-modules (Pi) connectivity scores - Adapted from Guo et al 2022
	seqdeg=degree(ig2.man)
	Nnodes=length(seqdeg)
	Z=seqdeg
	Z[]=0
	P=Z
	fc.man = cluster_edge_betweenness(ig2.man,weights =NULL)
	Membership=membership(fc.man)
	Seq=seq(1:Nnodes)
	for(i in 1:Nnodes){
		L=Membership==Membership[i]         
		neighbs=neighbors(ig2.man, i)               
		Kis=sum(L[neighbs])
		SUM=0
		SUMsq=0	
		SUMP=0
		Miv=Seq[L]
		for(j in 1:sum(L)){
			neighbsj=neighbors(ig2.man, Miv[j])
			Kjs=sum(L[neighbsj])
			SUM=SUM+Kjs
			SUMsq=SUMsq+Kjs^2
		}
		Z[i]=(Kis-SUM/sum(L))/sqrt(SUMsq/sum(L)-(SUM/sum(L))^2)
		if(Kis-SUM/sum(L)==0){Z[i]=0}
		for(k in 1:max(Membership)){
		Lp=Membership==k
		Kisp=sum(Lp[neighbs])
		SUMP=SUMP+(Kisp/seqdeg[i])^2}
		P[i]=1-SUMP
	}
	attr.man=cbind(degree=seqdeg,module=Membership,Pi=P,Zi=Z)

#subset data to Manrquis
	physeq.subset=subset_samples(physeq.merge, Cultivar == "Marquis")
##filter OTUs with low occurence across samples
	otu.cv.ts=otu_table(physeq.subset)
	taxa.cv.ts=tax_table(physeq.subset)
	sd.cv.ts=sample_data(physeq.subset)
	flt=filterTaxonMatrix(otu.cv.ts,minocc=6, keepSum = TRUE, return.filtered.indices = TRUE)
	otu.cv.ts.filtered=flt$mat
	taxa.cv.ts.filtered=taxa.cv.ts[setdiff(1:nrow(taxa.cv.ts),flt$filtered.indices),]
	dummyTaxonomy=c("k__","p__","c__","o__","f__","g__","s__")
	TAXA.filt=rbind(taxa.cv.ts.filtered, dummyTaxonomy)
	ColPhylum<-c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Basidiomycota", "Ascomycota")
	TAXA.filt[,"Phylum"] <- ifelse(! TAXA.filt[,"Phylum"] %in% ColPhylum, "Other", ifelse(TAXA.filt[,"Phylum"] == "Bacteroidota", "Bacteroidetes", ifelse(TAXA.filt[,"Phylum"] == "Firmicutes", "Firmicutes", ifelse (TAXA.filt[,"Phylum"] == "Actinobacteriota", "Actinobacteria", ifelse (TAXA.filt[,"Phylum"] == "Proteobacteria","Proteobacteria", ifelse (TAXA.filt[,"Phylum"] == "Ascomycota", "Ascomycota", "Basidiomycota" ))))))
	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"	
	rownames(otu.cv.ts.filtered)[nrow(otu.cv.ts.filtered)]="SUM"
	NewOTU=otu_table(otu.cv.ts.filtered, taxa_are_rows = TRUE)	
	NewTaxa=tax_table(TAXA.filt)
	physeq.filter=phyloseq(NewOTU, NewTaxa, sd.cv.ts)
#computing bacterial and fungal co-occurrences network 
print("##### computing newtwork using mb method #####")
	spiec.mar<-spiec.easi(physeq.filter, method="mb", verbose=TRUE, nlambda=99)
	ig2.mar <- adj2igraph(getRefit(spiec.mar), vertex.attr=list(name=taxa_names(physeq.filter)))
#computing within- (Zi) and among-modules (Pi) connectivity scores
	seqdeg=degree(ig2.mar)
	Nnodes=length(seqdeg)
	Z=seqdeg
	Z[]=0
	P=Z
	fc.mar = cluster_edge_betweenness(ig2.mar,weights =NULL)
	Membership=membership(fc.mar)
	Seq=seq(1:Nnodes)
	for(i in 1:Nnodes){
		L=Membership==Membership[i]         
		neighbs=neighbors(ig2.mar, i)               
		Kis=sum(L[neighbs])
		SUM=0
		SUMsq=0	
		SUMP=0
		Miv=Seq[L]
		for(j in 1:sum(L)){
			neighbsj=neighbors(ig2.mar, Miv[j])
			Kjs=sum(L[neighbsj])
			SUM=SUM+Kjs
			SUMsq=SUMsq+Kjs^2
		}
		Z[i]=(Kis-SUM/sum(L))/sqrt(SUMsq/sum(L)-(SUM/sum(L))^2)
		if(Kis-SUM/sum(L)==0){Z[i]=0}
		for(k in 1:max(Membership)){
		Lp=Membership==k
		Kisp=sum(Lp[neighbs])
		SUMP=SUMP+(Kisp/seqdeg[i])^2}
		P[i]=1-SUMP
	}
	attr.mar=cbind(degree=seqdeg,module=Membership,Pi=P,Zi=Z)


#subset data to Sweet Ann
	physeq.subset=subset_samples(physeq.merge, Cultivar == "SweetAnne")
##filter OTUs with low occurence across samples
	otu.cv.ts=otu_table(physeq.subset)
	taxa.cv.ts=tax_table(physeq.subset)
	sd.cv.ts=sample_data(physeq.subset)
	flt=filterTaxonMatrix(otu.cv.ts,minocc=6, keepSum = TRUE, return.filtered.indices = TRUE)
	otu.cv.ts.filtered=flt$mat
	taxa.cv.ts.filtered=taxa.cv.ts[setdiff(1:nrow(taxa.cv.ts),flt$filtered.indices),]
	dummyTaxonomy=c("k__","p__","c__","o__","f__","g__","s__")
	TAXA.filt=rbind(taxa.cv.ts.filtered, dummyTaxonomy)
	ColPhylum<-c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Basidiomycota", "Ascomycota")
	TAXA.filt[,"Phylum"] <- ifelse(! TAXA.filt[,"Phylum"] %in% ColPhylum, "Other", ifelse(TAXA.filt[,"Phylum"] == "Bacteroidota", "Bacteroidetes", ifelse(TAXA.filt[,"Phylum"] == "Firmicutes", "Firmicutes", ifelse (TAXA.filt[,"Phylum"] == "Actinobacteriota", "Actinobacteria", ifelse (TAXA.filt[,"Phylum"] == "Proteobacteria","Proteobacteria", ifelse (TAXA.filt[,"Phylum"] == "Ascomycota", "Ascomycota", "Basidiomycota" ))))))
	rownames(TAXA.filt)[nrow(TAXA.filt)]="SUM"	
	rownames(otu.cv.ts.filtered)[nrow(otu.cv.ts.filtered)]="SUM"
	NewOTU=otu_table(otu.cv.ts.filtered, taxa_are_rows = TRUE)	
	NewTaxa=tax_table(TAXA.filt)
	physeq.filter=phyloseq(NewOTU, NewTaxa, sd.cv.ts)
#computing bacterial and fungal co-occurrences network 
print("##### computing newtwork using mb method #####")
	spiec.sw<-spiec.easi(physeq.filter, method="mb", verbose=TRUE, nlambda=99)
	ig2.sw <- adj2igraph(getRefit(spiec.sw), vertex.attr=list(name=taxa_names(physeq.filter)))
#computing within- (Zi) and among-modules (Pi) connectivity scores
	seqdeg=degree(ig2.sw)
	Nnodes=length(seqdeg)
	Z=seqdeg
	Z[]=0
	P=Z
	fc.sw = cluster_edge_betweenness(ig2.sw,weights =NULL)
	Membership=membership(fc.sw)
	Seq=seq(1:Nnodes)
	for(i in 1:Nnodes){
		L=Membership==Membership[i]         
		neighbs=neighbors(ig2.sw, i)               
		Kis=sum(L[neighbs])
		SUM=0
		SUMsq=0	
		SUMP=0
		Miv=Seq[L]
		for(j in 1:sum(L)){
			neighbsj=neighbors(ig2.sw, Miv[j])
			Kjs=sum(L[neighbsj])
			SUM=SUM+Kjs
			SUMsq=SUMsq+Kjs^2
		}
		Z[i]=(Kis-SUM/sum(L))/sqrt(SUMsq/sum(L)-(SUM/sum(L))^2)
		if(Kis-SUM/sum(L)==0){Z[i]=0}
		for(k in 1:max(Membership)){
		Lp=Membership==k
		Kisp=sum(Lp[neighbs])
		SUMP=SUMP+(Kisp/seqdeg[i])^2}
		P[i]=1-SUMP
	}
	attr.sw=cbind(degree=seqdeg,module=Membership,Pi=P,Zi=Z)

	DT.man=data.table(attr.man, keep.rownames=T, key="rn" )
	DT.man$cultivar="Manresa"
	
	DT.mar=data.table(attr.mar, keep.rownames=T, key="rn" )
	DT.mar$cultivar="Marquis"
	
	DT.sw=data.table(attr.sw, keep.rownames=T, key="rn" )
	DT.sw$cultivar="SweetAnn"
	
	DT=rbind(DT.man,DT.mar,DT.sw)
	DT$names=ifelse(DT$rn == "FASV66","FASV66 (Mp)","")
	DT.tax=data.table( as(tax_table(physeq.merge),"matrix"), keep.rownames=T, key="rn" )
	DT=merge(DT,DT.tax)
	
	pp=ggplot(DT, aes(Pi, Zi, color=cultivar, shape=Kingdom))+geom_point()
	
	pp=ggplot(DT, aes(Pi, Zi, color=cultivar, shape=Kingdom))
	p1=pp+geom_hline(yintercept = 2.5, linetype="longdash") + 
	geom_vline(xintercept = 0.6, linetype="longdash") +
	geom_point(size=2.5)+
	geom_text(aes(label=names))+
	theme_bw()+theme_change+
	scale_shape_manual(values=c(15,16))+
	scale_colour_manual(values=c("#4f7b95","#ffa500","#155832"))+
	ylab("Within-module connectivity Zi")+
	xlab("Among-module connectivity Pi")
	
#	pdf("figure_ZiPi.pdf",useDingbats=FALSE)
	gridExtra::grid.arrange(p1, nrow=2, ncol=2)
#	dev.off()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	


