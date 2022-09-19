#require packages for the analysis the analysis
pkg=c( "ggplot2", "phyloseq", "data.table", "reshape2", "seqtime" ,"SpiecEasi", "huge")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])
set.seed(131)
Palette.phylum <-c("#EF5656","#47B3DA","#F7A415","#88888a","#2BB065")

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
break()
#computing bacterial and fungal co-occurrences network 
print("##### computing newtwork using mb method #####")
	spiec.out.mb<- spiec.easi(physeq.filter, method="mb", verbose=TRUE, nlambda=99)
	ig2.mb <- adj2igraph(getRefit(spiec.out.mb),  vertex.attr=list(name=taxa_names(physeq.filter)))
	p1<-plot_network(ig2.mb, physeq.filter, type='taxa', color="Phylum", shape="Kingdom", label=NULL)

#	seqdeg=degree(igraph.group1)
#	Nnodes=length(seqdeg)
#	Z=seqdeg
#	Z[]=0
#	P=Z
#	print(paste0("edge_betweenness =",modularity(fgc)))


print("##### network properties #####")
Mat.mb=as.matrix(symBeta(getOptBeta(spiec.out.mb)))
print(paste0("Nbr_nodes=", vcount(ig2.mb)))
print(paste0("Nbr_edges=", ecount(ig2.mb)))
print(paste0("Nbr_total_edges=",length(Mat.mb[Mat.mb!=0])/2))
print(paste0("Nbr_pos_edges=",length(Mat.mb[Mat.mb>0])/2))
print(paste0("Nbr_neg_edges=",length(Mat.mb[Mat.mb<0])/2))
print(paste0("clus. coef. glo.=",transitivity(ig2.mb, type = "global")))
print(paste0("clus. coef. ave.=",transitivity(ig2.mb, type = "average")))
print(paste0("graph_density=", graph.density(ig2.mb)))
print(paste0("graph_diam=", diameter(ig2.mb)))
print(paste0("ave_path_len=", average.path.length(ig2.mb)))

	otu.ids=colnames(spiec.out.mb$est$data)
	edges=E(ig2.mb)
	edge.colors=c()
	for(e.index in 1:length(edges)){
		adj.nodes=ends(ig2.mb,edges[e.index])
		xindex=which(otu.ids==adj.nodes[1])
		yindex=which(otu.ids==adj.nodes[2])
		beta=Mat.mb[xindex,yindex]
		if(beta>0){
			edge.colors=append(edge.colors,"darkgreen")
		}else if(beta<0){
			edge.colors=append(edge.colors,"darkred")
		}
	}
	E(ig2.mb)$color=edge.colors
	graph.mb=ig2.mb
	nodenames=V(graph.mb)$name
	V(graph.mb)$name=getTaxonomy(nodenames,TAXA.filt, level="phylum", useRownames=TRUE)
	E(graph.mb)$arrow.size=5
	V(graph.mb)$color="white"
	V(graph.mb)$frame.color="black"

	V(graph.mb)$color <- ifelse(V(graph.mb)$name == "Other", "#000000", ifelse(V(graph.mb)$name == "Actinobacteria","#EF5656",ifelse(V(graph.mb)$name=="Bacteroidetes","#47B3DA",ifelse(V(graph.mb)$name =="Firmicutes", "#F7A415",ifelse(V(graph.mb)$name =="Proteobacteria","#2BB065", ifelse(V(graph.mb)$name =="Ascomycota", "#4e8e8e", "#8e8e4e"))))))

	V(graph.mb)$name=getTaxonomy(nodenames,TAXA.filt, level="domain", useRownames=TRUE)
	V(graph.mb)$shape <- ifelse(V(graph.mb)$name == "Bacteria", "circle" , "square")

	V(graph.mb)$name=getTaxonomy(nodenames,TAXA.filt, level="genus" ,useRownames=TRUE)

#	pdf("network_ASV_mb_rhizo.pdf", useDingbats=FALSE)
	plot.igraph(graph.mb, vertex.size=3,layout=layout.fruchterman.reingold,vertex.label=NA) 
#	dev.off()
	
# 	save network metrics	
	metrics.mb <- data.frame ( 
	deg.mb=degree(ig2.mb), 
	bet.mb=betweenness(ig2.mb), 
	clo.mb=closeness(ig2.mb), 
	eig.mb=evcent(ig2.mb) 
	)
	DT_mb=data.table(metrics.mb, keep.rownames=T, key="rn")
	DT_mb=DT_mb[,1:5]
	DT_hs=data.table(data.frame(hubScore=hub_score(ig2.mb)$vector), keep.rownames=T, key="rn")
	DT_taxa=data.table(TAXA.filt, keep.rownames=T, key="rn")

	DT.m=DT_mb[DT_hs,][DT_taxa,]



	TAXA.filt[,7]<-rownames(TAXA.filt)
	V(graph.mb)$name=getTaxonomy(nodenames,TAXA.filt, level="species" ,useRownames=TRUE)
	Selec <- V(graph.mb)[name =="FASV66"]
#	Selec <- V(graph.mb)[name %in% c("FASV66","FASV310","FASV675","FASV917","FASV3800")]

	selecV <- ego(graph.mb, order=1, nodes = Selec, mode = "all", mindist = 0)
	selecG <- induced_subgraph(graph.mb,unlist(selecV))

#	pdf("network_Mc_mb_rhizo.pdf", useDingbats=FALSE)
	plot.igraph(selecG, vertex.size=6, vertex.label.dist=0.5, layout=layout.fruchterman.reingold, vertex.label.color="black", vertex.label.cex=0.6)
#	dev.off()

print(paste0("Nbr_nodes=", vcount(selecG)))
print(paste0("Nbr_edges=", ecount(selecG)))
print(paste0("clus. coef. glo.=",transitivity(selecG, type = "global")))
print(paste0("clus. coef. ave.=",transitivity(selecG, type = "average")))
print(paste0("graph_density=", graph.density(selecG)))
print(paste0("graph_diam=", diameter(selecG)))
print(paste0("ave_path_len=", average.path.length(selecG)))

# 	save subnetwork metrics	
	metrics.sub <- data.frame ( 
	deg.sub=degree(selecG), 
	bet.sub=betweenness(selecG), 
	clo.sub=closeness(selecG), 
	eig.sub=evcent(selecG) 
	)
	DT_sub=data.table(metrics.sub, keep.rownames=T, key="rn")
	DT_sub=DT_sub[,1:5]
	DT_hs=data.table(data.frame(hubScore=hub_score(selecG)$vector), keep.rownames=T, key="rn")
	
	DT_sub=DT_sub[DT_hs,][DT_taxa,]
	DT_sub=DT_sub[DT_sub$deg.sub != "NA"]
	


	DT_ra=data.table(melt(NewOTU[DT_sub$rn,]),keep.rownames=TRUE, key=c("Var1","Var2"))
	
	DT_sd=data.table(as(sample_data(physeq.subset),"data.frame"), keep.rownames=TRUE, key=c("rn"))

	DT_m1=merge(DT_ra, DT_sd, by.x="Var2", by.y="rn")
	DT_m2=merge(DT_m1, DT_sub, by.x="Var1", by.y="rn")
	
	DT_net=data.table(as_data_frame(selecG),keep.rownames=TRUE, key=c("rn"))

	DT_neg=DT_net[DT_net$color == "darkred"]
	DT_neg=DT_neg[DT_neg$from == "FASV66" | DT_neg$to == "FASV66"]

	DT_pos=DT_net[DT_net$color != "darkred"]
	DT_pos=DT_pos[DT_pos$from == "FASV66" | DT_pos$to == "FASV66"]

	DT_66=DT_net[DT_net$from == "FASV66" | DT_net$to == "FASV66"]

	DT_66$col=ifelse(DT_66$to == "FASV66", DT_66$from, DT_66$to)

	DT_m3=DT_m2[DT_m2$Var1 %in% DT_66$col]

	DT_m3=merge(DT_m3, DT_66, by.x="Var1",by.y="col")

	write.table(DT.m, "sweetann_metrics_full.txt", sep="\t")
	write.table(DT_sub, "sweetann_metrics_sub.txt", sep="\t")
	write.table(as_data_frame(graph.mb), "sweetann_network_full.txt", sep="\t")
	write.table(as_data_frame(selecG), "sweetann_network_sub", sep="\t")
	write.table(DT_m3, "sweetann_sub_relative.txt", sep="\t")

#	write.table(otu_table(NewOTU), "sweetann_Species_OTU_tab.txt", sep="\t")













