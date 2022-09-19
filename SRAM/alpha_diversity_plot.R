#***require packages for the analysis the analysis
pkg=c("ggplot2", "phyloseq", "ape", "picante", "data.table","PMCMR")
lapply(pkg, library, character.only = TRUE)

rm(list = ls()[grep("*.*", ls())])

alphaInd = c("Shannon", "Observed")


color_palette1<-c("#4f7b95","#ffa500","#000000","#155832")

print(" R scripts to compute alpha-diversity measures displayed in Fig 1A, Fig 1B, Supp Fig 1A & Supp Fig 1B ");
set.seed(123)
theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
#	legend.position="none",
	axis.text.x = element_text(angle=90, vjust=1, size=11, color="black"),
	)
theme_new <- theme (
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	axis.title.y=element_blank(),
	axis.title.x=element_blank(),
#	legend.position="none",
	axis.text.x = element_blank(),
	)
	
Order1=c("Bulk.none.none","Bulk.none.Mac","Rhizo.Manresa.none","Rhizo.Manresa.Mac","Rhizo.Marquis.none",
"Rhizo.Marquis.Mac","Rhizo.SweetAnne.none","Rhizo.SweetAnne.Mac","Root.Manresa.none","Root.Manresa.Mac",
"Root.Marquis.none","Root.Marquis.Mac","Root.SweetAnne.none","Root.SweetAnne.Mac")

##upload and prepare phyloseq objects bacteria data
	tmp<-get(load("/home/mahassani/Biodata/NewHaven/SRAM/phyloseq_nochim_silva.RData"))
	OTU=tmp@otu_table
	TAXA=tmp@tax_table
	SD=tmp@sam_data
	TREE=tmp@phy_tree
	physeq=phyloseq(OTU, TAXA, SD,TREE) 
	
	physeq.BA=subset_taxa( physeq, Kingdom %in% c("Archaea","Bacteria"))
	physeq.BA=subset_taxa( physeq.BA, Genus != "Mitochondria")
	physeq1=subset_samples( physeq.BA, Field == "Field1" )
		
#Rarefication to even depth based on smallest sample size
	print(paste0("minimum sequencing depth for bacteria/archaea 1000"))
	physeq.rrf=rarefy_even_depth(physeq1, 1000, replace=TRUE, rngseed = 131)

	sample_data(physeq.rrf)$Concat<-paste0(sample_data(physeq.rrf)$SampleType,".",
	sample_data(physeq.rrf)$Cultivar,".",sample_data(physeq.rrf)$Pathogen)	
	
#Ploting species richness	
	ppBac=plot_richness(physeq.rrf,"Concat","Cultivar" ,measures=alphaInd)
	ppBac$data$Concat <- ordered(ppBac$data$Concat, levels=Order1 )
	ppBac$layers <- ppBac$layers[-1]
	pBac=ppBac+geom_boxplot(data=ppBac$data, aes(x=Concat, y=value, color=Cultivar))+
	theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	scale_colour_manual(values=color_palette1)+theme_new
	
	
	p1 = ggplot(data=pBac$data[pBac$data$variable=="Observed",], aes(x=Concat, y=value))+
	geom_boxplot(aes(color=Cultivar), outlier.shape = 1, outlier.size= 1)+
	geom_point(size=.5)+theme_bw()+ theme_new +
	scale_color_manual(values=color_palette1) + 
	facet_wrap(~variable) 
	ylab("Observed ASVs")
#testing significance in observed CTUs
kruskal.test(data=p1$data[p1$data$SampleType == "Bulk",], value ~ Pathogen)
kruskal.test(data=p1$data[p1$data$SampleType == "Rhizo" & p1$data$Cultivar == "Manresa",], value ~ Pathogen)
kruskal.test(data=p1$data[p1$data$SampleType == "Rhizo" & p1$data$Cultivar == "Marquis",], value ~ Pathogen)
kruskal.test(data=p1$data[p1$data$SampleType == "Rhizo" & p1$data$Cultivar == "SweetAnne",], value ~ Pathogen)
kruskal.test(data=p1$data[p1$data$SampleType == "Root" & p1$data$Cultivar == "Manresa",], value ~ Pathogen)
kruskal.test(data=p1$data[p1$data$SampleType == "Root" & p1$data$Cultivar == "Marquis",], value ~ Pathogen)
kruskal.test(data=p1$data[p1$data$SampleType == "Root" & p1$data$Cultivar == "SweetAnne",], value ~ Pathogen)


	p2 = ggplot(data=pBac$data[pBac$data$variable=="Shannon",], aes(x=Concat, y=value))+
	geom_boxplot(aes(color=Cultivar), outlier.shape = 1, outlier.size= 1)+
	geom_point(size=.5)+theme_bw()+ theme_new +
	scale_color_manual(values=color_palette1) + 
	facet_wrap(~variable)+
	ylab("Shannon index")
#testing significance in Shannon index
kruskal.test(data=p2$data[p2$data$SampleType == "Bulk",], value ~ Pathogen)
kruskal.test(data=p2$data[p2$data$SampleType == "Rhizo" & p2$data$Cultivar == "Manresa",], value ~ Pathogen)
kruskal.test(data=p2$data[p2$data$SampleType == "Rhizo" & p2$data$Cultivar == "Marquis",], value ~ Pathogen)
kruskal.test(data=p2$data[p2$data$SampleType == "Rhizo" & p2$data$Cultivar == "SweetAnne",], value ~ Pathogen)
kruskal.test(data=p2$data[p2$data$SampleType == "Root" & p2$data$Cultivar == "Manresa",], value ~ Pathogen)
kruskal.test(data=p2$data[p2$data$SampleType == "Root" & p2$data$Cultivar == "Marquis",], value ~ Pathogen)
kruskal.test(data=p2$data[p2$data$SampleType == "Root" & p2$data$Cultivar == "SweetAnne",], value ~ Pathogen)

      
##upload and prepare phyloseq objects fungi data
	tmp<-get(load("/home/mahassani/Biodata/NewHaven/SRAM/phyloseq_nochim_unite.RData"))
	OTU=tmp@otu_table
	TAXA=tmp@tax_table
	SD=tmp@sam_data
	physeq=phyloseq(OTU, TAXA, SD) 

#Rarefication to even depth based on smallest sample size
	physeq.Fun=subset_taxa( physeq, Kingdom == "k__Fungi")
	physeq1=subset_samples( physeq.Fun, Field == "Field1" )	
	print(paste0("minimum sequencing depth for fungi 1000"))
	physeq.rrf=rarefy_even_depth(physeq1, 1000, replace=TRUE, rngseed = 131)

	sample_data(physeq.rrf)$Concat<-paste0(sample_data(physeq.rrf)$SampleType,".",
	sample_data(physeq.rrf)$Cultivar,".",sample_data(physeq.rrf)$Pathogen)	
	
#Ploting species richness	
	
	ppFun=plot_richness(physeq.rrf,"Concat","Cultivar" ,measures=alphaInd)
	ppFun$data$Concat <- ordered(ppFun$data$Concat, levels=Order1 )
	ppFun$layers <- ppFun$layers[-1]
	pFun=ppFun+geom_boxplot(data=ppFun$data, aes(x=Concat, y=value, color=Cultivar))+
	theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
	scale_colour_manual(values=color_palette1)+theme_new
	
	p3 = ggplot(data=pFun$data[pFun$data$variable=="Observed",], aes(x=Concat, y=value))+
	geom_boxplot(aes(color=Cultivar), outlier.shape = 1, outlier.size= 1)+
	geom_point(size=.5)+theme_bw()+ theme_new +
	scale_color_manual(values=color_palette1) + 
	facet_wrap(~variable) 
	ylab("Observed ASVs")
#testing significance in observed CTUs
kruskal.test(data=p3$data[p3$data$SampleType == "Bulk",], value ~ Pathogen)
kruskal.test(data=p3$data[p3$data$SampleType == "Rhizo" & p3$data$Cultivar == "Manresa",], value ~ Pathogen)
kruskal.test(data=p3$data[p3$data$SampleType == "Rhizo" & p3$data$Cultivar == "Marquis",], value ~ Pathogen)
kruskal.test(data=p3$data[p3$data$SampleType == "Rhizo" & p3$data$Cultivar == "SweetAnne",], value ~ Pathogen)
kruskal.test(data=p3$data[p3$data$SampleType == "Root" & p3$data$Cultivar == "Manresa",], value ~ Pathogen)
kruskal.test(data=p3$data[p3$data$SampleType == "Root" & p3$data$Cultivar == "Marquis",], value ~ Pathogen)
kruskal.test(data=p3$data[p3$data$SampleType == "Root" & p3$data$Cultivar == "SweetAnne",], value ~ Pathogen)


	p4 = ggplot(data=pFun$data[pFun$data$variable=="Shannon",], aes(x=Concat, y=value))+
	geom_boxplot(aes(color=Cultivar), outlier.shape = 1, outlier.size= 1)+
	geom_point(size=.5)+theme_bw()+ theme_new +
	scale_color_manual(values=color_palette1) + 
	facet_wrap(~variable)+
	ylab("Shannon index")

#testing significance in Shannon index
kruskal.test(data=p4$data[p4$data$SampleType == "Bulk",], value ~ Pathogen)
kruskal.test(data=p4$data[p4$data$SampleType == "Rhizo" & p4$data$Cultivar == "Manresa",], value ~ Pathogen)
kruskal.test(data=p4$data[p4$data$SampleType == "Rhizo" & p4$data$Cultivar == "Marquis",], value ~ Pathogen)
kruskal.test(data=p4$data[p4$data$SampleType == "Rhizo" & p4$data$Cultivar == "SweetAnne",], value ~ Pathogen)
kruskal.test(data=p4$data[p4$data$SampleType == "Root" & p4$data$Cultivar == "Manresa",], value ~ Pathogen)
kruskal.test(data=p4$data[p4$data$SampleType == "Root" & p4$data$Cultivar == "Marquis",], value ~ Pathogen)
kruskal.test(data=p4$data[p4$data$SampleType == "Root" & p4$data$Cultivar == "SweetAnne",], value ~ Pathogen)
	

#	pdf('figure_1.pdf', useDingbats=FALSE)
	gridExtra::grid.arrange(p2,p4,p1, p3, nrow=2, ncol=2)
#	dev.off()
























       

