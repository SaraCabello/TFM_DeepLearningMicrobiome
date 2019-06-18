library(phyloseq)
library(biomformat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

#load complete OTU table and mapping data from Walter et al., 2018
otu_table = import_biom("./data/original_data/otu_table_Waltersetall.biom")
metadata = import_qiime_sample_data("./data/original_data/mapping_file_Walteretall.txt")
original_data = merge_phyloseq(otu_table, metadata)

#load Walters et al., 2018 filtered data subset (717 OTUs)
percent80 = import_biom("./data/original_data/otu_table_shared80perc.biom")
percent80_OTUs=as.vector(taxa_names(percent80))

percent80_OTUs_p=gsub('^New.ReferenceOTU\\d+', "", percent80_OTUs)

OTUs=NULL
for (i in 1:length(percent80_OTUs_p)){
  if (percent80_OTUs_p[i] != ""){
    OTUs=c(OTUs, percent80_OTUs_p[i])
  }
}

#filter complete OTU table to obtain the 717 OTUs
authors_data=prune_taxa(OTUs, original_data)
#delete samples from bulk soil
data_nobulk=subset_samples(authors_data,  Description1=="rhizosphere")

#stratified multivariate sampling
samples_original=as.data.frame(sample_data(data_nobulk))
samples2save=samples_original %>%
  group_by(AGE, INBREDS, Temperature, Precipitation3Days) %>%
  sample_frac(.07)
ids2save=as.character(samples2save$X.SampleID)
data2save=subset_samples(data_nobulk, X.SampleID %in% ids2save) #2save = for predictive model
data2use=subset_samples(data_nobulk, !X.SampleID %in% ids2save) #2use = for the autoencoder

#Save biom objects to files
#pruned data = 717 OTUs and 4724 samples
data=as.data.frame(otu_table(data_nobulk))
taxas=as.data.frame(tax_table(data_nobulk))
otuids=rownames(taxas)
taxas=cbind(otuids, taxas)
metadata=as.data.frame(sample_data(data_nobulk))
biom_data=make_biom(data)
otuids=rownames(data)
data=cbind(otuids, data)
write.table(data, "./data/original_data/otu_table_all_80.csv", sep="\t", row.names = FALSE, col.names = TRUE)
write_biom(biom_data, "./data/original_data/otu_table_all_80.biom")
write.table(metadata, "./data/original_data/metadata_table_all_80.csv", sep="\t", row.names = FALSE, col.names = TRUE)
write.table(taxas, "./data/original_data/tax_table_all_80.csv", sep="\t", row.names = FALSE, col.names = TRUE)

#predictive model subset = 717 OTUs and 426 samples
data=as.data.frame(otu_table(data2save))
taxas=as.data.frame(tax_table(data2save))
otuids=rownames(taxas)
taxas=cbind(otuids, taxas)
metadata=as.data.frame(sample_data(data2save))
biom_data=make_biom(data)
otuids=rownames(data)
data=cbind(otuids, data)
write_biom(biom_data, "./data/original_data/otu_table2save_80.biom")
write.table(data, "./data/original_data/otu_table2save_80.csv", sep="\t", row.names = FALSE, col.names = TRUE)
write.table(metadata, "./data/original_data/metadata_table2save_80.csv", sep="\t", row.names = FALSE, col.names = TRUE)
write.table(taxas, "./data/original_data/tax_table2save_80.csv", sep="\t", row.names = FALSE, col.names = TRUE)

#autoencoder subset = 717 OTUs and 4298 samples
data=as.data.frame(otu_table(data2use))
taxas=as.data.frame(tax_table(data2use))
otuids=rownames(taxas)
taxas=cbind(otuids, taxas)
metadata=as.data.frame(sample_data(data2use))
biom_data=make_biom(data)
otuids=rownames(data)
data=cbind(otuids, data)
write_biom(biom_data, "./data/original_data/otu_table2use_80.biom")
write.table(data, "./data/original_data/otu_table2use_80.csv", sep="\t", row.names = FALSE, col.names = TRUE)
write.table(metadata, "./data/original_data/metadata_table2use_80.csv", sep="\t", row.names = FALSE, col.names = TRUE)
write.table(taxas, "./data/original_data/tax_table2use_80.csv", sep="\t", row.names = FALSE, col.names = TRUE)

#Varify mapping variables distribution
data_original=read.table("./data/original_data/metadata_table_all_80.csv", sep='\t', header=TRUE)
p1 <- ggplot(data_original, aes(x=as.numeric(data_original$AGE))) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), color="black", fill="grey") +
  labs(x = "Age")  + theme(axis.title.y=element_blank())
p2 <- ggplot(data_original) + 
  geom_bar(aes(x=data_original$INBREDS, y= (..count..)/sum(..count..)), color="black", fill="grey") + 
  labs (x = "Inbred", y = 'Frequency') + 
  theme(axis.text.x=element_text(angle=80,hjust=1,size='5'))
p3 <- ggplot(data_original, aes(x=as.numeric(data_original$Temperature))) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), color="black", fill="grey")+
  labs (x = "Temperature")  + theme(axis.title.y=element_blank())
p4 <- ggplot(data_original, aes(x=as.numeric(data_original$Precipitation3Days))) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), color="black", fill="grey") +
  labs(x = "Precipitation3Days")  + theme(axis.title.y=element_blank())
title1=textGrob("A) ORIGINAL SUBSET", gp=gpar(fontface="bold", fontsize=15))
plot_original=grid.arrange(p2, p1, p3, p4, nrow=1, top = title1 )

data2save=read.table("./data/original_data/metadata_table2save_80.csv", sep='\t', header=TRUE)
p5 <- ggplot(data2save, aes(x=as.numeric(data2save$AGE))) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), color="black", fill="grey") +
  labs(x = "Age")  + theme(axis.title.y=element_blank())
p6 <- ggplot(data2save) + 
  geom_bar(aes(x=data2save$INBREDS, y= (..count..)/sum(..count..)), color="black", fill="grey") + 
  labs (x = "Inbred", y = 'Frequency') + 
  theme(axis.text.x=element_text(angle=80,hjust=1,size='5'))
p7 <- ggplot(data2save, aes(x=as.numeric(data2save$Temperature))) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), color="black", fill="grey")+
  labs (x = "Temperature")  + theme(axis.title.y=element_blank())
p8 <- ggplot(data2save, aes(x=as.numeric(data2save$Precipitation3Days))) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), color="black", fill="grey") +
  labs(x = "Precipitation3Days")  + theme(axis.title.y=element_blank())
title2=textGrob("C) PREDICTIVE MODEL SUBSET", gp=gpar(fontface="bold", fontsize=15)) 
plot2save=grid.arrange(p6, p5, p7, p8, nrow=1, top = title2 )

data2use=read.table("./data/original_data/metadata_table2use_80.csv", sep='\t', header=TRUE)
p9 <- ggplot(data2use, aes(x=as.numeric(data2use$AGE))) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), color="black", fill="grey") +
  labs(x = "Age")  + theme(axis.title.y=element_blank())
p10 <- ggplot(data2use) + 
  geom_bar(aes(x=data2use$INBREDS, y= (..count..)/sum(..count..)), color="black", fill="grey") + 
  labs (x = "Inbred", y = 'Frequency') + 
  theme(axis.text.x=element_text(angle=80,hjust=1,size='5'))
p11 <- ggplot(data2use, aes(x=as.numeric(data2use$Temperature))) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), color="black", fill="grey")+
  labs (x = "Temperature")  + theme(axis.title.y=element_blank())
p12 <- ggplot(data2use, aes(x=as.numeric(data2use$Precipitation3Days))) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), color="black", fill="grey") +
  labs(x = "Precipitation3Days")  + theme(axis.title.y=element_blank())
title2=textGrob("B) AUTOENCODER SUBSET", gp=gpar(fontface="bold", fontsize=15)) 
plot2use=grid.arrange(p10, p9, p11, p12, nrow=1, top = title2 )
grid.arrange(plot_original, plot2use, plot2save)
