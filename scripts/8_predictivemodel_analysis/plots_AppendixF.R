library(phyloseq)
library(biomformat)
library(ggplot2)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)
library(hydroGOF)
library(data.table)
`%not_in%` <- purrr::negate(`%in%`)

grid_arrange_shared_legend <- function(...) {
  #pre-built function (public online) to use only one color legend in a composed plot
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

load_physeq_data = function (){
  #loading metadata
  metadata=import_qiime_sample_data("./data/original_data/metadata_table_all_80.csv")
  
  #load taxonomic information
  taxtable=read.table("./data/original_data/tax_table_all_80.csv", sep="\t")
  colnames(taxtable) = c("otuids", "Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7")
  taxtable=taxtable[2:134,]
  row.names(taxtable) = taxtable$otuids
  taxtable$otuids = NULL
  taxtable=as.matrix(taxtable)
  
  #load original otutable from test subset of predictive model training stage
  original_otutable_test=read.csv("./data/original_data/data_train_predictivemodel/2otu_table2save_80_ytest20_norm.csv", sep='\t')
  row.names(original_otutable_test) = original_otutable_test$X.SampleID
  original_otutable_test$X.SampleID = NULL
  original_otutable_test$X.SampleID.1 = NULL
  colnames(original_otutable_test) = gsub ("X", "", colnames(original_otutable_test))
  
  #load transformed data using the best predictive model: 4-4-20-60-100-1000
  transformed_otutable_test=read.csv("./data/transformed_data/predictive_model/best_model/bestmodel_test.csv", sep='\t', header=TRUE)
  row.names(transformed_otutable_test) = transformed_otutable_test$X
  transformed_otutable_test$X = NULL
  colnames(transformed_otutable_test) = gsub ("X", "", colnames(transformed_otutable_test))
  
  #select metadata of test subset
  metadata_test = metadata[row.names(metadata) %in% row.names(original_otutable_test),]
  OTU_or = otu_table(as.matrix(original_otutable_test), taxa_are_rows = FALSE)
  OTU_tr = otu_table(as.matrix(transformed_otutable_test), taxa_are_rows = FALSE)
  TAX = tax_table(taxtable)
  
  #global variables
  physeq_or <<- phyloseq(OTU_or, TAX)
  physeq_or <<- merge_phyloseq(physeq_or, metadata_test)
  
  physeq_tr <<- phyloseq(OTU_tr, TAX)
  physeq_tr <<- merge_phyloseq(physeq_tr, metadata)
}

plot_taxa_samples= function (data_or, data_tr, tax_rank){
  #function that plots OTU composition of individual samples from test subset
  #comparing original vs. transformed using predictive model
  #the bars are colored according to the taxonomic level indicated as 'tax_rank0
  
  #looking for the name of the taxonomic level
  ranks=c("PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES", "none")
  rank_name=ifelse(grepl("2", toString(tax_rank)), "PHYLUM", ifelse(grepl("3", toString(tax_rank)),"CLASS", ifelse(grepl("4", toString(tax_rank)), "ORDER", 
                                                                                                                   ifelse(grepl("5", toString(tax_rank)), "FAMILY", ifelse(grepl("6", toString(tax_rank)), "GENUS", ifelse(grepl("7", toString(tax_rank)),"SPECIES", "OTUs"))))))
  
  if (rank_name != "OTUs") {
    #group data according to its taxonomic level (at the level of tax_rank)
    data_or.rank <- data_or %>% 
      tax_glom (taxrank = tax_rank) 
    data_tr.rank <- data_tr %>% 
      tax_glom (taxrank = tax_rank) 
  }
  
  else{
    data_or.rank=data_or 
    data_tr.rank=data_tr 
  }
  
  data_or.rank_df=data_or.rank %>%
    psmelt()
  
  data_tr.rank_df=data_tr.rank %>%
    psmelt()
  
  #check that every sample is in both original and transformed data
  samples_id=sort(unique(data_or.rank_df$Sample))
  samples2delete=NULL
  for (i in (1: length(samples_id))){
    sample=samples_id[i]
    selected_rows=data_or.rank_df[data_or.rank_df$Sample == sample,]
    suma = sum(selected_rows$Abundance)
    if (suma==0){
      samples2delete=c(samples2delete, sample)
    }
  }
  
  samples_subset_or=data_or.rank_df
  samples_subset_or=samples_subset_or[samples_subset_or$Sample %not_in% samples2delete,]
  samples_subset_tr=data_tr.rank_df
  samples_subset_tr=samples_subset_tr[samples_subset_tr$Sample %in% samples_subset_or$Sample, ]
  
  #plot
  x_text=element_text(size = 0)
  #original data
  p1 <- ggplot(samples_subset_or, aes(x = samples_subset_or$Sample, y = Abundance, fill = samples_subset_or[,tax_rank]))+
    geom_bar(aes(), stat="identity", position="fill") + 
    guides(fill = guide_legend(title=rank_name)) +
    scale_fill_hue(c=80, l=70) +
    ylab("Relative abundance") +
    #xlab("Samples")+
    ggtitle("Original data") +
    theme(axis.title.x = x_text ,axis.text.x = x_text)
  
  #transformed data
  p2 <- ggplot(samples_subset_tr, aes(x = samples_subset_tr$Sample, y = Abundance, fill = samples_subset_tr[,tax_rank]))+
    geom_bar(aes(), stat="identity", position="fill") + 
    guides(fill = guide_legend(title=rank_name)) +
    scale_fill_hue(c=80, l=70) +
    ylab("Relative abundance") +
    xlab("Samples")+
    ggtitle("Transformed data") +
    theme(axis.text.x = x_text)
  
  text_title=paste(rank_name, " COMPOSITION" )
  title=textGrob(text_title, gp=gpar(fontface="bold", fontsize=20))
  plot=grid_arrange_shared_legend(p1, p2)
  grid.arrange(plot, top = title)
}

load_physeq_data()

#example
plot_taxa_samples(physeq_or, physeq_tr, 'Rank3')


