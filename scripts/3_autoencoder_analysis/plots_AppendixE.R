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

load_physeq_data = function (){
  #loading metadata
  metadata=import_qiime_sample_data("./data/original_data/metadata_table_all_80.csv")
  
  #loading taxa information
  taxtable=read.table("./data/original_data/tax_table_all_80.csv", sep="\t")
  colnames(taxtable) = c("otuids", "Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7")
  taxtable=taxtable[2:134,]
  row.names(taxtable) = taxtable$otuids
  taxtable$otuids = NULL
  taxtable=as.matrix(taxtable)
  
  #loading original OTU table (test subset)
  original_otutable_test=read.csv("./data/original_data/data_train_autoencoder/test_subset_norm.csv", sep='\t')
  row.names(original_otutable_test) = original_otutable_test$X
  original_otutable_test$X = NULL
  colnames(original_otutable_test) = gsub ("X", "", colnames(original_otutable_test))
  
  #loading transformed OTU table (test subset) using autoencoder
  #here are loaded the results of the best autoencoder
  transformed_otutable_test=read.csv("./data/transformed_data/autoencoder/architecture5/test_autoencoder2_5_10-6-100-1000.hdf5.csv", sep='\t', header=TRUE)
  row.names(transformed_otutable_test) = transformed_otutable_test$Unnamed..0
  transformed_otutable_test$Unnamed..0 = NULL
  colnames(transformed_otutable_test) = gsub ("X", "", colnames(transformed_otutable_test))
  
  #selecting rows from metadata the correspond to the test subset
  metadata_test = metadata[row.names(metadata) %in% row.names(original_otutable_test),]
  OTU_or = otu_table(as.matrix(original_otutable_test), taxa_are_rows = FALSE)
  OTU_tr = otu_table(as.matrix(transformed_otutable_test), taxa_are_rows = FALSE)
  TAX = tax_table(taxtable)
  
  #global variables, used out of this function
  physeq_or <<- phyloseq(OTU_or, TAX)
  physeq_or <<- merge_phyloseq(physeq_or, metadata_test)
  
  physeq_tr <<- phyloseq(OTU_tr, TAX)
  physeq_tr <<- merge_phyloseq(physeq_tr, metadata_test)
}

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

plot_taxa_subset = function (data_or, data_tr, trait, trait_level, tax_rank){
  #function to plot comparative barplots between original (upper row) and transformed data by autoencoder (bottom row)
  #color according to tax_rank
  #samples are shown individually and only the ones that belong to a specific class ('trait_level') of a variable ('trait) are plotted
  
  #tax rank name definition
  ranks=c("PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES", "none")
  rank_name=ifelse(grepl("2", toString(tax_rank)), "PHYLUM", ifelse(grepl("3", toString(tax_rank)),"CLASS", ifelse(grepl("4", toString(tax_rank)), "ORDER", 
                                                                                                                   ifelse(grepl("5", toString(tax_rank)), "FAMILY", ifelse(grepl("6", toString(tax_rank)), "GENUS", ifelse(grepl("7", toString(tax_rank)),"SPECIES", "OTUs"))))))
  
  if (rank_name != "OTUs") {
    #group data according to tax rank (as a parameter)
    data_or.rank <- data_or %>% 
      tax_glom (taxrank = tax_rank) 
    data_tr.rank <- data_tr %>% 
      tax_glom (taxrank = tax_rank) 
  }
  
  else{
    data_or.rank=data_or 
    data_tr.rank=data_tr 
  }
  
  #show all the possible values that the indicated variable can take
  print(unique(sample_data(data_or.rank)[,trait]))
  
  data_or.rank_df=data_or.rank %>%
    psmelt()
  
  data_tr.rank_df=data_tr.rank %>%
    psmelt()
  
  samples_id=sort(unique(data_or.rank_df$Sample))
  
  #check that every sample is in both original and transformed data
  samples2delete=NULL
  for (i in (1: length(samples_id))){
    sample=samples_id[i]
    selected_rows=data_or.rank_df[data_or.rank_df$Sample == sample,]
    suma = sum(selected_rows$Abundance)
    if (suma==0){
      samples2delete=c(samples2delete, sample)
    }
  }
  
  #select intereseting samples
  samples_subset_or=subset(data_or.rank_df, data_or.rank_df[[trait]] == trait_level)
  samples_subset_or=samples_subset_or[samples_subset_or$Sample %not_in% samples2delete,]
  samples_subset_tr=subset(data_tr.rank_df, data_tr.rank_df[[trait]] == trait_level)
  samples_subset_tr=samples_subset_tr[samples_subset_tr$Sample %in% samples_subset_or$Sample, ]
  
  #plot
  x_text=element_text(size = 0)
  #original plot
  p1 <- ggplot(samples_subset_or, aes(x = samples_subset_or$Sample, y = Abundance, fill = samples_subset_or[,tax_rank]))+
    geom_bar(aes(), stat="identity", position="fill") + 
    guides(fill = guide_legend(title=rank_name)) +
    scale_fill_hue(c=80, l=70) +
    ylab("Relative abundance") +
    #xlab("Samples")+
    ggtitle("Original data") +
    theme(axis.title.x = x_text ,axis.text.x = element_text(angle=60,hjust=1,size='8'))
  
  #transformed plot
  p2 <- ggplot(samples_subset_tr, aes(x = samples_subset_tr$Sample, y = Abundance, fill = samples_subset_tr[,tax_rank]))+
    geom_bar(aes(), stat="identity", position="fill") + 
    guides(fill = guide_legend(title=rank_name)) +
    scale_fill_hue(c=80, l=70) +
    ylab("Relative abundance") +
    xlab("Samples")+
    ggtitle("Transformed data") +
    theme(axis.text.x = element_text(angle=60,hjust=1,size='8'))
  
  text_title=paste(rank_name, " COMPOSITION: ", trait, " = ", trait_level, sep = "")
  title=textGrob(text_title, gp=gpar(fontface="bold", fontsize=20))
  plot=grid_arrange_shared_legend(p1, p2)
  grid.arrange(plot, top = title)
}

#loading all necessary data
load_physeq_data()

#plot example
plot_taxa_subset(physeq_or, physeq_tr,  'Maize_Line',  'Popcorn', 'Rank5')


