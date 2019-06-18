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

plot_taxa_group = function(data_or, data_tr, trait, tax_rank, w){
  #function to plot comparative barplots between original (upper row) and transformed data by autoencoder (bottom row)
  #color according to tax_rank
  #samples are group according to the different values found in a particular variable (indicated in 'trait' parameter)
  #w indicates the width of the bars, if w=NULL --> deafult value
  
  #rank name definition
  ranks=c("KINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES", "none")
  rank_name=ifelse(grepl("1", toString(tax_rank)), "KINGDOM", ifelse(grepl("2", toString(tax_rank)), "PHYLUM", ifelse(grepl("3", toString(tax_rank)),"CLASS", ifelse(grepl("4", toString(tax_rank)), "ORDER", 
                                                                                                                                                                     ifelse(grepl("5", toString(tax_rank)), "FAMILY", ifelse(grepl("6", toString(tax_rank)), "GENUS", ifelse(grepl("7", toString(tax_rank)),"SPECIES", "OTUs")))))))
  if (rank_name != "OTUs") {
    #group data according to tax rank (as a parameter)
    data_or.rank <- data_or %>% 
      tax_glom (taxrank = tax_rank) %>%
      transform_sample_counts(function(x) x / sum(x) ) %>%
      psmelt() 
    
    data_tr.rank <- data_tr %>% 
      tax_glom (taxrank = tax_rank) %>%
      transform_sample_counts(function(x) x / sum(x) ) %>%
      psmelt()
  }
  else{
    data_or.rank=data_or
    data_tr.rank=data_tr
  }
  
  #plot using original data
  p1 <- ggplot(data_or.rank, aes(x = data_or.rank[,trait], y = Abundance, fill = data_or.rank[,tax_rank]))+
    geom_bar(aes(), stat="identity", position="fill", width=w) + 
    guides(fill = guide_legend(title=rank_name)) +
    scale_fill_hue(c=80, l=70) +
    theme(axis.title.y=element_blank()) +
    theme(axis.title.x=element_blank()) +
    ggtitle("Original data")
  
  #plot with transformed data
  p2 <- ggplot(data_tr.rank, aes(x = data_tr.rank[,trait], y = Abundance, fill = data_tr.rank[,tax_rank]))+
    geom_bar(aes(), stat="identity", position="fill", width=w) + 
    guides(fill = guide_legend(title=rank_name)) +
    scale_fill_hue(c=80, l=70) +
    theme(axis.title.y =element_blank()) +
    xlab(trait)+
    ggtitle("Transformed data")
  
  text_title=paste(rank_name, "COMPOSITION", sep = " ")
  title=textGrob(text_title, gp=gpar(fontface="bold", fontsize=20))
  plot=grid_arrange_shared_legend(p1, p2)
  grid.arrange(plot, top = title)
}

#loading all necessary data
load_physeq_data()

#plot example
plot_taxa_group(physeq_or, physeq_tr, 'Maize_Line', 'Rank3', w=NULL)
