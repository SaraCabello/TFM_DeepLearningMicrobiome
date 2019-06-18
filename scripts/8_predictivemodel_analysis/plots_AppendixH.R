library(dplyr)
library(ggplot2)
library(scales)
`%not_in%` <- purrr::negate(`%in%`)

plot_best_predictivemodel = function (data, error_type, tax_rank){
  #function that plots a pie chart that represent taxonomic classification of OTUs better predicted 
  #using predictive model than the deafult predictor not included in the ones better predicted than 
  #defautl by the direct predictor
  
  if (error_type == 'MSE') {
    data_plot = data %>%
      select (id, tr1_mse, tax_rank)
    colnames(data_plot) = c('id', 'error', 'taxa')
    suffix ='_mse'
  } else if (error_type == 'SMAPE'){
    data_plot = data %>%
      select (id, tr1_smape, tax_rank)
    colnames(data_plot) = c('id', 'error', 'taxa')
    suffix='_smape'
  }
  
  ranks=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "none")
  rank_name=ifelse(grepl("1", toString(tax_rank)), "Kingdom", ifelse(grepl("2", toString(tax_rank)), "Phylum", ifelse(grepl("3", toString(tax_rank)),"Class", ifelse(grepl("4", toString(tax_rank)), "Order", 
                    ifelse(grepl("5", toString(tax_rank)), "Family", ifelse(grepl("6", toString(tax_rank)), "Genus", ifelse(grepl("7", toString(tax_rank)),"Species", "OTUs")))))))
  
  path='./data/transformed_data/predictive_model/best_model/'
  file_better = paste('OTUs_betterwithcode', suffix, '.csv', sep='')
  file_overlap = paste('OTUs_overlap', suffix, '.csv', sep='')
  better = read.table(paste(path, file_better, sep=''), sep='\t', header = TRUE)
  overlap = read.table(paste(path, file_overlap, sep=''), sep='\t', header = TRUE)
  data_plot_better = data_plot[data_plot$id %in% better$id,]
  data_plot_better = data_plot_better[data_plot_better$id %not_in% overlap$id,]
  table = as.data.frame(table(data_plot_better$taxa))
  suma = sum(table$Freq)
  table$Freq = (table$Freq/suma) * 100
  #plot
  p1 <- ggplot(table, aes(x='', y=Freq, fill=Var1))+
    geom_bar(width=1, stat='identity') +
    coord_polar('y', start=0) +
    guides(fill = guide_legend(title=rank_name)) +

    theme(axis.title=element_blank(), plot.title = element_text(size=15, face="bold"))   +
    labs (title=paste(rank_name, ' distribution of best predicted OTUs using the predictive model', sep=''),
          subtitle = paste('According to ', error_type, ' criterion', sep='')) 
    
  p1
  }

#load data
path='./data/transformed_data/predictive_model/best_model/'
file='/table2save_otustest.csv'
data = read.table(paste(path, file, sep=''), sep='\t', header = TRUE)

#example
plot_best_predictivemodel(data, 'SMAPE', 'Rank2')

#all OTUs: phylum level
table = as.data.frame(table(data$Rank2))
suma = sum(table$Freq)
table$Freq = (table$Freq/suma) * 100
p1 <- ggplot(table, aes(x='', y=Freq, fill=Var1))+
  geom_bar(width=1, stat='identity') +
  coord_polar('y', start=0) +
  guides(fill = guide_legend(title='Phylum')) +
  theme(axis.title=element_blank(), plot.title = element_text(size=15, face="bold"))   +
  labs (title='Phylum distribution of all OTUs') 

p1