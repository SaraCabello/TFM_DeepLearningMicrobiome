library(ggplot2)

better_than_default = function(df){
  binary_data = df[,3]
  better = NULL
  for (i in (1:length(binary_data))){
    if (binary_data[i]==1){
      better = c(better, row.names(df)[i])
    }
  }
  return(better)
}

worse_than_default = function(df){
  binary_data = df[,3]
  worse = NULL
  for (i in (1:length(binary_data))){
    if (binary_data[i]==0){
      worse = c(worse, row.names(df)[i])
    }
  }
  return(worse)
}

otus_mse=read.csv("./data/transformed_data/autoencoder/results_bestautoencoder/mse_otus_comparison.csv", sep='\t')
row.names(otus_mse) = otus_mse$otu_id
otus_mse$otu_id = NULL

otus_smape=read.csv("./data/transformed_data/autoencoder/results_bestautoencoder/smape_otus_comparison.csv", sep='\t')
row.names(otus_smape) = otus_smape$otu_id
otus_smape$otu_id = NULL

better_otus_mse = better_than_default(otus_mse)
better_otus_smape = better_than_default(otus_smape)
worse_otus_mse = worse_than_default(otus_mse)
worse_otus_smape = worse_than_default(otus_smape)

taxtable=read.table("./data/original_data/tax_table_all_80.csv", sep="\t")
colnames(taxtable) = c("otuids", "Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7")
taxtable=taxtable[2:718,]
row.names(taxtable) = taxtable$otuids
taxtable$otuids = NULL

taxa_worse_mse = taxtable[row.names(taxtable) %in% worse_otus_mse,]

data=as.data.frame(table(taxa_worse_mse$Rank2))

ggplot(data, aes(x=Var1, y=Freq)) + 
  geom_bar(stat='identity', fill='grey', color ='black') + 
  labs (x = "Phylum", y = 'Count', title='Phylum distribution of the worst predicted OTUs (MSE criterion)') + 
  theme(axis.text.x=element_text(angle=60,hjust=1,size='12'),
        plot.title = element_text(size=20, face="bold"))

taxa_worse_smape = taxtable[row.names(taxtable) %in% worse_otus_smape,]
data=as.data.frame(table(taxa_worse_smape$Rank2))
ggplot(data, aes(x=Var1, y=Freq)) + 
  geom_bar(stat='identity', fill='grey', color ='black') + 
  labs (x = "Phylum", y = 'Count', title='Phylum distribution of the worst predicted OTUs (SMAPE criterion)') + 
  theme(axis.text.x=element_text(angle=60,hjust=1,size='12'),
        plot.title = element_text(size=20, face="bold"))

