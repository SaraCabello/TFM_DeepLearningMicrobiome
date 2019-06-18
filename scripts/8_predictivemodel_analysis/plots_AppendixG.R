library(ggplot2)
library(dplyr)

plot_samples_predictivemodel=function(data, error_type, color_variable, w){
  #function that plot the mean and std error according to MSE o SMAPE (='error_type')
  #whithin the samples in the test subset using the predictive model
  #with each of the three models used:
  #code = using predictive model
  #nocode = using direct preditor
  #def/default = using default predictor
  
  if (error_type=='MSE'){ #select interesting columns from the df accoding to the type of error selected 
    data_plot=data %>% 
      select(id, tr1_mse, color_variable)
    colnames(data_plot) = c('id', 'error', 'variable')
      
  }else if (error_type=='SMAPE'){
    data_plot=data %>% 
      select(id, tr1_smape, color_variable)
    colnames(data_plot) = c('id', 'error', 'variable')
  }
  
  #plot
  p1 <- ggplot(data_plot, aes(x = variable, y = error)) +
    geom_boxplot(notch=FALSE, fill='grey') +#, aes(group = cut_width(variable, w))) +
    labs(title = paste(error_type, ' values within samples grouped according to ', color_variable, sep='')) +
    labs (x = color_variable, y = error_type) +
    theme(plot.title = element_text(size=20, face="bold"), axis.text.x=element_text(angle=60,hjust=1,size='10')) 
  p1
}


#load data (produced by script table_AppendixI.R)
path='./data/transformed_data/predictive_model/best_model/'
file='/table2save_samplestest.csv'
data = read.table(paste(path, file, sep=''), sep='\t', header = TRUE)

#example
plot_samples_predictivemodel(data, 'MSE', 'Maize_Line', NULL)

