library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)

plot_error=function(data, error_type, variables, representation_type='n'){
  #function that produces 3 plots comparing the error in each OTU predicted by three models:
  #predictive model, direct predictor and default predictor
  #the plots can be represented using logarithmic scale if representation_type is set to 'log'
  
  #selecting relevant information according to 'error_type', that can be MSE or SMAPE
  if (error_type=='MSE'){
    data_plot=select(data, id, tr1_mse, tr2_mse, def_mse)
    colnames(data_plot)=c('id', 'code', 'no-code', 'default')
  }else if (error_type=='SMAPE'){
    data_plot=select(data, id, tr1_smape, tr2_smape, def_smape)
    colnames(data_plot)=c('id', 'code', 'no-code', 'default')
  }
  
  #logarithmic representation
  if (representation_type=='log'){
    p1 = ggplot(data_plot, aes(x=log(data_plot$code), y=log(data_plot$default))) + 
      geom_point(size=0.75) +
      labs (x = "Using code", y = "Default") + 
      geom_abline(intercept = 0, slope = 1, color="red") + 
      theme(legend.position = "none")
    
    p2 = ggplot(data_plot, aes(x=log(data_plot$`no-code`), y=log(data_plot$default))) + 
      geom_point(size=0.75) +
      labs (x = "Without code", y = "Default" ) + 
      geom_abline(intercept = 0, slope = 1, color="red")  + 
      theme(legend.position = "none")
    
    p3 = ggplot(data_plot, aes(x=log(data_plot$code), y=log(data_plot$`no-code`))) + 
      geom_point(size=0.75) +
      labs (x = "Using code", y = "Without code" ) + 
      geom_abline(intercept = 0, slope = 1, color="red")  + 
      theme(legend.position = "none")
  } else { #normal scale 
    p1 = ggplot(data_plot, aes(x=(data_plot$code), y=(data_plot$default))) + 
      geom_point(size=0.75) +
      labs (x = "Using code", y = "Default") + 
      geom_abline(intercept = 0, slope = 1, color="red") + 
      theme(legend.position = "none")
    
    p2 = ggplot(data_plot, aes(x=(data_plot$`no-code`), y=(data_plot$default))) + 
      geom_point(size=0.75) +
      labs (x = "Without code", y = "Default" ) + 
      geom_abline(intercept = 0, slope = 1, color="red")  + 
      theme(legend.position = "none")
    
    p3 = ggplot(data_plot, aes(x=(data_plot$code), y=(data_plot$`no-code`))) + 
      geom_point(size=0.75) +
      labs (x = "Using code", y = "Without code" ) + 
      geom_abline(intercept = 0, slope = 1, color="red")  + 
      theme(legend.position = "none")
  }
  
  #building the title
  if (representation_type!='n'){
    text=paste(representation_type, error_type, " comparison between predictive strategies.",sep='', collapse = NULL)
  } else{
    text=paste(error_type, " comparison between three predictive strategies", sep='', collapse = NULL)
  }
  title=textGrob(text, gp=gpar(fontface="bold", fontsize=20))
  grid.arrange(p1, p2, p3, nrow=1, top = title,   bottom = textGrob(
    paste("Variables used: ", variables, "  ", sep=''),
    gp = gpar(fontface = 3, fontsize = 10),
    hjust = 1,
    x = 1
  ))
  
}

#load data (produced by script table_AppendixI.R)
path='./data/transformed_data/predictive_model/best_model/'
file='/table2save_otustest.csv'
data = read.table(paste(path, file, sep=''), sep='\t', header = TRUE)

#possibilities
plot_error(data, 'MSE', 'all')
plot_error(data, 'MSE', 'all', 'log')
plot_error(data, 'SMAPE', 'all')
plot_error(data, 'SMAPE', 'all', 'log')