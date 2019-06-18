library(rowr)
library(tidyr)
library(dplyr)

load_transformed_data=function(path, option){
  #load transformed data: train and test subsets, using both direct predictor (nocode) and predictive model (code)
  filename1=paste(path, 'codetr_', option, '_train.csv', sep='', collapse = NULL)
  otus_tr1_train<<-read.table(filename1, sep = '\t', header = TRUE)
  rownames(otus_tr1_train)<<- otus_tr1_train$X
  otus_tr1_train$X<<- NULL
  filename2=paste(path, 'nocodetr_', option, '_train.csv', sep='', collapse = NULL)
  otus_tr2_train<<-read.table(filename2, sep = '\t', header = TRUE)
  rownames(otus_tr2_train) <<- otus_tr2_train$X
  otus_tr2_train$X<<-NULL
  filename3=paste(path, 'codetr_', option, '_test.csv', sep='', collapse = NULL)
  otus_tr1_test<<-read.table(filename3, sep = '\t', header = TRUE)
  rownames(otus_tr1_test)<<- otus_tr1_test$X
  otus_tr1_test$X <<- NULL
  filename4=paste(path, 'nocodetr_', option, '_test.csv', sep='', collapse = NULL)
  otus_tr2_test<<-read.table(filename4, sep = '\t', header = TRUE)
  rownames(otus_tr2_test) <<- otus_tr2_test$X
  otus_tr2_test$X <<- NULL
}

create_groups=function(df){
  #creates groups according to the precission in the prediction of OTUs and samples
  
  diff_mse = as.numeric(as.character(df$diff_mse))
  diff_smape = as.numeric(as.character(df$diff_smape))
  group_mse = numeric(length(diff_mse))
  group_smape = numeric(length(diff_smape))
  for (k in 1:length(diff_mse)){
    if (diff_mse[k] < -0.0001){
      group_mse[k] = 'better'
    } else if (diff_mse[k] > 0.0001){
      group_mse[k] = 'worse'
    } else {
      group_mse[k] = 'similar'
    }
  }
  
  for (k in 1:length(diff_smape)){
    if (diff_smape[k] < -1){
      group_smape[k] = 'better'
    } else if (diff_smape[k] > 1){
      group_smape[k] = 'worse'
    } else {
      group_smape[k] = 'similar'
    }
  }
  
  df_new=cbind(df, group_mse, group_smape)
  return (df_new)
}

calculations_error = function (train_or, test_or, otus_tr1_train, otus_tr2_train, otus_tr1_test, otus_tr2_test){
  #calculates error for samples and OTUs using the three data: original, transformed by predictive model and 
  #transformed using direct predictor
  #MSE and SMAPE measurements are produced
  
  #formatting column names
  columns = colnames (otus_tr1_train)
  new_columns= gsub("X", "", columns)
  colnames(train_or) = new_columns
  colnames(test_or) = new_columns
  colnames(otus_tr1_train) = new_columns
  colnames(otus_tr1_test) = new_columns
  colnames(otus_tr2_train) = new_columns
  colnames(otus_tr2_test) = new_columns
  
  #error calculation for each cell of the four dataframes and the original values
  mse_tr1_train=(train_or-otus_tr1_train)^2
  mse_tr2_train=(train_or-otus_tr2_train)^2
  mse_tr1_test=(test_or-otus_tr1_test)^2
  mse_tr2_test=(test_or-otus_tr2_test)^2
  
  smape_tr1_train=abs(train_or-otus_tr1_train)/(abs(train_or)+abs(otus_tr1_train))
  smape_tr2_train=abs(train_or-otus_tr2_train)/(abs(train_or)+abs(otus_tr2_train))
  smape_tr1_test=abs(test_or-otus_tr1_test)/(abs(test_or)+abs(otus_tr1_test))
  smape_tr2_test=abs(test_or-otus_tr2_test)/(abs(test_or)+abs(otus_tr2_test))
  #na is replace by 0 because na presence is caused by df_or[i,j] + df_tr[i,j] = 0 --> their equal
  smape_tr1_train[is.na(smape_tr1_train)] <- 0
  smape_tr2_train[is.na(smape_tr2_train)] <- 0
  smape_tr1_test[is.na(smape_tr1_test)] <- 0
  smape_tr2_test[is.na(smape_tr2_test)] <- 0
  
  
  #error grouped by SAMPLE
  #train subset
  error_samples_train=NULL
  samples_train=row.names(mse_tr1_train)
  
  for (j in 1:length(samples_train)){
    id_sample=samples_train[j]
    mse_value_tr1=sum(mse_tr1_train[id_sample,])/ncol(mse_tr1_train)
    mse_value_tr2=sum(mse_tr2_train[id_sample,])/ncol(mse_tr2_train)
    smape_value_tr1=(sum(smape_tr1_train[id_sample,])/ncol(smape_tr1_train))*100
    smape_value_tr2=(sum(smape_tr2_train[id_sample,])/ncol(smape_tr2_train))*100
    new_column=c(id_sample, mse_value_tr1, smape_value_tr1, mse_value_tr2, smape_value_tr2)
    error_samples_train=cbind(error_samples_train, new_column)
  }
  
  colnames(error_samples_train)=error_samples_train[1,]
  row.names(error_samples_train) = c('id', 'tr1_mse', 'tr1_smape', 'tr2_mse', 'tr2_smape')
  error_samples_train=data.frame(t(error_samples_train))
  
  #test_subset
  error_samples_test=NULL
  samples_test=row.names(mse_tr1_test)
  
  for (j in 1:length(samples_test)){
    id_sample=samples_test[j]
    mse_value_tr1=sum(mse_tr1_test[id_sample,])/ncol(mse_tr1_test)
    mse_value_tr2=sum(mse_tr2_test[id_sample,])/ncol(mse_tr2_test)
    smape_value_tr1=(sum(smape_tr1_test[id_sample,])/ncol(smape_tr1_test))*100
    smape_value_tr2=(sum(smape_tr2_test[id_sample,])/ncol(smape_tr2_test))*100
    new_column=c(id_sample, mse_value_tr1, smape_value_tr1, mse_value_tr2, smape_value_tr2)
    error_samples_test=cbind(error_samples_test, new_column)
  }
  colnames(error_samples_test)=error_samples_test[1,]
  row.names(error_samples_test) = c('id', 'tr1_mse', 'tr1_smape', 'tr2_mse', 'tr2_smape')
  error_samples_test=data.frame(t(error_samples_test))
  
  #error grouped by OTU
  error_otus_train=NULL
  error_otus_test=NULL
  otu_names=colnames(mse_tr1_test)
  
  for (j in 1:length(otu_names)){
    id_otu=otu_names[j]
    mse_value_train_tr1=sum(mse_tr1_train[,id_otu])/nrow(mse_tr1_train)
    mse_value_train_tr2=sum(mse_tr2_train[,id_otu])/nrow(mse_tr2_train)
    mse_value_test_tr1=sum(mse_tr1_test[,id_otu])/nrow(mse_tr1_test)
    mse_value_test_tr2=sum(mse_tr2_test[,id_otu])/nrow(mse_tr2_test)
    smape_value_train_tr1=sum(smape_tr1_train[,id_otu])/nrow(smape_tr1_train)*100
    smape_value_train_tr2=sum(smape_tr2_train[,id_otu])/nrow(smape_tr2_train)*100
    smape_value_test_tr1=sum(smape_tr1_test[,id_otu])/nrow(smape_tr1_test)*100
    smape_value_test_tr2=sum(smape_tr2_test[,id_otu])/nrow(smape_tr2_test)*100
    new_column_train=c(id_otu, mse_value_train_tr1, smape_value_train_tr1, mse_value_train_tr2, smape_value_train_tr2)
    new_column_test=c(id_otu, mse_value_test_tr1, smape_value_test_tr1, mse_value_test_tr2, smape_value_test_tr2)
    
    error_otus_train=cbind(error_otus_train, new_column_train)
    error_otus_test=cbind(error_otus_test, new_column_test)
  }
  colnames(error_otus_train)=error_otus_train[1,]
  row.names(error_otus_train) = c('id', 'tr1_mse', 'tr1_smape', 'tr2_mse', 'tr2_smape')
  error_otus_train=data.frame(t(error_otus_train))
  
  colnames(error_otus_test)=error_otus_test[1,]
  row.names(error_otus_test) = c('id', 'tr1_mse', 'tr1_smape', 'tr2_mse', 'tr2_smape')
  error_otus_test=data.frame(t(error_otus_test))
  
  #calculate differences to know with transformation is more accurate for the four subsets
  diff_mse=NULL
  diff_smape=NULL
  diff_mse=as.numeric(as.character(error_samples_train$tr1_mse)) - as.numeric(as.character(error_samples_train$tr2_mse))
  diff_smape=as.numeric(as.character(error_samples_train$tr1_smape)) - as.numeric(as.character(error_samples_train$tr2_smape))
  error_samples_train=cbind(error_samples_train, diff_mse, diff_smape)
  error_samples_train=create_groups(error_samples_train)
  
  diff_mse=NULL
  diff_smape=NULL
  diff_mse=as.numeric(as.character(error_samples_test$tr1_mse)) - as.numeric(as.character(error_samples_test$tr2_mse))
  diff_smape=as.numeric(as.character(error_samples_test$tr1_smape)) - as.numeric(as.character(error_samples_test$tr2_smape))
  error_samples_test=cbind(error_samples_test, diff_mse, diff_smape)
  error_samples_test=create_groups(error_samples_test)
  
  diff_mse=NULL
  diff_smape=NULL
  diff_mse=as.numeric(as.character(error_otus_train$tr1_mse)) - as.numeric(as.character(error_otus_train$tr2_mse))
  diff_smape=as.numeric(as.character(error_otus_train$tr1_smape)) - as.numeric(as.character(error_otus_train$tr2_smape))
  error_otus_train=cbind(error_otus_train, diff_mse, diff_smape)
  error_otus_train=create_groups(error_otus_train)
  
  diff_mse=NULL
  diff_smape=NULL
  diff_mse=as.numeric(as.character(error_otus_test$tr1_mse)) - as.numeric(as.character(error_otus_test$tr2_mse))
  diff_smape=as.numeric(as.character(error_otus_test$tr1_smape)) - as.numeric(as.character(error_otus_test$tr2_smape))
  error_otus_test=cbind(error_otus_test, diff_mse, diff_smape)
  error_otus_test=create_groups(error_otus_test)
  
  return(cbind.fill(error_samples_train, error_samples_test, error_otus_train, error_otus_test, fill=NA))
}

load_metadata=function(data, metadata_path, variable){
  #load metadata table and select the interesting variable
  # returns a dataframe that contains sample id + metadata variable value
  
  metadata = read.table(metadata_path, sep = '\t', header = TRUE)
  rownames(metadata) = metadata$X.SampleID
  metadata$X.SampleID = NULL
  
  new_column=NULL
  
  for (i in 1:nrow(data)){
    sample_id=row.names(data)[i]
    new_column=c(new_column, as.character(metadata[sample_id, variable]))
  }
  
  new_cols=c(colnames(data), variable)
  data=cbind(data, new_column)
  colnames(data)=new_cols
  return(data)
}

load_taxa=function(data, taxa_path, taxrank){
  #load taxonomic information and select the interesting taxonomic level (taxrank)
  #similar to the previous function
  
  taxa = read.table(taxa_path, sep = '\t', header = TRUE)
  rownames(taxa) = taxa$otuids
  taxa$otuids = NULL
  
  new_column=NULL
  
  for (i in 1:nrow(data)){
    otu_id=row.names(data)[i]
    new_column=c(new_column, as.character(taxa[otu_id, taxrank]))
  }
  
  new_cols=c(colnames(data), taxrank)
  data=cbind(data, new_column)
  colnames(data)=new_cols
  return(data)
}

complete_information = function (data, metadata_variables, taxa_ranks) {
  #obtain a df that contains all the available information for samples or OTUs
  
  samples_train = data[,1:9] %>% drop_na()
  row.names(samples_train) = samples_train$id
  samples_test = data[,10:18] %>% drop_na()
  row.names(samples_test) = samples_test$id
  otus_train = data[,19:27] %>% drop_na()
  row.names(otus_train) = otus_train$id
  otus_test = data[,28:36] %>% drop_na()
  row.names(otus_test) = otus_test$id
  
  #add metadata info to the sample data
  for (i in 1:length(metadata_variables)){
    samples_train=load_metadata(samples_train, metadata_path, metadata_variables[i])
    samples_test=load_metadata(samples_test, metadata_path, metadata_variables[i])
  }
  
  #add taxonomy to otus data
  for (i in 1:length(taxa_ranks)){
    otus_train=load_taxa(otus_train,  taxa_path, taxa_ranks[i])
    otus_test=load_taxa(otus_test,  taxa_path, taxa_ranks[i])
  }
  
  return(cbind.fill(samples_train, samples_test, otus_train, otus_test, fill=NA))
  
}

better_worse_df = function(diff){
  #function that counts how many OTUs/samples are better or worse predicted according to one of the models
  pos = 0
  neg = 0
  for (i in 1:length(diff)){
    if (diff[i] > 0){
      pos = pos +1
    }
    if (diff[i] < 0){
      neg = neg + 1
    }
  }
  result=c(neg, pos)
  return(result)
}

better_worse = function (df1, df2, df3, df4){
  #obtain the number of better/worse predicted samples/OTUs for each of the 4 subsets 
  #and according to MSE and SMAPE
  result1 = better_worse_df(as.numeric(as.character(df1$diff_mse)))
  result2 = better_worse_df(as.numeric(as.character(df1$diff_smape)))
  result3 = better_worse_df(as.numeric(as.character(df2$diff_mse)))
  result4 = better_worse_df(as.numeric(as.character(df2$diff_smape)))
  result5 = better_worse_df(as.numeric(as.character(df3$diff_mse)))
  result6 = better_worse_df(as.numeric(as.character(df3$diff_smape)))
  result7 = better_worse_df(as.numeric(as.character(df4$diff_mse)))
  result8 = better_worse_df(as.numeric(as.character(df4$diff_smape)))
  return(c(result1, result3, result5, result7, result2, result4, result6, result8))
}

load_default_transformation=function(path){
  #load transformed data by default predictor
  filename1=paste(path, 'ddefault_otutable_test.csv', sep='', collapse = NULL)
  otus_def_train<<-read.table(filename1, sep = '\t', header = TRUE)
  rownames(otus_def_train)<<- otus_def_train$X
  otus_def_train$X<<- NULL
  
  filename3=paste(path, 'default_otutable_train.csv', sep='', collapse = NULL)
  otus_def_test<<-read.table(filename3, sep = '\t', header = TRUE)
  rownames(otus_def_test)<<- otus_def_test$X
  otus_def_test$X <<- NULL
}

calculations_error_def = function (train_or, test_or, otus_def_train, otus_def_test){
  #calculate error measurements of default prediction compared to the original data
  #formatting column names
  columns = colnames (otus_def_train)
  new_columns= gsub("X", "", columns)
  colnames(train_or) = new_columns
  colnames(test_or) = new_columns
  colnames(otus_def_train) = new_columns
  colnames(otus_def_test) = new_columns
  
  #error calculation for each cell of the 2 dataframes and the original values
  mse_def_train=(train_or-otus_def_train)^2
  mse_def_test=(test_or-otus_def_test)^2
  smape_def_train=(abs(train_or-otus_def_train)/(abs(train_or)+abs(otus_def_train)))*100
  smape_def_test=abs(test_or-otus_def_test)/(abs(test_or)+abs(otus_def_test))*100
  smape_def_train[is.na(smape_def_train)] = 0
  smape_def_test[is.na(smape_def_test)] = 0
  
  #error grouped by SAMPLE
  #train subset
  error_samples_train=NULL
  samples_train=row.names(mse_def_train)
  
  for (j in 1:length(samples_train)){
    id_sample=samples_train[j]
    mse_value_def=sum(mse_def_train[id_sample,])/ncol(mse_def_train)
    smape_value_def=sum(smape_def_train[id_sample,])/ncol(smape_def_train)
    new_column=c(id_sample, mse_value_def, smape_value_def)
    error_samples_train=cbind(error_samples_train, new_column)
  }
  
  colnames(error_samples_train)=error_samples_train[1,]
  row.names(error_samples_train) = c('id', 'def_mse', 'def_smape')
  error_samples_train=data.frame(t(error_samples_train))
  
  #test_subset
  error_samples_test=NULL
  samples_test=row.names(mse_def_test)
  
  for (j in 1:length(samples_test)){
    id_sample=samples_test[j]
    mse_value_def=sum(mse_def_test[id_sample,])/ncol(mse_def_test)
    smape_value_def=sum(smape_def_test[id_sample,])/ncol(smape_def_test)
    new_column=c(id_sample, mse_value_def, smape_value_def)
    error_samples_test=cbind(error_samples_test, new_column)
  }
  
  colnames(error_samples_test)=error_samples_test[1,]
  row.names(error_samples_test) = c('id', 'def_mse', 'def_smape')
  error_samples_test=data.frame(t(error_samples_test))
  
  #error grouped by OTU
  error_otus_train=NULL
  error_otus_test=NULL
  otu_names=colnames(mse_def_test)
  
  for (j in 1:length(otu_names)){
    id_otu=otu_names[j]
    mse_value_train=sum(mse_def_train[,id_otu])/nrow(mse_def_train)
    mse_value_test=sum(mse_def_test[,id_otu])/nrow(mse_def_test)
    smape_value_train=sum(smape_def_train[,id_otu])/nrow(smape_def_train)
    smape_value_test=sum(smape_def_test[,id_otu])/nrow(smape_def_test)
    new_column_train=c(id_otu, mse_value_train, smape_value_train)
    new_column_test=c(id_otu, mse_value_test, smape_value_test)
    
    error_otus_train=cbind(error_otus_train, new_column_train)
    error_otus_test=cbind(error_otus_test, new_column_test)
  }
  colnames(error_otus_train)=error_otus_train[1,]
  row.names(error_otus_train) = c('id', 'def_mse', 'def_smape')
  error_otus_train=data.frame(t(error_otus_train))
  
  colnames(error_otus_test)=error_otus_test[1,]
  row.names(error_otus_test) = c('id', 'def_mse', 'def_smape')
  error_otus_test=data.frame(t(error_otus_test))
  
  return(cbind.fill(error_samples_train, error_samples_test, error_otus_train, error_otus_test, fill=NA))
}

default_analysis = function (original_otutable_train, original_otutable_test){
  #function that builds the predicted OTU tables (train and test) by the default predictor,
  #and uses them to call the funcion that calculates error measurements
  
  #default OTU table train subset calculations
  default_otutable_train=NULL
  for (i in (1:ncol(original_otutable_train))){
    otu_abundance=as.numeric(original_otutable_train[,i])
    mean_value=mean(otu_abundance)
    new_column=rep.int(mean_value, nrow(original_otutable_train))
    default_otutable_train=cbind(default_otutable_train, new_column)
  }
  colnames(default_otutable_train) = colnames(original_otutable_train)
  row.names(default_otutable_train) = row.names(original_otutable_train)
  
  #the same for the test subset
  default_otutable_test = NULL
  for (i in (1:ncol(original_otutable_test))){
    otu_id=colnames(original_otutable_test)[i]
    mean_value= mean(as.numeric(default_otutable_train[,i]))
    new_column=rep.int(mean_value, nrow(original_otutable_test))
    default_otutable_test=cbind(default_otutable_test, new_column)
  }
  colnames(default_otutable_test) = colnames(original_otutable_test)
  row.names(default_otutable_test) = row.names(original_otutable_test)
  
  #calculate error measurements
  def_error = calculations_error_def(train_or, test_or, default_otutable_train, default_otutable_test)
  return(def_error)
}

compared2default_values = function(df){
  #function that compares the result of the prediction by predictive model or direct predictor to the results
  #obtained using the default predictor
  #returns a df wiht the count of the better (than default) predicted OTUs/samples using each method and according to MSE and SMAPE criterion 
  
  result=NULL
  rows=as.character(df$id)
  result=matrix(,length(rows), 4)
  colnames(result)=c('better_mse_codetr', 'better_mse_directtr', 'better_smape_codetr', 'better_smape_directtr')
  rownames(result)=rows
  
  for (i in(1:length(rows))) {
    result[i, 'better_mse_codetr'] = as.numeric(as.character(df$def_mse))[i] - as.numeric(as.character(df$tr1_mse))[i]
    result[i, 'better_mse_directtr'] = as.numeric(as.character(df$def_mse))[i] - as.numeric(as.character(df$tr2_mse))[i]
    result[i, 'better_smape_codetr'] = as.numeric(as.character(df$def_smape))[i] - as.numeric(as.character(df$tr1_smape))[i]
    result[i, 'better_smape_directtr'] = as.numeric(as.character(df$def_smape))[i] - as.numeric(as.character(df$tr2_smape))[i] 
  }
  return(result)
} 

write_individual_files = function (train_or, test_or, metadata_variables, metadata_path, taxa_ranks, taxa_path, type, def_error){
  #function that writes individual files for each model that contain all the analysis carried out
  
  path=paste('./from_metadata/', type, sep='')
  dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
  
  #load transformed data
  load_transformed_data('./transformed_data/mapping_variables_study/', type)
  
  #error calculations
  error_data = calculations_error (train_or, test_or, otus_tr1_train, otus_tr2_train, otus_tr1_test, otus_tr2_test)
  
  #writing files with completed information for plotting
  error_data_completed = complete_information(error_data, metadata_variables, taxa_ranks) 
  
  samples_train = error_data_completed[,1:15] %>% drop_na()
  rownames(samples_train) = samples_train$id
  samples_train_def = def_error[,1:3] %>% drop_na()
  rownames(samples_train_def) = samples_train_def$id
  samples_train_c <<- merge(samples_train, samples_train_def, by.x='id', by.y='id')
  write.table(samples_train_c, paste(path, '/table2save_samplestrain.csv', sep=''), sep='\t',  row.names = FALSE, col.names = TRUE)
  
  samples_test = error_data_completed[,16:30] %>% drop_na()
  rownames(samples_test) = samples_test$id
  samples_test_def = def_error[,4:6] %>% drop_na()
  rownames(samples_test_def) = samples_test_def$id
  samples_test_c <<- merge(samples_test, samples_test_def, by.x='id', by.y='id')
  write.table(samples_test_c, paste(path, '/table2save_samplestest.csv', sep=''), sep='\t',  row.names = FALSE, col.names = TRUE)
  
  otus_train = error_data_completed[,31:46] %>% drop_na()
  rownames(otus_train) = otus_train$id
  otus_train_def = def_error[,7:9] %>% drop_na()
  rownames(otus_train_def) = otus_train_def$id
  otus_train_c <<- merge(otus_train, otus_train_def, by.x='id', by.y='id')
  write.table(otus_train_c, paste(path, '/table2save_otustrain.csv', sep=''), sep='\t',  row.names = FALSE, col.names = TRUE)
  
  otus_test <<- error_data_completed[,47:62] %>% drop_na()
  rownames(otus_test) = otus_test$id
  otus_test_def = def_error[,10:12] %>% drop_na()
  rownames(otus_test_def) = otus_test_def$id
  otus_test_c <<- merge(otus_test, otus_test_def, by.x='id', by.y='id')
  write.table(otus_test_c, paste(path, '/table2save_otustest.csv', sep=''), sep='\t',  row.names = FALSE, col.names = TRUE)
  
  #compared to default predictor
  samplestrain_vs_default=compared2default_values(samples_train_c)
  samplestrain_vs_default=cbind(row.names(samplestrain_vs_default), samplestrain_vs_default)
  write.table(samplestrain_vs_default,  paste(path, '/samplestrain_vs_default_values.csv', sep=''), row.names=FALSE, col.names=TRUE, sep='\t')
  
  samplestest_vs_default=compared2default_values(samples_test_c)
  samplestest_vs_default=cbind(row.names(samplestest_vs_default), samplestest_vs_default)
  write.table(samplestest_vs_default,  paste(path, '/samplestest_vs_default_values.csv', sep=''), row.names=FALSE, col.names=TRUE, sep='\t')
  
  otustrain_vs_default=compared2default_values(otus_train_c)
  otustrain_vs_default=cbind(row.names(otustrain_vs_default), otustrain_vs_default)
  write.table(otustrain_vs_default,  paste(path, '/otustrain_vs_default_values.csv', sep=''), row.names=FALSE, col.names=TRUE, sep='\t')
  
  otustest_vs_default=compared2default_values(otus_test_c)
  otustest_vs_default=cbind(row.names(otustest_vs_default), otustest_vs_default)
  write.table(otustest_vs_default,  paste(path, '/otustest_vs_default_values.csv', sep=''), row.names=FALSE, col.names=TRUE, sep='\t')
  
}

setwd('/Volumes/USB\ DISK')

#load original data
train_or = read.table('./data/original_data/data_train_predictivemodel/2otu_table2save_80_ytrain80_norm.csv', sep='\t', header = TRUE)
rownames(train_or) = train_or$X.SampleID
train_or$X.SampleID = NULL
train_or$X.SampleID.1 = NULL

test_or = read.table('./data/original_data/data_train_predictivemodel/2otu_table2save_80_ytest20_norm.csv', sep='\t', header = TRUE)
rownames(test_or) = test_or$X.SampleID
test_or$X.SampleID = NULL
test_or$X.SampleID.1 = NULL

#contains all the possible types of models build according to the variables used as input
#-all --> all variables used
#-AGE --> all variables escept AGE
possible_types=c('-all', '-elevation', '-AGE', '-Temperature', '-Precipitation3Days', '-INBREDS', '-Maize_Line', '-elevation&AGE', '-elevation&Temperature', '-elevation&Precipitation3Days', 
                 '-elevation&INBREDS', '-elevation&Maize_Line', '-AGE&Temperature', '-AGE&Precipitation3Days', '-AGE&INBREDS', '-AGE&Maize_Line', '-Temperature&Precipitation3Days',
                 '-Temperature&INBREDS', '-Temperature&Maize_Line', '-Precipitation3Days&INBREDS', '-Precipitation3Days&Maize_Line', '-INBREDS&Maize_Line', '-elevation&AGE&Temperature',
                 '-elevation&AGE&Precipitation3Days', '-elevation&AGE&INBREDS', '-elevation&AGE&Maize_Line', '-elevation&Temperature&Precipitation3Days', '-elevation&Temperature&INBREDS',
                 '-elevation&Temperature&Maize_Line', '-elevation&Precipitation3Days&INBREDS', '-elevation&Precipitation3Days&Maize_Line', '-elevation&INBREDS&Maize_Line',
                 '-AGE&Temperature&Precipitation3Days', '-AGE&Temperature&INBREDS', '-AGE&Temperature&Maize_Line', '-AGE&Precipitation3Days&INBREDS', '-AGE&Precipitation3Days&Maize_Line',
                 '-AGE&INBREDS&Maize_Line', '-Temperature&Precipitation3Days&INBREDS', '-Temperature&Precipitation3Days&Maize_Line', '-Temperature&INBREDS&Maize_Line',
                 '-Precipitation3Days&INBREDS&Maize_Line', '-elevation&AGE&Temperature&Precipitation3Days', '-elevation&AGE&Temperature&INBREDS', '-elevation&AGE&Temperature&Maize_Line',
                 '-elevation&AGE&Precipitation3Days&INBREDS', '-elevation&AGE&Precipitation3Days&Maize_Line', '-elevation&AGE&INBREDS&Maize_Line', '-elevation&Temperature&Precipitation3Days&INBREDS',
                 '-elevation&Temperature&Precipitation3Days&Maize_Line', '-elevation&Temperature&INBREDS&Maize_Line', '-elevation&Precipitation3Days&INBREDS&Maize_Line',
                 '-AGE&Temperature&Precipitation3Days&INBREDS', '-AGE&Temperature&Precipitation3Days&Maize_Line', '-AGE&Temperature&INBREDS&Maize_Line', '-AGE&Precipitation3Days&INBREDS&Maize_Line',
                 '-Temperature&Precipitation3Days&INBREDS&Maize_Line', '-elevation&AGE&Temperature&Precipitation3Days&INBREDS', '-elevation&AGE&Temperature&Precipitation3Days&Maize_Line',
                 '-elevation&AGE&Temperature&INBREDS&Maize_Line', '-elevation&AGE&Precipitation3Days&INBREDS&Maize_Line', '-elevation&Temperature&Precipitation3Days&INBREDS&Maize_Line',
                 '-AGE&Temperature&Precipitation3Days&INBREDS&Maize_Line')

#variables included as input for the models, related to the previuos array
variables=c("Elevation, Age, Temperature, Precipitation3Days, Inbreds, Maize_Line",
            "Age, Temperature, Precipitation3Days, Inbreds, Maize_Line",
            "Elevation, Temperature, Precipitation3Days, Inbreds, Maize_Line",
            "Elevation, Age, Precipitation3Days, Inbreds, Maize_Line",
            "Elevation, Age, Temperature, Inbreds, Maize_Line",
            "Elevation, Age, Temperature, Precipitation3Days, Maize_Line",
            "Elevation, Age, Temperature, Precipitation3Days, Inbreds",
            "Temperature, Precipitation3Days, Inbreds, Maize_Line",
            "Age,  Precipitation3Days, Inbreds, Maize_Line",
            " Age, Temperature, Inbreds, Maize_Line",
            " Age, Temperature, Precipitation3Days,  Maize_Line",
            "Age, Temperature, Precipitation3Days, Inbreds",
            "Elevation, Precipitation3Days, Inbreds, Maize_Line",
            "Elevation, Temperature, Inbreds, Maize_Line",
            "Elevation, Temperature, Precipitation3Days, Maize_Line",
            "Elevation, Temperature, Precipitation3Days, Inbreds",
            "Elevation, Age, Inbreds, Maize_Line",
            "Elevation, Age, Precipitation3Days,  Maize_Line",
            "Elevation, Age, Precipitation3Days, Inbreds",
            "Elevation, Age, Temperature, Maize_Line",
            "Elevation, Age, Temperature, Inbreds",
            "Elevation, Age, Temperature, Precipitation3Days",
            "Precipitation3Days, Inbreds, Maize_Line",
            "Temperature, Inbreds, Maize_Line",
            "Temperature, Precipitation3Days, Maize_Line",
            "Temperature, Precipitation3Days, Inbreds",
            "Age, Inbreds, Maize_Line",
            "Age, Precipitation3Days, Maize_Line",
            "Age, Precipitation3Days, Inbreds",
            "Age, Temperature, Maize_Line",
            "Age, Temperature, Inbreds",
            "Age, Temperature, Precipitation3Days",
            "Elevation, Inbreds, Maize_Line",
            "Elevation, Precipitation3Days, Maize_Line",
            "Elevation, Precipitation3Days, Inbreds",
            "Elevation, Temperature, Maize_Line",
            "Elevation, Temperature, Inbreds",
            "Elevation, Temperature, Precipitation3Days",
            "Elevation, Age, Maize_Line",
            "Elevation, Age, Inbreds",
            "Elevation, Age, Precipitation3Days",
            "Elevation, Age, Temperature",
            "Inbreds, Maize_Line",
            "Precipitation3Days, Maize_Line",
            "Precipitation3Days, Inbreds",
            "Temperature, Maize_Line",
            "Temperature, Inbreds",
            "Temperature, Precipitation3Days",
            "Age, Maize_Line",
            "Age, Inbreds",
            "Age, Precipitation3Days",
            "Age, Temperature",
            "Elevation, Maize_Line",
            "Elevation, Inbreds",
            "Elevation, Precipitation3Days",
            "Elevation, Temperature",
            "Elevation, Age",
            "Maize_Line",
            "Inbreds",
            "Precipitation3Days",
            "Temperature",
            "Age",
            "Elevation")

#matrix that will save how many samples/OTUs are better/worse predicted by predictive model and direct predictor
betterworse=matrix(,length(variables),16)
colnames(betterworse)=c('mse_samples_train_betterwithcode', 'mse_samples_train_worsewithcode', 'mse_samples_test_betterwithcode', 'mse_samples_test_worsewithcode', 
                        'mse_otus_train_betterwithcode', 'mse_otus_train_worsewithcode', 'mse_otus_test_betterwithcode', 'mse_otus_test_worsewithcode',
                        'smape_samples_train_betterwithcode', 'smape_samples_train_worsewithcode', 'smape_samples_test_betterwithcode', 'smape_samples_test_worsewithcode', 
                        'smape_otus_train_betterwithcode', 'smape_otus_train_worsewithcode', 'smape_otus_test_betterwithcode', 'smape_otus_test_worsewithcode') 

#table of Annex I
final_table=matrix(,length(variables),4)
colnames(final_table)=c('mse_nocode', 'smape_nocode', 'mse_code', 'smape_code')

metadata_path='./data/original_data/metadata_table2save_80.csv'
taxa_path='./data/original_data/tax_table2save_80.csv'
metadata_variables=c('elevation', 'AGE', 'Temperature', 'Precipitation3Days', 'INBREDS', 'Maize_Line')
taxa_ranks=c("Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7")

#running default predictor
def_error = default_analysis(train_or, test_or)

i=1
for (type in possible_types){   #for each possible combination of input variables
  #calculations for each model
  write_individual_files(train_or, test_or, metadata_variables, metadata_path, taxa_ranks, taxa_path, type, def_error)
  
  #counting better/worse OTUs/samples
  betterworse[i,]=better_worse(samples_train_c, samples_test_c, otus_train_c, otus_test_c)
  
  #MSE and SMAPE calculation for transformations using predictive model and direct predictor
  code_mse=as.numeric(as.character(samples_test_c$tr1_mse))  
  code_smape=as.numeric(as.character(samples_test_c$tr1_smape))  
  nocode_mse=as.numeric(as.character(samples_test_c$tr2_mse))  
  nocode_smape=as.numeric(as.character(samples_test_c$tr2_smape))  
  #mean and std for each model and error measurements are saved
  final_table[i,]=c(paste(mean(nocode_mse), '±', sd(nocode_mse), sep=' '), paste(mean(nocode_smape), '±', sd(nocode_smape), sep=' '), paste(mean(code_mse), '±', sd(code_mse), sep=' '), paste(mean(code_smape), '±', sd(code_smape), sep=' '))
  i=i+1
}

#MSE and SMAPE calculations for default predictor
mse_def = as.numeric(as.character(samples_test_c$def_mse))
smape_def = as.numeric(as.character(samples_test_c$def_smape))
string_mse_def = paste(mean(mse_def), '±', sd(mse_def), sep=' ')
string_smape_def = paste(mean(smape_def), '±', sd(smape_def), sep=' ')
mse_def = rep(string_mse_def, i-1)
smape_def = rep(string_smape_def, i-1)

#addition of default predictor results to the table
final_table=cbind(final_table, mse_def, smape_def)

#write global files
final_table=cbind(variables, final_table)
write.table(final_table, './Results/maize80perc/model2/from_metadata/table3.3_700.csv', sep='\t',  row.names = FALSE, col.names = TRUE)

betterworse=cbind(variables, betterworse)
write.table(betterworse, './Results/maize80perc/model2/from_metadata/betterworse_2hn.csv', sep='\t',  row.names = FALSE, col.names = TRUE)

