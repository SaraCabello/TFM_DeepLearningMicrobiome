default_predictor = function (original_otutable_train, original_otutable_test){
  #model that assigns to each OTU predicted value the mean value within all samples in train subset
  
  #transforming train subset using default predictor
  default_otutable_train=NULL
  for (i in (1:ncol(original_otutable_train))){
    otu_abundance=as.numeric(original_otutable_train[,i])
    mean_value=mean(otu_abundance)
    new_column=rep.int(mean_value, nrow(original_otutable_train))
    default_otutable_train=cbind(default_otutable_train, new_column)
  }
  colnames(default_otutable_train) = colnames(original_otutable_train)
  row.names(default_otutable_train) = row.names(original_otutable_train)
  
  #transforming test subset using default predictor
  default_otutable_test = NULL
  for (i in (1:ncol(original_otutable_test))){
    otu_id=colnames(original_otutable_test)[i]
    mean_value= mean(as.numeric(default_otutable_train[,i]))
    new_column=rep.int(mean_value, nrow(original_otutable_test))
    default_otutable_test=cbind(default_otutable_test, new_column)
  }
  colnames(default_otutable_test) = colnames(original_otutable_test)
  row.names(default_otutable_test) = row.names(original_otutable_test)
  
  #write transformed data to files
  default_otutable_train = cbind(row.names(default_otutable_train), default_otutable_train)
  write.table(default_otutable_train, './data/transformed_data/default_predictor/default_otutable_train.csv', sep='\t', col.names=TRUE, row.names = FALSE)
  
  default_otutable_test = cbind(row.names(default_otutable_test), default_otutable_test)
  write.table(default_otutable_test, './data/transformed_data/default_predictor/default_otutable_test.csv', sep='\t', col.names=TRUE, row.names = FALSE)
  
}

setwd('/Volumes/USB\ DISK/TFM_DeepLearningMicrobiome')

#load original data
train_or = read.table('./data/original_data/data_train_predictivemodel/2otu_table2save_80_ytrain80_norm.csv', sep='\t', header = TRUE)
rownames(train_or) = train_or$X.SampleID
train_or$X.SampleID = NULL
train_or$X.SampleID.1 = NULL

test_or = read.table('./data/original_data/data_train_predictivemodel/2otu_table2save_80_ytest20_norm.csv', sep='\t', header = TRUE)
rownames(test_or) = test_or$X.SampleID
test_or$X.SampleID = NULL
test_or$X.SampleID.1 = NULL

default_predictor(train_or, test_or)
