#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 12:46:42 2019

@author: saracabellopinedo
"""

# import modules
import pandas as pd
from sklearn.model_selection import train_test_split


def load_data(path):
    #load predictive model subset
    #input: path where the data is found
    #output: df containing OTU abundance table
    
    original_data = pd.read_csv(path, sep='\t')
    original_data.index=original_data['otuids'].tolist()
    original_data=original_data.drop('otuids', axis=1)
    df_transposed=original_data.T
    return(df_transposed)
    

def shuffle_records(metadata, df_transposed):
    #shuffle data records 
    #returns the same two df but with a different row order
    
    n_col=len(metadata.columns)
    df_joined=metadata.set_index('X.SampleID').join(df_transposed)  
    df_joined=df_joined.sample(frac=1)

    metadata=df_joined.iloc[:, 0:n_col-1]
    df_transposed=df_joined.iloc[:, n_col-1:len(df_joined.columns)+1]
        
    return metadata, df_transposed

def train_test_files(metadata, supervised_data, size):
    #write files with train and test subsets
    #size = size of the test subset (range between 0 and 1)
    
    if 'X.SampleID' in metadata:
        metadata.drop(['X.SampleID'], axis='columns', inplace=True)
        
    X_train, X_test, y_train, y_test = train_test_split(metadata, supervised_data, test_size=size)
    X_train.to_csv("./data/original_data/data_train_predictivemodel/2otu_table2save_80_xtrain80.csv", sep='\t')
    X_test.to_csv("./data/original_data/data_train_predictivemodel/2otu_table2save_80_xtest20.csv", sep='\t')
    y_train.to_csv("./data/original_data/data_train_predictivemodel/2otu_table2save_80_ytrain80.csv", sep='\t')
    y_test.to_csv("./data/original_data/data_train_predictivemodel/2otu_table2save_80_ytest20.csv", sep='\t')

#load OTU table from predictive model subset
original_file = load_data ("./data/original_data/otu_table2save_80.csv")

metadata = pd.read_csv("./data/original_data/metadata_table2save_80_selectedvariables.csv", sep='\t')

metadata, df_otus = shuffle_records(metadata, original_file)

train_test_files(metadata, df_otus, 0.2)
