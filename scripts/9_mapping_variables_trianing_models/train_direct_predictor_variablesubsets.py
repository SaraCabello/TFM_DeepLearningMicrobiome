#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 19:16:17 2019

@author: saracabellopinedo
"""

#import modules
from numpy.random import seed
seed(1)
from tensorflow import set_random_seed
set_random_seed(2)

import re
import pandas as pd
import numpy as np
from keras.models import Model
from keras.layers import Input, Dense
from sklearn.metrics import mean_squared_error


def load_subsets(original_path):
     #function to load train and test subsets for x and y:
    #x = environmental variables
    #y = OTU table (that in this models acts as the class to be predicted)
    
    regex = re.compile(".biom$")
    path = re.sub(regex, "", original_path)
    
    xtrain_file=path+"_xtrain80.csv"
    X_train = pd.read_csv(xtrain_file, sep='\t')
    ytrain_file=path+"_ytrain80.csv"
    y_train = pd.read_csv(ytrain_file, sep='\t')
    xtest_file=path+"_xtest20.csv"
    X_test = pd.read_csv(xtest_file, sep='\t')
    ytest_file=path+"_ytest20.csv"
    y_test = pd.read_csv(ytest_file, sep='\t')

    return X_train, y_train, X_test, y_test

def normalize_interval(data):
    #function for normalizing OTU table within [0, 1] interval
    #input: df with original values
    #output: df with normalized values
    
    #create empty df to allocate normalized data
    normalized_data=pd.DataFrame(columns=data.columns) #[1:len(data.columns)+1])
    normalized_data['X.SampleID']=data['X.SampleID']
    normalized_data.index=normalized_data['X.SampleID']
    data.index=data['X.SampleID']
    
    for row in data.itertuples(index=True, name='Pandas'):
        id_sample=row[1]
        values=[]
        values=np.array(data.loc[id_sample][1:len(row)+2])
        min_val=np.min(values)
        max_val=np.max(values)

        if (min_val != max_val):
            try:
                new_values=[(x - min_val)/ (max_val- min_val) for x in values]
            except:
                print("Se han encontrado strings")
                new_values=values
        else:
            new_values=values
        

        new_values.insert(0, id_sample)
        normalized_data.at[id_sample]=new_values
    
    return normalized_data

def normalize_maximum(data):
    #normalize over the maximum value for each mapping variable
    normalized_data=pd.DataFrame(columns=data.columns)
    normalized_data['X.SampleID']=data['X.SampleID']
    for colname in data:
        if colname != 'X.SampleID':
            values=data[colname].tolist()
            maximum=max(values)
            if maximum > 1:
                try:
                    new_values=[x / maximum for x in values]
                except:
                    print("Se han encontrado strings")
                    new_values=values
            else:
                new_values=values
            normalized_data[colname]=new_values
    return normalized_data

def categorical2numbers(X_train, X_test):
    #transform each value of a categorical variable into a number
    n_train=len(X_train)
    metadata=X_train.append(X_test)
    #substitution by numbers
    if 'INBREDS' in metadata.columns:
        #dictionary that stores all inbreds possible values and its corresponding number
        dict_inbreds={'Oh7B':0, 'P39':1, 'CML333':2, 'Il14H':3, 'MS71':4, 'Oh43': 5, 'CML52':6, 
                      'Mo18W':7, 'M37W':8, 'Mo17':9, 'CML103':10, 'CML69':11, 'CML228':12,
                      'CML322':13, 'Tzi8':14, 'B73':15, 'NC350':16, 'CML247':17, 'B97':18,
                      'Ki3':19, 'Ky21':20, 'CML277':21, 'M162W':22, 'Hp301':23, 'NC358':24, 
                      'Tx303':25, 'Ki11':26}
        
        for key in dict_inbreds.keys():
            try:
                metadata=metadata.replace(key, dict_inbreds[key])
            except:
                pass
    
    if 'Maize_Line' in metadata.columns:
        #dictionary that stores maize line possible values and its corresponding number
        dict_maizeline={'Non_Stiff_Stalk':0, 'Sweet_Corn':1, 'Tropical':2, 'Mixed':3, 'Stiff_Stalk':4, 'Popcorn':5}   
        for key in dict_maizeline.keys():
            try:
                metadata=metadata.replace(key, dict_maizeline[key])
            except:
                pass
    X_train = metadata.iloc[0:n_train,:]
    X_test = metadata.iloc[n_train:len(metadata)+1,:]
    return X_train, X_test

def categorical2columns(X_train_df, X_test_df, n):
     #the n values found in a particular categorical variable are transformed into n columns 
    #which are filled using 0 (= not belonging to the variable level) or 1 (= belonging)
    
    n_train=len(X_train_df)
    metadata=X_train_df.append(X_test_df)

    #create columns for each value
    df_new=metadata.iloc[:, 0:n+1]
    metadata.index=metadata['X.SampleID']
    df_new.index=df_new['X.SampleID']
    
    if 'INBREDS' in metadata.columns:
        new_columns1 = ['INBREDS_Oh7B', 'INBREDS_P39', 'INBREDS_CML333', 'INBREDS_Il14H', 'INBREDS_MS71', 'INBREDS_Oh43', 'INBREDS_CML52', 
                      'INBREDS_Mo18W', 'INBREDS_M37W', 'INBREDS_Mo17', 'INBREDS_CML103', 'INBREDS_CML69', 'INBREDS_CML228',
                      'INBREDS_CML322', 'INBREDS_Tzi8', 'INBREDS_B73', 'INBREDS_NC350', 'INBREDS_CML247', 'INBREDS_B97',
                      'INBREDS_Ki3', 'INBREDS_Ky21', 'INBREDS_CML277', 'INBREDS_M162W', 'INBREDS_Hp301', 'INBREDS_NC358', 
                      'INBREDS_Tx303', 'INBREDS_Ki11']
        for column in new_columns1:
            df_new[column]=0
        
        for row in metadata.itertuples(index=True, name='Pandas'):
            id_sample=row[1]
            inbred=metadata.loc[id_sample, 'INBREDS']
            inbred_regexp=re.compile(str(".*" + str(inbred)))
            column_inbred=list(filter(inbred_regexp.match, new_columns1))[0]
            df_new.at[id_sample, column_inbred] = 1
            
    
    if 'Maize_Line' in metadata.columns:
        new_columns2 = ['LINE_Non_Stiff_Stalk', 'LINE_Sweet_Corn', 'LINE_Tropical', 'LINE_Mixed', 
                      'LINE_Stiff_Stalk', 'LINE_Popcorn']
        for column in new_columns2:
            df_new[column]=0
         
        for row in metadata.itertuples(index=True, name='Pandas'):
            id_sample=row[1]
            maize_line=metadata.loc[id_sample, 'Maize_Line']
            maizeline_regexp=re.compile(str(".*" + str(maize_line)))
            column_maizeline=list(filter(maizeline_regexp.match, new_columns2))[0]
            df_new.at[id_sample, column_maizeline] = 1
            
    metadata=df_new   
    X_train = metadata.iloc[0:n_train,:]
    X_test = metadata.iloc[n_train:len(metadata)+1,:]
    return X_train, X_test

def categorical2binary(X_train, X_test, n):
    #each value of a particular categorical variable is codified using binary codification
    #and n new columns are created were n = number of necessary binary digits
    
    n_train=len(X_train)
    metadata=X_train.append(X_test)

    #create columns for each value
    df_new=metadata.iloc[:, 0:n+1]
    metadata.index=metadata['X.SampleID']
    df_new.index=df_new['X.SampleID']
    
    if 'INBREDS' in metadata.columns:
        #column names created in the presence of variable Inbreds
        dict_inbreds={'Oh7B':'00001', 'P39':'00010', 'CML333':'00011', 'Il14H':'00100', 'MS71':'00101', 'Oh43': '00110', 'CML52':'00111', 
                      'Mo18W':'01000', 'M37W':'01001', 'Mo17':'01010', 'CML103':'01011', 'CML69':'01100', 'CML228':'01101',
                      'CML322':'01110', 'Tzi8':'01111', 'B73':'10000', 'NC350':'10001', 'CML247':'10010', 'B97':'10011',
                      'Ki3':'10100', 'Ky21':'10101', 'CML277':'10110', 'M162W':'10111', 'Hp301':'11000', 'NC358':'11001', 
                      'Tx303':'11010', 'Ki11':'11011'}
        new_columns1 = ['INBREDS_4', 'INBREDS_3', 'INBREDS_2', 'INBREDS_1', 'INBREDS_0']
        
        for column in new_columns1:
            df_new[column]=0

        for row in metadata.itertuples(index=True, name='Pandas'):
            id_sample=row[1]
            inbred=metadata.loc[id_sample, 'INBREDS']
            code_inbred=list(dict_inbreds[inbred])
            df_new.at[id_sample, 'INBREDS_4'] = code_inbred[0]
            df_new.at[id_sample, 'INBREDS_3'] = code_inbred[1]
            df_new.at[id_sample, 'INBREDS_2'] = code_inbred[2]
            df_new.at[id_sample, 'INBREDS_1'] = code_inbred[3]
            df_new.at[id_sample, 'INBREDS_0'] = code_inbred[4]   
            
    if 'Maize_Line' in metadata.columns:
        #column names created in the presence of variable Maize_Line
        dict_maizeline={'Non_Stiff_Stalk':'001', 'Sweet_Corn':'010', 'Tropical':'011', 'Mixed':'100', 'Stiff_Stalk':'101', 'Popcorn':'110'}   
        new_columns2 =[ 'LINE_2', 'LINE_1', 'LINE_0']
        
        for column in new_columns2:
            df_new[column]=0
    
        for row in metadata.itertuples(index=True, name='Pandas'):
            id_sample=row[1]
            maize_line=metadata.loc[id_sample, 'Maize_Line']
            code_maizeline=list(dict_maizeline[maize_line])
            df_new.at[id_sample, 'LINE_2'] = code_maizeline[0]    
            df_new.at[id_sample, 'LINE_1'] = code_maizeline[1]    
            df_new.at[id_sample, 'LINE_0'] = code_maizeline[2]    
    
    metadata=df_new   
    
    X_train = metadata.iloc[0:n_train,:]
    X_test = metadata.iloc[n_train:len(metadata)+1,:]
    return X_train, X_test

def delete_variable(variable_names, X_train, X_test):
    #function that delete from dataframe with all original mapping variables (X_train and X_test) 
    #the variables not included in 'variable_names'
    #'variable_names' = array that contains the name of the columns to preserve in the data
    X_train=X_train.loc[:, variable_names]
    X_test=X_test.loc[:, variable_names]
    n=len(X_train.columns)
    if 'INBREDS' in X_train.columns and 'Maize_Line' in X_train.columns:
        n=n-2
    elif 'INBREDS' in X_train.columns and 'Maize_Line' not in X_train.columns:
        n=n-1
    elif 'INBREDS' not in X_train.columns and 'Maize_Line' in X_train.columns:
        n=n-1
    else:
        n=n
        
    return X_train, X_test, n

def format_data(X_train, y_train, X_test, y_test):
    #function to transform the 4 df required to train the model to fit the requirements of
    #training cycles of deep learning in keras
    
    #removing columns that contain the id of the sample, not interesting for the model
    if 'X.SampleID' in X_train:
        X_train.drop(['X.SampleID'], axis='columns', inplace=True)

    if 'X.SampleID' in X_test:
        X_test.drop(['X.SampleID'], axis='columns', inplace=True) 
        
    if 'X.SampleID' in y_train:
        y_train.drop(['X.SampleID'], axis='columns', inplace=True)
        
    if 'X.SampleID' in y_test:
        y_test.drop(['X.SampleID'], axis='columns', inplace=True) 
        
    #transforming df into list of lists (format requiered by keras models)
    if not isinstance(X_train, np.ndarray):
        X_train_list=[]
        for row in X_train.itertuples(index=True, name='Pandas'):
                X_train_list.append(row[1:])
        X_train = np.asarray(X_train)
    
    if not isinstance(X_test, np.ndarray):
        X_test_list=[]
        for row in X_test.itertuples(index=True, name='Pandas'):
                X_test_list.append(row[1:])  
        X_test = np.asarray(X_test_list)
        
    if not isinstance(y_train, np.ndarray):
        y_train_list=[]
        for row in y_train.itertuples(index=True, name='Pandas'):
            y_train_list.append(row[1:])
        y_train=np.asarray(y_train_list)
    
    if not isinstance(y_test, np.ndarray):
        y_test_list=[]
        for row in y_test.itertuples(index=True, name='Pandas'):
            y_test_list.append(row[1:])  
        y_test = np.asarray(y_test_list)
    
    return X_train, y_train, X_test, y_test

def train_model_nolayers(X_train_df, y_train_df, X_test_df, y_test_df, n_iter, activation_function):
    #train model that predicts OTU composition from environmental data witout code and with no hidden 
    #layers between them
    #returns the trained model and the history of error measurements during the training steps
    
    X_train, y_train, X_test, y_test = format_data(X_train_df, y_train_df, X_test_df, y_test_df)

    #   #model definition
    original_dim=len(X_train[0])
    encoding_dim = len (y_train[0])
    
    input_layer= Input(shape=(original_dim,))
    
    encoded = Dense (encoding_dim, activation=activation_function)(input_layer)
    model = Model (input_layer, encoded)
    model.compile(optimizer='adam', loss='mse')
    
    result=model.fit(X_train, y_train,
                    nb_epoch=n_iter,
                    batch_size=256,
                    shuffle=True,
                    validation_data=(X_test, y_test),
                    verbose=0)

    return model, result   


def calculate_OTUs_mse(X_test_df, y_test_df, model):
    #calculate mse between original OTU composition (y_test_df) and the predicted one by the model
    #X_test_df: environmental variables
    normalized_xtest_list=[]
    for row in X_test_df.itertuples(index=True, name='Pandas'):
        normalized_xtest_list.append(row[1:])
    
    normalized_xtest_list=np.asarray(normalized_xtest_list) 
        
    predicted_OTUs=model.predict(normalized_xtest_list)
        
    df_predictedOTUs = pd.DataFrame(data=predicted_OTUs,
                                  index=y_test_df.index,
                                  columns=y_test_df.columns.values.tolist())
    
    mse_otus = mean_squared_error(y_test_df, df_predictedOTUs)
    return mse_otus #returns mse global value

def calculate_smape(X_test_df, y_test_df, model):
    #calculate smape between original OTU composition (y_test_df) and the predicted one 
    #X_test_df: environmental variables
    normalized_xtest_list=[]
    for row in X_test_df.itertuples(index=True, name='Pandas'):
        normalized_xtest_list.append(row[1:])
    
    normalized_xtest_list=np.asarray(normalized_xtest_list) 
        
    predicted_OTUs=model.predict(normalized_xtest_list)
        
    df_predictedOTUs = pd.DataFrame(data=predicted_OTUs,
                                  index=y_test_df.index,
                                  columns=y_test_df.columns.values.tolist())
    
   #A=real value, F=forecast value
    A=np.array(y_test_df.values.flatten())
    F=np.array(df_predictedOTUs.values.flatten())
    return 100./len(A) * np.sum(np.abs(F - A) / (np.abs(A) + np.abs(F))) #returns smape global value

def transform_data(X_train_df, X_test_df, y_train_df, y_test_df, option, n):
    #function to transform original data. 2 types of transformation are performed here:
    # 1. categorical variables tranformation
    # 2. normalization of data to train the model: OTU table (to interval [0, 1]) and mapping variables (dividing into the maximum value)
    #output: 4 df
    
    if option==1:
        #substitution by numbers --> option 1
        X_train_df, X_test_df=categorical2numbers(X_train_df, X_test_df)
    elif option==2:
        #create colums for each value--> option 2
        X_train_df, X_test_df=categorical2columns(X_train_df, X_test_df, n-1)
    elif option==3:
        #create columns with binary coding --> option 3
        X_train_df, X_test_df=categorical2binary(X_train_df, X_test_df, n-1)
    
    X_test_df = normalize_maximum(X_test_df)
    X_train_df = normalize_maximum(X_train_df)
    y_test_df = normalize_interval(y_test_df)
    y_train_df = normalize_interval(y_train_df)
    
    return X_train_df, X_test_df, y_train_df, y_test_df

def save_results(X_train_df, y_train_df, X_test_df, y_test_df, n_iter, variable_names, df_results):
    #activation function between input layer and predicted OTUs 
    activation_function='sigmoid'

    #train model 
    model, result = train_model_nolayers(X_train_df, y_train_df, X_test_df, y_test_df, n_iter, activation_function)
    #calculate error measurements
    mse = calculate_OTUs_mse(X_test_df, y_test_df, model)
    smape = calculate_smape(X_test_df, y_test_df, model)
    
    #add results to the df
    ind = len(df_results.index) + 1
    df_results.loc[ind] = [str(variable_names), mse, smape]
    
    return model, result, df_results




#df to save error measurements of the trained models
df_results=pd.DataFrame(columns=['variables', 'mse', 'smape'])
   
     
option=1 #option parameter regulates the type of transformation of categorical data
#load data for building direct predictor
X_train_df, y_train_df, X_test_df, y_test_df = load_subsets("./data/original_data/2otu_table2save_80.biom")
variable_names=X_train_df.columns.tolist()
n=len(X_train_df.columns) - 2 #n -> necessary for categorical attribute transformation

#substitute categorical data and transform to relative frequency
X_train_df, X_test_df, y_train_df, y_test_df = transform_data(X_train_df, X_test_df, y_train_df, y_test_df, option, n)

#number of iterations to train the model
n_iter=30000

#train model with all mapping variables and save model and results
model, result, df_results = save_results(X_train_df, y_train_df, X_test_df, y_test_df, n_iter, variable_names, df_results)
model.save("./models/mapping_variables_study/direct_predictor/metadata2OTUs_maize2_all.hdf5", overwrite=True)
df_results.to_csv('./summary.csv', sep='\t')

#train model deleting 1 mapping variable in each training cycle and save each model
X_train_df, y_train_df, X_test_df, y_test_df = load_subsets("./data/original_data/2otu_table2save_80.biom")
variable_names=X_train_df.columns.tolist()
n=len(X_train_df.columns) - 2

for variable in variable_names:
    if variable != 'X.SampleID':
        variable_names.remove(variable)
        X_train_df, X_test_df, n = delete_variable(variable_names, X_train_df, X_test_df)
        X_train_df, X_test_df, y_train_df, y_test_df = transform_data(X_train_df, X_test_df, y_train_df, y_test_df, option, n)
        model, result, df_results = save_results(X_train_df, y_train_df, X_test_df, y_test_df, n_iter, variable_names, df_results)
        print(df_results)
        file_model='./models/mapping_variables_study/direct_predictor/metadata2OTUs_maize2_-'+str(variable)+'.hdf5'
        print(file_model)
        model.save(file_model, overwrite=True)
    X_train_df, y_train_df, X_test_df, y_test_df= load_subsets("./data/original_data/2otu_table2save_80.biom")
    variable_names=X_train_df.columns.tolist()

df_results.to_csv('./summary.csv', sep='\t')

#train model deleting 2 mapping variables in each training cycle and save each model
X_train_df, y_train_df, X_test_df, y_test_df = load_subsets("./data/original_data/2otu_table2save_80.biom")
variable_names=X_train_df.columns.tolist()
n=len(X_train_df.columns) - 2
n_iter=30000

for i in range(1, len(variable_names)):
    variable1=variable_names[i]  
    variable_names1=variable_names.copy()
    variable_names1.remove(variable1)
    for j in range(i+1,len(variable_names)):
        variable2=variable_names[j]
        variable_names2=variable_names1.copy()
        variable_names2.remove(variable2)
        X_train_df, X_test_df, n = delete_variable(variable_names2, X_train_df, X_test_df)
        X_train_df, X_test_df, y_train_df, y_test_df = transform_data(X_train_df, X_test_df, y_train_df, y_test_df, option, n)
        model, result, df_results = save_results(X_train_df, y_train_df, X_test_df, y_test_df, n_iter, variable_names2, df_results)
        print(df_results)
        file_model='./models/mapping_variables_study/direct_predictor/metadata2OTUs_maize2_-'+str(variable1)+'&'+str(variable2)+'.hdf5'
        model.save(file_model, overwrite=True)
        X_train_df, y_train_df, X_test_df, y_test_df = load_subsets("./data/original_data/2otu_table2save_80.biom")
        variable_names=X_train_df.columns.tolist()

df_results.to_csv('./summary.csv', sep='\t')

#train model deleting 3 mapping variables in each training cycle and save each model
X_train_df, y_train_df, X_test_df, y_test_df = load_subsets("./data/original_data/2otu_table2save_80.biom")
variable_names=X_train_df.columns.tolist()
n=len(X_train_df.columns) - 2
n_iter=30000

for i in range(1, len(variable_names)):
    variable1=variable_names[i]  
    variable_names1=variable_names.copy()
    variable_names1.remove(variable1)
    for j in range(i+1,len(variable_names)):
        variable2=variable_names[j]
        variable_names2=variable_names1.copy()
        variable_names2.remove(variable2)
        for k in range(j+1,len(variable_names)):
            variable3=variable_names[k]
            variable_names3=variable_names2.copy()
            variable_names3.remove(variable3)
            X_train_df, X_test_df, n = delete_variable(variable_names3, X_train_df, X_test_df)
            X_train_df, X_test_df, y_train_df, y_test_df = transform_data(X_train_df, X_test_df, y_train_df, y_test_df, option, n)
            model, result, df_results = save_results(X_train_df, y_train_df, X_test_df, y_test_df, n_iter, variable_names3, df_results)
            print(df_results)
            file_model='./models/mapping_variables_study/direct_predictor/metadata2OTUs_maize2_-'+str(variable1)+'&'+str(variable2)+'&'+str(variable3)+'.hdf5'
            model.save(file_model, overwrite=True)
            X_train_df, y_train_df, X_test_df, y_test_df = load_subsets("./data/original_data/2otu_table2save_80.biom")
            variable_names=X_train_df.columns.tolist()

df_results.to_csv('./summary.csv', sep='\t')

#train model deleting 4 mapping variables in each training cycle and save each model
X_train_df, y_train_df, X_test_df, y_test_df = load_subsets("./data/original_data/2otu_table2save_80.biom")
variable_names=X_train_df.columns.tolist()
n=len(X_train_df.columns) - 2
n_iter=30000

for i in range(1, len(variable_names)):
    variable1=variable_names[i]  
    variable_names1=variable_names.copy()
    variable_names1.remove(variable1)
    for j in range(i+1,len(variable_names)):
        variable2=variable_names[j]
        variable_names2=variable_names1.copy()
        variable_names2.remove(variable2)
        for k in range(j+1,len(variable_names)):
            variable3=variable_names[k]
            variable_names3=variable_names2.copy()
            variable_names3.remove(variable3)
            for l in range(k+1, len(variable_names)):
                variable4=variable_names[l]
                variable_names4=variable_names3.copy() 
                variable_names4.remove(variable4)
                X_train_df, X_test_df, n = delete_variable(variable_names4, X_train_df, X_test_df)
                X_train_df, X_test_df, y_train_df, y_test_df = transform_data(X_train_df, X_test_df, y_train_df, y_test_df, option, n)
                model, result, df_results = save_results(X_train_df, y_train_df, X_test_df, y_test_df, n_iter, variable_names4, df_results)
                print(df_results)
                file_model='./models/mapping_variables_study/direct_predictor/metadata2OTUs_maize2_-'+str(variable1)+'&'+str(variable2)+'&'+str(variable3)+'&'+str(variable4)+'.hdf5'
                model.save(file_model, overwrite=True)
                X_train_df, y_train_df, X_test_df, y_test_df = load_subsets("./data/original_data/2otu_table2save_80.biom")
                variable_names=X_train_df.columns.tolist()

df_results.to_csv('./summary.csv', sep='\t')

#train model deleting 5 mapping variables in each training cycle and save each model
X_train_df, y_train_df, X_test_df, y_test_df = load_subsets("./data/original_data/2otu_table2save_80.biom")
variable_names=X_train_df.columns.tolist()
n=len(X_train_df.columns) - 2
n_iter=30000

for i in range(1, len(variable_names)):
    variable1=variable_names[i]  
    variable_names1=variable_names.copy()
    variable_names1.remove(variable1)
    for j in range(i+1,len(variable_names)):
        variable2=variable_names[j]
        variable_names2=variable_names1.copy()
        variable_names2.remove(variable2)
        for k in range(j+1,len(variable_names)):
            variable3=variable_names[k]
            variable_names3=variable_names2.copy()
            variable_names3.remove(variable3)
            for l in range(k+1, len(variable_names)):
                variable4=variable_names[l]
                variable_names4=variable_names3.copy() 
                variable_names4.remove(variable4)
                for m in range(l+1, len(variable_names)):
                    variable5=variable_names[m]
                    variable_names5=variable_names4.copy() 
                    variable_names5.remove(variable5)
                    X_train_df, X_test_df, n = delete_variable(variable_names5, X_train_df, X_test_df)
                    X_train_df, X_test_df, y_train_df, y_test_df = transform_data(X_train_df, X_test_df, y_train_df, y_test_df, option, n)
                    model, result, df_results = save_results(X_train_df, y_train_df, X_test_df, y_test_df, n_iter, variable_names5, df_results)
                    print(df_results)
                    file_model='./models/mapping_variables_study/direct_predictor/metadata2OTUs_maize2_-'+str(variable1)+'&'+str(variable2)+'&'+str(variable3)+'&'+str(variable4)+'&'+str(variable5)+'.hdf5'
                    model.save(file_model, overwrite=True)
                    X_train_df, y_train_df, X_test_df, y_test_df = load_subsets("./data/original_data/2otu_table2save_80.biom")
                    variable_names=X_train_df.columns.tolist()

df_results.to_csv('./summary.csv', sep='\t')
