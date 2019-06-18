#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 11:58:03 2019

@author: saracabellopinedo
"""

#importing necessary modules
from numpy.random import seed
seed(1)
from tensorflow import set_random_seed
set_random_seed(2)

import pandas as pd
import numpy as np
from keras.layers import Input, Dense
from keras.models import Model
from keras import regularizers
import matplotlib.pyplot as plt


def normalize_interval(data):
    #function for normalizing OTU table within [0, 1] interval
    #input: df with original values
    #output: df with normalized values
    normalized_data=pd.DataFrame(index=data.columns[1:len(data.columns)+1])
    
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
            
        normalized_data[id_sample]=new_values
    
    return normalized_data.T

def train_test_split(data):
    #split data into train and test subsets, and write both subsets to files
    #input: df with normalized data
    #otuput: two lists of lists, one for each subset
    data_train, data_test = np.split(data.sample(frac=1), [int(.8*len(data))])
    data_train.to_csv('./data/data_train_autoencoder/train_subset_norm.csv', sep='\t')
    data_test.to_csv('./data/data_train_autoencoder/test_subset_norm.csv', sep='\t')
        
    #adapting input data to the format required to keras models trining: list of lists
    list_train_maize=[]
    for row in data_train.itertuples(index=True, name='Pandas'):
        list_train_maize.append(row[1:])
    
    list_test_maize=[]
    for row in data_test.itertuples(index=True, name='Pandas'):
        list_test_maize.append(row[1:])
    
    list_train_maize=np.asarray(list_train_maize)
    list_test_maize=np.asarray(list_test_maize)
    
    return list_train_maize, list_test_maize

def train_model_arch0(list_train_maize, list_test_maize, encoding_dim, n_iter):
    #train models from architecture 0
    #takes as argument the size of the coding layer and the number of iterations
    #returns the model itself and the values obtained during the training cycle
    
    #model definition
    original_dim=len(list_train_maize[0])

    #input placeholder
    input_layer= Input(shape=(original_dim,))
    
    #encoded representation of the input
    encoded = Dense (encoding_dim, activation='tanh')(input_layer)
    
    #recontruction of the input from the encoded version
    decoded = Dense (original_dim, activation='sigmoid')(encoded)
    
    #model that maps an input to its reconstruction --> whole autoencoder
    autoencoder = Model (input_layer, decoded)
    autoencoder.compile(optimizer='adam', loss='mse')
    
    #setting random seed to make the assays reproducible
    seed(1)
    set_random_seed(2)
    history=autoencoder.fit(list_train_maize, list_train_maize,
                    nb_epoch=n_iter,
                    batch_size=256,
                    shuffle=False,
                    validation_data=(list_test_maize, list_test_maize))
    return autoencoder, history.history

def train_model_arch1(list_train_maize, list_test_maize, encoding_dim2, encoding_dim1, n_iter):
    #train models from architecture 1
    #takes as argument the size of the different hidden layers included in the model and the number of iterations
    #returns the model itself and the values obtained during the training cycle
    
    #model definition
    original_dim=len(list_train_maize[0])

    #input placeholder
    input_layer= Input(shape=(original_dim,))
    
    #encoded representation of the input
    encoded1 = Dense (encoding_dim2, activation='tanh')(input_layer)
    encoded2 = Dense (encoding_dim1, activation='tanh')(encoded1)
    
    #recontruction of the input from the encoded version
    decoded = Dense (original_dim, activation='sigmoid')(encoded2)
    
    #model that maps an input to its reconstruction --> whole autoencoder
    autoencoder = Model (input_layer, decoded)
    autoencoder.compile(optimizer='adam', loss='mse')
    
    seed(1)
    set_random_seed(2)
    history=autoencoder.fit(list_train_maize, list_train_maize,
                    nb_epoch=n_iter,
                    batch_size=256,
                    shuffle=False,
                    validation_data=(list_test_maize, list_test_maize))
    return autoencoder, history.history

def train_model_arch2(list_train_maize, list_test_maize, encoding_dim2, encoding_dim1, encoding_dim3, n_iter):
    #train models from architecture 2
    #takes as argument the size of the different hidden layers included in the model and the number of iterations
    #returns the model itself and the values obtained during the training cycle
    
    #model definition
    original_dim=len(list_train_maize[0])

    #input placeholder
    input_layer= Input(shape=(original_dim,))
    
    #encoded representation of the input
    encoded1 = Dense (encoding_dim2, activation='tanh')(input_layer)
    encoded2 = Dense (encoding_dim1, activation='tanh')(encoded1)
    
    #recontruction of the input from the encoded version
    decoded1 = Dense (encoding_dim3, activation='sigmoid')(encoded2)
    decoded2 = Dense (original_dim, activation='sigmoid')(decoded1)
    
    #model that maps an input to its reconstruction --> whole autoencoder
    autoencoder = Model (input_layer, decoded2)
    autoencoder.compile(optimizer='adam', loss='mse')


    seed(1)
    set_random_seed(2)
    history=autoencoder.fit(list_train_maize, list_train_maize,
                    nb_epoch=n_iter,
                    batch_size=256,
                    shuffle=False,
                    validation_data=(list_test_maize, list_test_maize))
    return autoencoder, history.history

def train_model_arch3(list_train_maize, list_test_maize, encoding_dim4, encoding_dim2, encoding_dim1, encoding_dim3, n_iter):
    #train models from architecture 3
    #takes as argument the size of the different hidden layers included in the model and the number of iterations
    #returns the model itself and the values obtained during the training cycle
    
    #model definition
    original_dim=len(list_train_maize[0])

    #input placeholder
    input_layer= Input(shape=(original_dim,))
    
    #encoded representation of the input
    encoded1 = Dense (encoding_dim4, activation='tanh')(input_layer)
    encoded2 = Dense (encoding_dim2, activation='tanh')(encoded1)
    encoded3 = Dense (encoding_dim1, activation='tanh')(encoded2)
    
    #recontruction of the input from the encoded version
    decoded1 = Dense (encoding_dim3, activation='sigmoid')(encoded3)
    decoded2 = Dense (original_dim, activation='sigmoid')(decoded1)
    
    #model that maps an input to its reconstruction --> whole autoencoder
    autoencoder = Model (input_layer, decoded2)
    autoencoder.compile(optimizer='adam', loss='mse')
    
#    n_iter=10000
    seed(1)
    set_random_seed(2)
    history=autoencoder.fit(list_train_maize, list_train_maize,
                    nb_epoch=n_iter,
                    batch_size=256,
                    shuffle=False,
                    validation_data=(list_test_maize, list_test_maize))
    return autoencoder, history.history

def train_model_arch4(list_train_maize, list_test_maize, encoding_dim4, encoding_dim2, encoding_dim1, encoding_dim3, encoding_dim5, n_iter):
    #train models from architecture 4
    #takes as argument the size of the different hidden layers included in the model and the number of iterations
    #returns the model itself and the values obtained during the training cycle
    
    #model definition
    original_dim=len(list_train_maize[0])

    #input placeholder
    input_layer= Input(shape=(original_dim,))
    
    #encoded representation of the input
    encoded1 = Dense (encoding_dim4, activation='tanh')(input_layer)
    encoded2 = Dense (encoding_dim2, activation='tanh')(encoded1)
    encoded3 = Dense (encoding_dim1, activation='tanh')(encoded2)
    
    #recontruction of the input from the encoded version
    decoded1 = Dense (encoding_dim3, activation='sigmoid')(encoded3)
    decoded2 = Dense (encoding_dim5, activation='sigmoid')(decoded1)
    decoded3 = Dense (original_dim, activation='sigmoid')(decoded2)
    
    #model that maps an input to its reconstruction --> whole autoencoder
    autoencoder = Model (input_layer, decoded3)
    autoencoder.compile(optimizer='adam', loss='mse')
    
#    n_iter=10000
    seed(1)
    set_random_seed(2)
    history=autoencoder.fit(list_train_maize, list_train_maize,
                    nb_epoch=n_iter,
                    batch_size=256,
                    shuffle=False,
                    validation_data=(list_test_maize, list_test_maize))
    return autoencoder, history.history

def train_model_arch5(list_train_maize, list_test_maize, encoding_dim2, encoding_dim1, encoding_dim3, encoding_dim5, n_iter):
    #train models from architecture 5
    #takes as argument the size of the different hidden layers included in the model and the number of iterations
    #returns the model itself and the values obtained during the training cycle
    
    #model definition
    original_dim=len(list_train_maize[0])

    #input placeholder
    input_layer= Input(shape=(original_dim,))
    
    #encoded representation of the input
    encoded2 = Dense (encoding_dim2, activation='tanh')(input_layer)
    encoded3 = Dense (encoding_dim1, activation='tanh')(encoded2)
    
    #recontruction of the input from the encoded version
    decoded1 = Dense (encoding_dim3, activation='sigmoid')(encoded3)
    decoded2 = Dense (encoding_dim5, activation='sigmoid')(decoded1)
    decoded3 = Dense (original_dim, activation='sigmoid')(decoded2)
    
    #model that maps an input to its reconstruction --> whole autoencoder
    autoencoder = Model (input_layer, decoded3)
    autoencoder.compile(optimizer='adam', loss='mse')
    
#    n_iter=10000
    seed(1)
    set_random_seed(2)
    history=autoencoder.fit(list_train_maize, list_train_maize,
                    nb_epoch=n_iter,
                    batch_size=256,
                    shuffle=False,
                    validation_data=(list_test_maize, list_test_maize))
    return autoencoder, history.history

def train_model_arch2_sparse(list_train_maize, list_test_maize, encoding_dim2, encoding_dim1, encoding_dim3, n_iter):
    #similar to train_model_arch2() but adding sparsity to the first layer
    
    #model definition
    original_dim=len(list_train_maize[0])

    #input placeholder
    input_layer= Input(shape=(original_dim,))
    
    #encoded representation of the input
    encoded1 = Dense (encoding_dim2, activation='tanh', activity_regularizer=regularizers.l1(0.001))(input_layer) #sparsity term
    encoded2 = Dense (encoding_dim1, activation='tanh')(encoded1)
    
    #recontruction of the input from the encoded version
    decoded1 = Dense (encoding_dim3, activation='sigmoid')(encoded2)
    decoded2 = Dense (original_dim, activation='sigmoid')(decoded1)
    
    #model that maps an input to its reconstruction --> whole autoencoder
    autoencoder = Model (input_layer, decoded2)
    autoencoder.compile(optimizer='adam', loss='mse')
    
#    n_iter=10000
    seed(1)
    set_random_seed(2)
    history=autoencoder.fit(list_train_maize, list_train_maize,
                    nb_epoch=n_iter,
                    batch_size=256,
                    shuffle=False,
                    validation_data=(list_test_maize, list_test_maize))
    return autoencoder, history.history

def train_model_arch3_sparse(list_train_maize, list_test_maize, encoding_dim4, encoding_dim2, encoding_dim1, encoding_dim3, n_iter):
    #similar to train_model_arch3() but adding sparsity to the first layer
    
    #model definition
    original_dim=len(list_train_maize[0])

    #input placeholder
    input_layer= Input(shape=(original_dim,))
    
    #encoded representation of the input
    encoded1 = Dense (encoding_dim4, activation='tanh', activity_regularizer=regularizers.l1(0.001))(input_layer) #sparsity term
    encoded2 = Dense (encoding_dim2, activation='tanh')(encoded1)
    encoded3 = Dense (encoding_dim1, activation='tanh')(encoded2)
    
    #recontruction of the input from the encoded version
    decoded1 = Dense (encoding_dim3, activation='sigmoid')(encoded3)
    decoded2 = Dense (original_dim, activation='sigmoid')(decoded1)
    
    #model that maps an input to its reconstruction --> whole autoencoder
    autoencoder = Model (input_layer, decoded2)
    autoencoder.compile(optimizer='adam', loss='mse')
    
#    n_iter=10000
    seed(1)
    set_random_seed(2)
    history=autoencoder.fit(list_train_maize, list_train_maize,
                    nb_epoch=n_iter,
                    batch_size=256,
                    shuffle=False,
                    validation_data=(list_test_maize, list_test_maize))
    return autoencoder, history.history

def train_model_arch4_sparse(list_train_maize, list_test_maize, encoding_dim4, encoding_dim2, encoding_dim1, encoding_dim3, encoding_dim5, n_iter):
    #similar to train_model_arch4() but adding sparsity to the first layer
    
    #model definition
    original_dim=len(list_train_maize[0])

    #input placeholder
    input_layer= Input(shape=(original_dim,))
    
    #encoded representation of the input
    encoded1 = Dense (encoding_dim4, activation='tanh', activity_regularizer=regularizers.l1(0.001))(input_layer) #sparsity term
    encoded2 = Dense (encoding_dim2, activation='tanh')(encoded1)
    encoded3 = Dense (encoding_dim1, activation='tanh')(encoded2)
    
    #recontruction of the input from the encoded version
    decoded1 = Dense (encoding_dim3, activation='sigmoid')(encoded3)
    decoded2 = Dense (encoding_dim5, activation='sigmoid')(decoded1)
    decoded3 = Dense (original_dim, activation='sigmoid')(decoded2)
    
    #model that maps an input to its reconstruction --> whole autoencoder
    autoencoder = Model (input_layer, decoded3)
    autoencoder.compile(optimizer='adam', loss='mse')
    
#    n_iter=10000
    seed(1)
    set_random_seed(2)
    history=autoencoder.fit(list_train_maize, list_train_maize,
                    nb_epoch=n_iter,
                    batch_size=256,
                    shuffle=False,
                    validation_data=(list_test_maize, list_test_maize))
    return autoencoder, history.history


def train_model_arch5_sparse(list_train_maize, list_test_maize, encoding_dim2, encoding_dim1, encoding_dim3, encoding_dim5, n_iter):
    #similar to train_model_arch5() but adding sparsity to the first layer
    
    #model definition
    original_dim=len(list_train_maize[0])

    #input placeholder
    input_layer= Input(shape=(original_dim,))
    
    #encoded representation of the input
    encoded2 = Dense (encoding_dim2, activation='tanh', activity_regularizer=regularizers.l1(0.001))(input_layer) #sparsity term
    encoded3 = Dense (encoding_dim1, activation='tanh')(encoded2)
    
    #recontruction of the input from the encoded version
    decoded1 = Dense (encoding_dim3, activation='sigmoid')(encoded3)
    decoded2 = Dense (encoding_dim5, activation='sigmoid')(decoded1)
    decoded3 = Dense (original_dim, activation='sigmoid')(decoded2)
    
    #model that maps an input to its reconstruction --> whole autoencoder
    autoencoder = Model (input_layer, decoded3)
    autoencoder.compile(optimizer='adam', loss='mse')
    
#    n_iter=10000
    seed(1)
    set_random_seed(2)
    history=autoencoder.fit(list_train_maize, list_train_maize,
                    nb_epoch=n_iter,
                    batch_size=256,
                    shuffle=False,
                    validation_data=(list_test_maize, list_test_maize))
    return autoencoder, history.history


#loading data and formatting it to fit the requirements of the functions
data_maize = pd.read_csv("./data/original_data/otu_table2use_80.csv", sep='\t')
data_maize.index=data_maize['otuids'].tolist()
data_maize=data_maize.drop('otuids', axis=1)
data_maize=data_maize.T
data_maize.insert(loc=0, column='X.SampleID', value=data_maize.index)

data_norm=normalize_interval(data_maize)

list_train, list_test = train_test_split(data_norm)

#examples of calling the different functions
autoencoder20_3, history20_3 = train_model_arch1(list_train, list_test, 20, 3, 10000)
autoencoder50_6_20, history50_6_20 = train_model_arch2(list_train, list_test, 50, 6, 20, 10000)
autoencoder50_10_5_100, history50_10_5_100= train_model_arch3(list_train, list_test, 50, 10, 5, 100, 10000)
autoencoder50_10_5_100_1000, history50_10_5_100_1000= train_model_arch4(list_train, list_test, 50, 10, 5, 100, 1000, 794)
autoencoder10_6_50_250, history10_6_50_250= train_model_arch5(list_train, list_test, 10, 6, 50, 250, 831)
autoencoder10_6_50s, history10_6_50s= train_model_arch2_sparse(list_train, list_test, 10, 6, 50, 12865)
autoencoder50_10_6_100s, history50_10_6_100s= train_model_arch3_sparse(list_train, list_test, 50, 10, 6, 100, 8553)
autoencoder50_10_6_50_500s, history50_10_6_50_500s= train_model_arch4_sparse(list_train, list_test, 50, 10, 6, 50, 500, 18934)
autoencoder10_6_100_1000s, history10_6_100_1000s= train_model_arch5_sparse(list_train, list_test, 10, 6, 100, 1000, 10378)

