# -*- coding: utf-8 -*-
"""
Created on Thu May 20 14:29:52 2021

Constructs an L2 penalized or an Elastic net model that predicts the RBP activity from BoW vectors containing 
the number of times a nucleotide seqence occurs according to the order of the vocabulary. After the model has been 
trained and tested, the most predict-decisiv RBP sequences are calculated using LIME.

@author: Mattias Wassbjer
"""

import random
import torch
import numpy as np
from torch.autograd import Variable
import matplotlib.pyplot as plt
from mlFunctions import process
from mlFunctions import fit
from mlFunctions import calcLime
from mlFunctions import plotResult
from mlFunctions import plotCost


elastic = True # If the model should be elastic net or L2 penalized only
if elastic:
    mode = "elastic_net"
else:
    mode = "L2_linear_regression"
    
ep = [] # Array with the containing the number of epoches. 1-20 in this case
lo = [] # Contains the loss from each epoch
errorDif = []
MSE_values =[]
RMSE_values = []
MAE_values = []
plotOrder = ["SEQ3\n4-mers", "SEQ3\n5-mers", "SEQ3\n6-mers", "SEQ3\n7-mers", "SEQ3\n8-mers", "SEQ5\n4-mers", "SEQ5\n5-mers", "SEQ5\n6-mers", "SEQ5\n7-mers", "SEQ5\n8-mers"]

trainprecentage = 0.66  # How much of the data that should be for training
batchSize = 100 # How many data samples each batch should contain
learningRate = 0.01
epochs = 20 # How many learning cycle the training should consist of
ep.extend(range(epochs))
wDecay = 0.1 # Value of the lambda for the L1 penalization
if elastic: # If the model should use L2 penalization
    lamda = 0.1
    print("Modeling Elastic Net")
else:
    lamda = 0
    print("Modeling Linear Regression with L2")
learningRateDecay = 0.95    # The precenge of which the learning rate decays

plt.rcParams.update({'font.size': 6})

# Loops through all the sequences and values of k
for i in range(2):
    seq = 3 + (i*2)
    for k in range(4,9):
        
        bowsPath = "UTR data/seq" + str(seq) + "/BoW/vector_BoW_kmers-" + str(k) + "_seq" + str(seq) + ".txt"
        proteinsPath = "UTR data/seq" + str(seq) + "/RBP_activity/activity_kmers-" + str(k) + "_seq" + str(seq) + ".txt"
        vocabularyPath = "UTR data/vocabulary/vocabulary_kmers-" + str(k) + ".txt"
        limePath = mode + "_results/LIME_values/lime_seq"+ str(seq) + "_" + str(k) + "-mer.txt"
        bowsFile = open(bowsPath, "r")
        proteinsFile = open(proteinsPath, "r")
        vocabularyFile = open(vocabularyPath, "r")
        
        bowsLines = bowsFile.readlines()
        bowsLines = [s.strip("\n") for s in bowsLines]
        bowsLines = [s.split(" ") for s in bowsLines] 
        bowsLines = [list(map(int, s)) for s in bowsLines] 
        proteinsLines = proteinsFile.readlines()
        proteinsLines = [s.strip("\n") for s in proteinsLines]
        proteinsLines = [s.split("\t") for s in proteinsLines]
        proteinsLines = [s[3:] for s in proteinsLines]
        for p in proteinsLines:
            for i in range(len(p)):
                    if p[i] == "NA": p[i] = "-1"
        proteinsLines = [list(map(float, i)) for i in proteinsLines]
        featureNames = vocabularyFile.readlines()
        featureNames = [s.strip("\n") for s in featureNames]
        bowsFile.close()
        proteinsFile.close()
        
        print("Files read")
        
        # Randomly shuffles the data 
        shuffled = list(zip(bowsLines, proteinsLines))
        random.shuffle(shuffled)
        bowsLines, proteinsLines = zip(*shuffled)
        bowsLines, proteinsLines = list(bowsLines), list(proteinsLines)

        trainingSize = int(len(bowsLines)*trainprecentage)
        testSize = len(bowsLines) - trainingSize
        
        # Devides the data for training and testing
        trSet = bowsLines[:trainingSize]
        teSet = bowsLines[trainingSize:len(bowsLines)]
        
        inDim = len(trSet[0]) + 1   # Number of input nodes
        outDim = 1  # Number of output nodes
        
        # Creates the linear model and puts in into the CUDA
        model = torch.nn.Linear(inDim, outDim)
        model.cuda()
        
        criterion = torch.nn.MSELoss()
        optimizer = torch.optim.SGD(model.parameters(), lr=learningRate, weight_decay=wDecay)
            
        print("Start model training")
        fit(epochs, model, trainingSize, batchSize, criterion, optimizer, trSet, proteinsLines, learningRateDecay, lo, lamda)
        print("Model training complete")
            
        model.eval()
        with torch.no_grad(): # Don't need gradients in the testing 
            preds = np.array([])
            targets = np.array([])
            for batch in range(int(testSize/batchSize) + round(((testSize%batchSize)/batchSize) + 0.5)):
                x_te, y_te = process(batch*batchSize, (batch+1)*batchSize, teSet, proteinsLines)
                x_te = torch.cuda.FloatTensor(x_te)
                y_te = torch.cuda.FloatTensor(y_te)   
                preds = np.append(preds, model(Variable(x_te.cuda())).cpu().data.numpy())   
                targets = np.append(targets, Variable(y_te.cuda()).cpu().data.numpy())
    
        difference = preds - targets
        print("Evaluation metrics:")
        MSE = np.mean(difference**2)
        MAE = np.mean(abs(difference))
        RMSE = np.sqrt(MSE)
        
        print("MSE:", MSE)
        print("RMSE:", RMSE)
        print("MAE:", MAE)
        
        MSE_values.append(MSE)
        RMSE_values.append(RMSE)
        MAE_values.append(MAE)
        errorDif.append((MAE/np.mean(targets)) *100)
        
        plotCost(plt, ep, lo, mode, seq, k)
        lo.clear()  # Clears the loss array for the next model training
        
        # Extracts the values of the beta parameters
        params = list(model.parameters())[0].tolist()[0]
        params.pop()
        params = np.array(params)
        
        print("Calculating lime values")
        limeFile = open(limePath, "w")
        
        kmerImportance = np.zeros(inDim-1)
        calcLime((inDim-1), np.array(bowsLines), kmerImportance, params)
        for v in range(len(kmerImportance)):
            limeFile.write(featureNames[v] + ": " + str(kmerImportance[v]/len(bowsLines)) + "\n")
        limeFile.close()
        
        print("Done with seq" + str(seq) + " " + str(k) + "-mer")

plotResult(plt, plotOrder, errorDif, mode, MSE_values, RMSE_values, MAE_values)