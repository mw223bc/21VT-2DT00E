# -*- coding: utf-8 -*-
"""
Created on Thu May 20 14:29:52 2021

Constructs an Bayesian model that predicts the RBP activity from BoW vectors containing 
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
from mlFunctions import plotResult
from mlFunctions import plotCost

mode = "bayesian_regression"
ep = [] # Array with the containing the number of epoches. 1-20 in this case
lo = [] # Contains the loss from each epoch
errorDif = []
MSE_values =[]
RMSE_values = []
MAE_values = []
plotOrder = ["SEQ3\n4-mers", "SEQ3\n5-mers", "SEQ3\n6-mers", "SEQ3\n7-mers", "SEQ3\n8-mers", "SEQ5\n4-mers", "SEQ5\n5-mers", "SEQ5\n6-mers", "SEQ5\n7-mers", "SEQ5\n8-mers"]

trainprecentage = 0.66  # How much of the data that should be for training
batchSize = 20  # How many data samples each batch should contain
learningRate = 0.001
epochs = 20 # How many learning cycle the training should consist of
ep.extend(range(epochs))
wDecay = 0.1 # Value of the lambda for the L1 penalization
lamda = 0.1 # Value of the lambda for the L2 penalization
learningRateDecay = 0.95    # The precenge of which the learning rate decays

plt.rcParams.update({'font.size': 6})

def model(x, single):   # Function that sepcifies how the model calculates their predictions. 
    beta = torch.abs(b_sd).cuda() * torch.randn(x.shape).cuda() + b_mu
    alpha = torch.abs(a_sd).cuda() * torch.randn(x.shape).cuda() + a_mu
    mu = beta * x + alpha
    sigma = torch.abs(sig).cuda() * torch.abs(torch.randn(x.shape)).cuda() + sig_mean
    y = torch.abs(sigma).cuda() * torch.randn(x.shape).cuda() + mu
    if single: # If the output should be calculated as a single value
        y = torch.mean(y, 1, True).cuda()
    return y

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
        
        
        inDim = len(trSet[0]) + 1 # Number of input features
        outDim = 1  # Number of output features
        
        # Creates the parameters of which should be improved after each learning cycle
        b_mu = Variable(torch.Tensor([2.0]).cuda(), requires_grad=True)
        b_sd = Variable(torch.Tensor([1.0]).cuda(), requires_grad=True)
        a_mu = Variable(torch.Tensor([0.1]).cuda(), requires_grad=True)
        a_sd = Variable(torch.Tensor([1.0]).cuda(), requires_grad=True)
        sig = Variable(torch.Tensor([1.0]).cuda(), requires_grad=True)
        sig_mean = Variable(torch.Tensor([1.0]).cuda(), requires_grad=True)
        
        optimizer = torch.optim.SGD(params=[b_mu, b_sd, a_mu, a_sd, sig, sig_mean], lr=learningRate, weight_decay=wDecay)
        criterion = torch.nn.MSELoss()
        
        print("Start model training")
        torch.cuda.empty_cache()
        for epoch in range(epochs):
            for batch in range(int(trainingSize/batchSize) + round(((trainingSize%batchSize)/batchSize) + 0.5)):
                optimizer.zero_grad()
                x, y = process(batch*batchSize, (batch+1)*batchSize, trSet, proteinsLines)
                x = torch.cuda.FloatTensor(x)
                y = torch.cuda.FloatTensor(y)
                
                pred = model(x, True)
                reg_loss = criterion(pred,y)
                regularization_loss = 0
                for g in optimizer.param_groups:    # Calculates the L1 penalization
                    for a in g['params']:
                        regularization_loss += abs(a)
                loss = reg_loss + (lamda * regularization_loss)
                
                loss.backward()
                optimizer.step()
                if not batch%10:
                    print("Batch: " + str(batch) + " from epoch " + str(epoch) + " done")
                
            for g in optimizer.param_groups:    # Reduces the learning rate
                g['lr'] *= learningRateDecay
                
                print('epoch {}, loss {}'.format(epoch, loss.item()))
                lo.append(loss.item())
    
        print("Model training complete")
            
        with torch.no_grad(): # Don't need gradients in the testing
            preds = np.array([])
            targets = np.array([])
            for batch in range(int(testSize/batchSize) + round(((testSize%batchSize)/batchSize) + 0.5)):
                x_te, y_te = process(batch*batchSize, (batch+1)*batchSize, teSet, proteinsLines)
                x_te = torch.cuda.FloatTensor(x_te)
                y_te = torch.cuda.FloatTensor(y_te)   
                preds = np.append(preds, model(Variable(x_te.cuda()), True).cpu().data.numpy())   
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
        print(errorDif)
        # LR, error acc, MSE RMSE MAE diff, LIME
        # Initialise the subplot function using number of rows and columns
        
        plotCost(plt, ep, lo, mode, seq, k)
        lo.clear()  # Clears the loss array for the next model training
        
        print("Calculating lime values")
        limeFile = open(limePath, "w")
        
        kmerImportance = np.zeros(inDim-1)
        for count, l in enumerate(bowsLines): 
            kmerNumber = [x+2 for x in l]
            kmerImportance += ((model(torch.cuda.FloatTensor(kmerNumber), False).cpu().detach().numpy())/2)**2
            if not count%100:
                print(str(count) + " of " + str(len(bowsLines)) + " done")
    
        for v in range(len(kmerImportance)):
            limeFile.write(featureNames[v] + ": " + str(kmerImportance[v]/len(bowsLines)) + "\n")
        limeFile.close()
        
        print("Done with seq" + str(seq) + " " + str(k) + "-mer")

plotResult(plt, plotOrder, errorDif, mode, MSE_values, RMSE_values, MAE_values)