# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 13:00:14 2021

A set of common functions used in the implementations of the models

@author: Mattias Wassbjer
"""
import torch
import numpy as np
from sklearn import preprocessing
from numba import jit

# Maps the RBP activity to the corresponding BoW vectors
def process(start, end, tSet, proteinsLines):
    x = []
    y = []
    if end > len(tSet): end = len(tSet)

    # Adds the condition terms to the BoW vectors, 'X', and the corresponding output RBP activiy to 'y'
    for i in range(start,end): 
        for count, a in enumerate(proteinsLines[i]):
            if a != -1:
                x.append(tSet[i] + [count])
                y.append(float(a))
                
    x = np.array(x, dtype=np.float32)
    x = preprocessing.normalize(x)  # Normalizes the input features
    y = np.array(y, dtype=np.float32)
    y = y.reshape(-1,1)
    return x, y

# Fits model onto the training data
def fit(num_epochs, model, dataSize, batchSize, criterion, opt, trSet, proteinsLines, learningRateDecay, lo, lamda):
    torch.cuda.empty_cache()    
    for epoch in range(num_epochs):
        for batch in range(int(dataSize/batchSize) + round(((dataSize%batchSize)/batchSize) + 0.5)):
            opt.zero_grad()
            x, y = process(batch*batchSize, (batch+1)*batchSize, trSet, proteinsLines)  # Processes the current batch
            x = torch.cuda.FloatTensor(x)
            y = torch.cuda.FloatTensor(y)
            
            loss = train(x,y, model, criterion, opt, lamda) # Trains the model and calculates its loss
            if not batch%10:
                print("Batch: " + str(batch) + " from epoch " + str(epoch) + " done")
            
        for g in opt.param_groups:  # Reduces the learning rate
            g['lr'] *= learningRateDecay
            
            print('epoch {}, loss {}'.format(epoch, loss.item()))
            lo.append(loss.item())
     
# Trains the model on the x and y data, and calculates the regression loss
def train(x, y, model, criterion, opt, lamda):
    pred = model(x) 
    regularization_loss = 0
    for param in model.parameters():    # Calculates the L1 penalization
        regularization_loss += torch.sum(abs(param))
 
    reg_loss = criterion(pred,y)
    
    loss = reg_loss + lamda * regularization_loss
    
    loss.backward()
    opt.step()
    return loss

# Calculates the LIME values with optimized machine code
@jit(nopython=True)
def calcLime(rge, bowsLines, kmerImportance, params):
    for a in range(rge):
            for l in bowsLines: 
                kmerNumber = l[a]
                kmerImportance[a] += ((params[a]*(kmerNumber+2))/2)**2
            if not a%100:
                    print(str(a) + " of " + str(rge) + " done")

# Plots and saves the model's learning curves
def plotCost(plt, ep, lo, mode, seq, k):
    plt.plot(ep, lo)
    plt.title(mode + " - Train Cost - Sequence " + str(seq) + ", " + str(k) + "-mer")
    plt.xlabel("Iterations")
    plt.ylabel("Cost")
    plt.savefig(mode + "_results/plots/CostPlot_seq" + str(seq) + "_k" + str(k), dpi=1200)
    plt.clf()

# Plots and saves the models' evaluation metrics
def plotResult(plt, plotOrder, errorDif, mode, MSE_values, RMSE_values, MAE_values):        
    plt.bar(plotOrder, errorDif)
    plt.title(mode + " - Mean error difference")
    plt.ylabel("Difference in %")
    plt.savefig(mode + "_results/plots/error_diff", dpi=1200)
    plt.clf()
    
    plt.bar(plotOrder, MSE_values)
    plt.title(mode + " - Mean squared error comparison")
    plt.ylabel('MSE')
    plt.savefig(mode + "_results/plots/MSE_diff", dpi=1200)
    plt.clf()
    
    plt.bar(plotOrder, RMSE_values)
    plt.title(mode + " - Root-mean squared error comparison")
    plt.ylabel('RMSE')
    plt.savefig(mode + "_results/plots/RMSE_diff", dpi=1200)
    plt.clf()
    
    plt.bar(plotOrder, MAE_values)
    plt.title(mode + " - Mean absolute error comparison")
    plt.ylabel('MAE')
    plt.savefig(mode + "_results/plots/MAE_diff", dpi=1200)
    plt.clf()