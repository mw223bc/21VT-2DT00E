# -*- coding: utf-8 -*-
"""
Created on Thu May  6 14:30:55 2021

Finds the longest sequences and stores them in a new file.

@author: Mattias Wassbjer
"""

for j in range(2):  #Loops through both sequence files
    seq = 3 + (j*2)
    readpath = "UTR data/seq" + str(seq) + "/seq" + str(seq) + ".txt"
    writepath = "UTR data/seq" + str(seq) + "/longest/longest_UTRs_seq" + str(seq) + ".txt"
    readFile = open(readpath, "r")
    
    lines = readFile.readlines()
    
    longest_UTRs = ["Empty"]
    numberOfUTRs = 0
    
    for l in lines: # Finds the largest sequences
        a = l.split("|")
        if a[0] == longest_UTRs[numberOfUTRs].split("|")[0]:
            longest_UTRs[numberOfUTRs] = max(l, longest_UTRs[numberOfUTRs] , key=len)
        else:
            longest_UTRs.append(l)
            numberOfUTRs += 1
            
    longest_UTRs.pop(0)
    print(len(longest_UTRs))
    
    writeFile = open(writepath, "w")
    for i in range(len(longest_UTRs)):
        writeFile.write(longest_UTRs[i])