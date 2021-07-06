# -*- coding: utf-8 -*-
"""
Created on Thu May  6 17:50:59 2021

Converts the longest UTR sequences to substrings of length k

@author: Mattias Wassbjer
"""

kmers = ""
for i in range(2):  # Loops through both seqences
    seq = 3 + (i*2)
    readpath = "UTR data/seq" + str(seq) + "/longest/longest_UTRs_seq" + str(seq) + ".txt"
    noOutputsPath = "UTR data/seq" + str(seq) + "/GenesWithNoOutputData_seq" + str(seq) + ".txt"
    readFile = open(readpath, "r")
    noOutputs = open(noOutputsPath, 'r')
    lines = readFile.readlines()
    lines = [s.strip("\n") for s in lines]
    noOuts = noOutputs.readlines()
    noOuts = [s.strip("\n") for s in noOuts]
    
    for k in range(4, 9):
        print("Start of k:" + str(k) + " seq: " + str(seq))
        writepath = "UTR data/seq" + str(seq) + "/kmer/kmers-" + str(k) + "_seq" + str(seq) + ".txt"
        writeFile = open(writepath, "w")
        for counter, l in enumerate(lines):
            a = l.split("\t")   
            
            if len(a[1]) >= k and l.split("|")[0] not in noOuts:    # If the sequence is longer than k and has an output
                while len(a[1]) >= k:
                    next = a[1][0:k]
                    if 'N' not in next: # Removes k-mers that contains 'any base' N
                        kmers += next + " "
                    a[1] = a[1][1:]
                kmers = kmers[:-1]
                if len(kmers) > 0:  # If the sequences is split into one or more substrings of length k
                    writeFile.write(a[0] + "\t" + kmers + "\n")
                kmers = ""
            if not counter%1000:
                print(str(counter) + " of " + str(len(lines)) + " done")
        writeFile.close()