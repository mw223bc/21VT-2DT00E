# -*- coding: utf-8 -*-
"""
Created on Fri May  7 14:57:35 2021

Converts the k-mers into Bag of Words vectors.

@author: Mattias Wassbjer
"""

# Loops through all sequences and values of k
for i in range(2):
    seq = 3 + (i*2)
    for k in range (4,9):
        readfilepath = "UTR data/seq" + str(seq) + "/kmer/kmers-" + str(k) + "_seq" + str(seq) + ".txt"
        vocabularypath = "UTR data/vocabulary/vocabulary_kmers-" + str(k) + ".txt"
        writefilepath = "UTR data/seq" + str(seq) + "/BoW/vector_BoW_kmers-" + str(k) + "_seq" + str(seq) + ".txt"
        readFile = open(readfilepath, "r")
        writeFile = open(writefilepath, "w")
        vocabularyFile = open(vocabularypath, "r")
        kmearLines = readFile.readlines()
        kmearLines = [s.strip("\n") for s in kmearLines]
        vocabulary = vocabularyFile.readlines()
        vocabulary = [s.strip("\n") for s in vocabulary]
        lines = 0   # Number of lines written
        print('Length of file: ' + str(len(kmearLines)))
        
        for counter, kl in enumerate(range(len(kmearLines))):
            data = ""
            BoW = [0] * len(vocabulary)
            words = kmearLines[kl].split('\t')[1].split(' ')
            for w in words: # Count the number of times all substrings occurs
                BoW[vocabulary.index(w)] += 1
                
            for j in range(len(BoW)):
                data += str(BoW[j]) + " "   # Converts the BoW array to a data string
            data = data[:-1] + "\n"
            writeFile.write(data)
            lines += 1
            if not counter%100:
                print(str(counter) + " of " + str(len(kmearLines)) + " done")
        print("Done with BoW \tseq: " + str(seq) + "\tk: " + str(k) + "\tdata size: " + str(lines))
        writeFile.close()