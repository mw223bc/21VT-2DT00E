# -*- coding: utf-8 -*-
"""
Created on Fri May  7 13:05:33 2021

Creates and stores all the possible k-mer substrings for all values of k between 4 and 8

@author: Mattias Wassbjer
"""

nuc = ['A','C','G','T'] # All types of nucleotides
index = 0

for k in range (4,9):   # Loops through all the values of k
    writefilepath = "UTR data/vocabulary/vocabulary_kmers-" + str(k) + ".txt"
    writeFile = open(writefilepath, "w")
    next = [0] * k  # Indexes of the next sequence
    for i in range(4**k):
        word = ""
        while 4 in next:    # As long as 4 is found in 'next' array, increase the upcoming index by 1 and reset the number 4 to 0
            index = next.index(4)
            next[index-1] += 1
            next[index] = 0 
        for j in range(len(next)):  # Construct the nucleotide sequence acording to the 'next' array
            word += nuc[next[j]]
        writeFile.write(word + "\n")
        next[len(next)-1] += 1
                
    print("Done with vocabulary\tk: " + str(k) + "\tvocabulary size: " + str(4**k))