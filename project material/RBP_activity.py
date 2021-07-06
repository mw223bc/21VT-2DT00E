# -*- coding: utf-8 -*-
"""
Created on Mon May 17 20:39:17 2021

Finds and stores the RBP activity for each and every gene in the k-mer document

@author: mwass
"""

# Loops through both sequences and all values of k
for i in range(2):
    seq = 3 + (i*2)
    for k in range(4,9):
        kmersPath = "UTR data/seq" + str(seq) + "/kmer/kmers-" + str(k) + "_seq" + str(seq) + ".txt"
        MassSpecPath = "UTR data/MassSpec/NCOMMS-19-7936188_MassSpec_scaled_abundances.txt"
        writePath = "UTR data/seq" + str(seq) + "/RBP_activity/activity_kmers-" + str(k) + "_seq" + str(seq) + ".txt"
        kmers = open(kmersPath, "r")
        MassSpec = open(MassSpecPath, "r")
        
        kmerLines = kmers.readlines()
        kmerLines = [s.strip("\n") for s in kmerLines]
        RBPsActivity = MassSpec.readlines()
        RBPsActivity = [s.strip("\n") for s in RBPsActivity]
        
        RBPsActivity.pop(0) # Removes columns header
        
        writeFile = open(writePath, "w")
        
        for counter, km in enumerate(kmerLines):    # Find the RBP activity for all genes in the k-mer document
            activityLine = next((s for s in RBPsActivity if km.split("|")[0] in s), None)
            writeFile.write(activityLine + "\n")
            if not counter%100:
                print(str(counter) + " of " + str(len(kmerLines)) + " done")
        
        print("Seq:" + str(seq) + "\tk: " + str(k) + " done")
writeFile.close()