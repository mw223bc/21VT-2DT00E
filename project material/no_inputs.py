# -*- coding: utf-8 -*-
"""
Created on Mon May 17 23:15:16 2021

Finds genes that is listed in the longest UTR documents but has no listed RBP activity

@author: Mattias Wassbjer
"""

missing = []    # Array containing the gene ids with missing RBP activity

for i in range(2):  # Loops through both sequences
    seq = 3 + (i*2)
    largest_UTRs_Path = "UTR data/seq" + str(seq) + "/longest/longest_UTRs_seq" + str(seq) + ".txt"
    MassSpec_Path = "UTR data/MassSpec/NCOMMS-19-7936188_MassSpec_scaled_abundances.txt"
    write_Path = "UTR data/seq" + str(seq) + "/GenesWithNoOutputData_seq" + str(seq) + ".txt"
    largest_UTRs = open(largest_UTRs_Path, "r")
    MassSpec = open(MassSpec_Path, "r")
    
    UTRs = largest_UTRs.readlines()
    UTRs = [s.strip("\n") for s in UTRs]
    RBPs = MassSpec.readlines()
    RBPs = [s.strip("\n") for s in RBPs]
    
    RBPs.pop(0) # Removes column names
    writeFile = open(write_Path, "w")
    c = 0   # Number of non-listed genes
    for counter, u in enumerate(UTRs):
        if not any(u.split("|")[0] in s for s in RBPs): # If the gene isn't found in both documents
            c += 1
            writeFile.write(u.split("|")[0] + "\n")
        if not counter%1000:
            print(str(counter) + " of " + str(len(UTRs)) + " done")
            
    print(str(c))
    print("Seq:" + str(seq) + " done")
writeFile.close()