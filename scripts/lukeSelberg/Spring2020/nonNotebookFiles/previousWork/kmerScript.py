# Imports
import readFasta
from itertools import product
import pandas as pd
import numpy as np

# Functions 
def possibleKMERS(n):
    return [p for p in product(["A","G","T","C"], repeat = n)]

def kmerCounter(test,n):
    pk = possibleKMERS(n)

    zeros = np.zeros((len(test),len(pk)))

    for row in range(0, len(test)):
        for column in range(0, len(test.iloc[0])-n+1):
            if test.iloc[row].iloc[column] is None:
                break;
            zeroIndex = 0
            for kmer in pk:
                if (kmer == tuple(test.iloc[row].iloc[column:column+n])):
                    zeros[row][zeroIndex] += 1
                    break
                zeroIndex += 1
            
    return pd.DataFrame(columns=pk, index=test.index, data = zeros)

def frequencyMaker(count, total):
    return count/total

def kmerCountToFrequency(test):
    df = pd.DataFrame()
    testTotal = test.total
    for row in range(0, len(test)):
        df = df.append(test.iloc[row].apply(frequencyMaker, args=(test.iloc[row].total,)))
    df.total = testTotal
    return df

# Script
virus1 = readFasta.readFasta("data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta")

kmerCount = kmerCounter(virus1,2)
kmerCount["total"] = kmerCount.apply(np.sum,axis=1)

kmerCount.to_csv("kmerCount.csv")
pd.read_csv("kmerCount.csv")

kmerFreq = kmerCountToFrequency(kmerCount)
kmerFreq.to_csv("kmerfreq.csv")
pd.read_csv("kmerfreq.csv")