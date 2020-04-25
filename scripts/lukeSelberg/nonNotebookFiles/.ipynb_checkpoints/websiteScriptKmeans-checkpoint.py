import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import sys
import matplotlib.patches as mpatches
import warnings

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import MeanShift

# REPLACE these paths for your computer. this data trains the model.
path1 = "../data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta"
path01 = "../data/mycovirus_genbank_all_refseq_nucleotide_unique.fasta"

# Contains four models:

# kmeans with PCA visulization
def kmeansPCA(fastaList):
    # read in fasta files
    virus1 = parseFasta(path1)
    virus01 = parseFasta(path01)
    virus01 = virus01.append(virus1)
    virus01 = virus01.drop_duplicates(keep="last")

    # make kmer tables
    kmer7Table1 = kmerXTable(virus1, 7,7)
    kmer7Table01 = kmerXTable(virus01, 7,7)

    # only use columns that confirmed virus killers have no zeros in for kmer length 7
    cols = kmer7Table1.loc[:, (kmer7Table1 == 0).any(axis=0) != True].columns

    #Kmeans
    km = KMeans(random_state = 42, n_clusters = 2)
    km.fit(kmer7Table01[cols])
    
    output = []
    for i in range(0, len(fastaList)):
        inputData = parseFasta(fastaList[i])
        inputData["Sequence"] = inputData["Sequence"].apply(lambda x: x.replace("-", ""))
        kmer7TableInput = kmerXTable(inputData, 7,7)
        
        # generate array of labels
        y_hat = km.predict(kmer7TableInput[cols])

        #PCA
        if len(inputData) > 1:
            embedding = PCA()
            embedding.fit(kmer7TableInput[cols])
            show = pd.DataFrame(embedding.transform(kmer7TableInput[cols]))
            # show kmeans clustering
            show["labels"] = y_hat
            ax = show[show["labels"]==1].plot.scatter(x=0, y=1, style="o", color="red", s=2)
            show[show["labels"]==0].plot.scatter(x=0, y=1, style="o", color="blue", s=2, ax=ax)
            red = mpatches.Patch(color='red', label='Fungus Killers')
            blue = mpatches.Patch(color='blue', label='Fungus Non-Killers')
            plt.legend(handles=[red, blue])
            plt.title('PCA Visualization for KMeans Clustering\n' + fastaList[i])
            plt.xlabel('First Principal Component')
            plt.ylabel('Second Principal Component')
            plt.savefig('kmeansPCA' + str(i) + '.png', bbox_inches='tight')
            plt.close()

        inputData["Labels"] = y_hat
        output.append(inputData)
    return output

# kmeans with tSNE visulization
def kmeansTSNE(fastaList):
    # read in fasta files
    virus1 = parseFasta(path1)
    virus01 = parseFasta(path01)
    virus01 = virus01.append(virus1)
    virus01 = virus01.drop_duplicates(keep="last")

    # make kmer tables
    kmer7Table1 = kmerXTable(virus1, 7,7)
    kmer7Table01 = kmerXTable(virus01, 7,7)

    # only use columns that confirmed virus killers have no zeros in for kmer length 7
    cols = kmer7Table1.loc[:, (kmer7Table1 == 0).any(axis=0) != True].columns

    #Kmeans
    km = KMeans(random_state = 42, n_clusters = 2)
    km.fit(kmer7Table01[cols])
    
    output = []
    for i in range(0, len(fastaList)):
        inputData = parseFasta(fastaList[i])
        inputData["Sequence"] = inputData["Sequence"].apply(lambda x: x.replace("-", ""))
        kmer7TableInput = kmerXTable(inputData, 7,7)
        
        # generate array of labels
        y_hat = km.predict(kmer7TableInput[cols])

        #tSNE
        if len(inputData) > 1:
            tSNEembedding = TSNE(n_components= 2, random_state = 0)
            tSNEembedding_low = tSNEembedding.fit_transform(kmer7TableInput[cols])
            show = pd.DataFrame(tSNEembedding_low)
            # show kmeans clustering
            show["labels"] = y_hat
            ax = show[show["labels"]==1].plot.scatter(x=0, y=1, style="o", color="red", s=2)
            show[show["labels"]==0].plot.scatter(x=0, y=1, style="o", color="blue", s=2, ax=ax)
            red = mpatches.Patch(color='red', label='Fungus Killers')
            blue = mpatches.Patch(color='blue', label='Fungus Non-Killers')
            plt.legend(handles=[red, blue])
            plt.title('tSNE Visualization for KMeans Clustering\n' + fastaList[i])
            plt.xlabel('First Component')
            plt.ylabel('Second Component')
            plt.savefig('kmeansTSNE' + str(i) + '.png', bbox_inches='tight')
            plt.close()

        inputData["Labels"] = y_hat
        output.append(inputData)
    return output

# kmeans and meanshift with PCA visulization
def kmeansMeanshiftPCA(fastaList):
    # read in fasta files
    virus1 = parseFasta(path1)
    virus01 = parseFasta(path01)
    virus01 = virus01.append(virus1)
    virus01 = virus01.drop_duplicates(keep="last")

    # make kmer tables
    kmer7Table1 = kmerXTable(virus1, 7,7)
    kmer7Table01 = kmerXTable(virus01, 7,7)

    # only use columns that confirmed virus killers have no zeros in for kmer length 7
    cols = kmer7Table1.loc[:, (kmer7Table1 == 0).any(axis=0) != True].columns
    
    #meanshift
    ms = MeanShift()
    embedding = PCA()
    ms.fit(pd.DataFrame(embedding.fit_transform(kmer7Table01[cols])))
    cluster_centers = ms.cluster_centers_

    #Kmeans
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        km = KMeans(init = cluster_centers, n_clusters = len(cluster_centers))
        km.fit(kmer7Table01[cols])
    newLabelClusters = km.predict(kmer7Table01[cols])[343:]
    
    output = []
    for i in range(0, len(fastaList)):
        inputData = parseFasta(fastaList[i])
        inputData["Sequence"] = inputData["Sequence"].apply(lambda x: x.replace("-", ""))
        kmer7TableInput = kmerXTable(inputData, 7,7)
        
        # generate array of labels
        y_hat = km.predict(kmer7TableInput[cols])
        newLabels = []
        for label in y_hat:
            if label in newLabelClusters:
                newLabels.append(1)
            else:
                newLabels.append(0)  

        #PCA
        if len(inputData) > 1:
            embedding = PCA()
            embedding.fit(kmer7TableInput[cols])
            show = pd.DataFrame(embedding.transform(kmer7TableInput[cols]))
            # show kmeans clustering
            show["labels"] = newLabels
            ax = show[show["labels"]==1].plot.scatter(x=0, y=1, style="o", color="red", s=2)
            show[show["labels"]==0].plot.scatter(x=0, y=1, style="o", color="blue", s=2, ax=ax)
            red = mpatches.Patch(color='red', label='Fungus Killers')
            blue = mpatches.Patch(color='blue', label='Fungus Non-Killers')
            plt.legend(handles=[red, blue])
            plt.title('PCA Visualization for KMeans Clustering with Mean Shift\n' + fastaList[i])
            plt.xlabel('First Principal Component')
            plt.ylabel('Second Principal Component')
            plt.savefig('kmeansMeanshiftPCA' + str(i) + '.png', bbox_inches='tight')
            plt.close()

        inputData["Labels"] = newLabels
        output.append(inputData)
    return output

# kmeans and meanshift with tSNE visulization
def kmeansMeanshiftTSNE(fastaList):
    # read in fasta files
    virus1 = parseFasta(path1)
    virus01 = parseFasta(path01)
    virus01 = virus01.append(virus1)
    virus01 = virus01.drop_duplicates(keep="last")

    # make kmer tables
    kmer7Table1 = kmerXTable(virus1, 7,7)
    kmer7Table01 = kmerXTable(virus01, 7,7)

    # only use columns that confirmed virus killers have no zeros in for kmer length 7
    cols = kmer7Table1.loc[:, (kmer7Table1 == 0).any(axis=0) != True].columns
    
    #meanshift
    ms = MeanShift()
    embedding = PCA()
    ms.fit(pd.DataFrame(embedding.fit_transform(kmer7Table01[cols])))
    cluster_centers = ms.cluster_centers_

    #Kmeans
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        km = KMeans(init = cluster_centers, n_clusters = len(cluster_centers))
        km.fit(kmer7Table01[cols])
    newLabelClusters = km.predict(kmer7Table01[cols])[343:]
    
    output = []
    for i in range(0, len(fastaList)):
        inputData = parseFasta(fastaList[i])
        inputData["Sequence"] = inputData["Sequence"].apply(lambda x: x.replace("-", ""))
        kmer7TableInput = kmerXTable(inputData, 7,7)
        
        # generate array of labels
        y_hat = km.predict(kmer7TableInput[cols])
        newLabels = []
        for label in y_hat:
            if label in newLabelClusters:
                newLabels.append(1)
            else:
                newLabels.append(0)  

        #PCA
        if len(inputData) > 1:
            tSNEembedding = TSNE(n_components= 2, random_state = 0)
            tSNEembedding_low = tSNEembedding.fit_transform(kmer7TableInput[cols])
            show = pd.DataFrame(tSNEembedding_low)
            # show kmeans clustering
            show["labels"] = newLabels
            ax = show[show["labels"]==1].plot.scatter(x=0, y=1, style="o", color="red", s=2)
            show[show["labels"]==0].plot.scatter(x=0, y=1, style="o", color="blue", s=2, ax=ax)
            red = mpatches.Patch(color='red', label='Fungus Killers')
            blue = mpatches.Patch(color='blue', label='Fungus Non-Killers')
            plt.legend(handles=[red, blue])
            plt.title('tSNE Visualization for KMeans Clustering with Mean Shift\n' + fastaList[i])
            plt.xlabel('First Component')
            plt.ylabel('Second Component')
            plt.savefig('kmeansMeanshiftTSNE' + str(i) + '.png', bbox_inches='tight')
            plt.close()

        inputData["Labels"] = newLabels
        output.append(inputData)
    return output

# Parse datafile. Contains an index with sequence names and a column with sequences
# Pass the fasta file as the data parameter
# Example: virus1 = parseFasta("data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta")
def parseFasta(data):
    d = {fasta.id : str(fasta.seq) for fasta in SeqIO.parse(data, "fasta")}
    pd.DataFrame([d])

    s = pd.Series(d, name='Sequence')
    s.index.name = 'ID'
    s.reset_index()
    return pd.DataFrame(s)

# Make kmer table with kmers of length a to b
# Pass a dataframe from parseFasta method as the s parameter, and the min to max kmer range as a and b respectively
# Example: kmer7Table1 = kmerXTable(virus1, 7,7)
def kmerXTable(s, a, b):
    tfid_vector = TfidfVectorizer(analyzer='char', ngram_range=(a,b))
    s_hat = tfid_vector.fit_transform(s.Sequence)
    kmerNames = tfid_vector.get_feature_names()
    kmers = s_hat.toarray()
    return pd.DataFrame(kmers,columns=kmerNames, index = s.index)

## Below shows example usage of this script. Paths are relative to my computer.
## Uncomment below to check functionality if file paths are given in the code

# x = kmeansPCA(["../data/mycovirus_genbank_all_refseq_nucleotide_unique.fasta",  "../data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta"])
# for df in x:
#     print(df.head())
    
# x = kmeansTSNE(["../data/mycovirus_genbank_all_refseq_nucleotide_unique.fasta",  "../data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta"])
# for df in x:
#     print(df.head())
    
# x = kmeansMeanshiftPCA(["../data/mycovirus_genbank_all_refseq_nucleotide_unique.fasta",  "../data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta"])
# for df in x:
#     print(df.head())
    
# x = kmeansMeanshiftTSNE(["../data/mycovirus_genbank_all_refseq_nucleotide_unique.fasta",  "../data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta"])
# for df in x:
#     print(df.head())
    

###########################
# Old model:
# Pass a list of fasta file locations
# Example: websiteScriptKmeans(["../data/mycovirus_genbank_all_refseq_nucleotide_unique.fasta",  "../data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta"])
def websiteScriptKmeans(fastaList):
    # read in fasta files
    virus1 = parseFasta(path1)
    virus01 = parseFasta(path01)
    virus01 = virus01.append(virus1)
    virus01 = virus01.drop_duplicates(keep="last")

    # make kmer tables
    kmer7Table1 = kmerXTable(virus1, 7,7)
    kmer7Table01 = kmerXTable(virus01, 7,7)

    # only use columns that confirmed virus killers have no zeros in for kmer length 7
    cols = kmer7Table1.loc[:, (kmer7Table1 == 0).any(axis=0) != True].columns

    #Kmeans
    km = KMeans(random_state = 42, n_clusters = 2)
    km.fit(kmer7Table01[cols])
    
    output = []
    for i in range(0, len(fastaList)):
        inputData = parseFasta(fastaList[i])
        inputData["Sequence"] = inputData["Sequence"].apply(lambda x: x.replace("-", ""))
        kmer7TableInput = kmerXTable(inputData, 7,7)
        
        # generate array of labels
        y_hat = km.predict(kmer7TableInput[cols])

        #PCA
        if len(inputData) > 1:
            embedding = PCA()
            embedding.fit(kmer7TableInput[cols])
            show = pd.DataFrame(embedding.transform(kmer7TableInput[cols]))
            # show kmeans clustering
            ax = show.plot.scatter(x=0, y=1, style="o", c=y_hat, cmap = "viridis", s=2)
            plt.title('PCA Visualization for KMeans Clustering')
            plt.xlabel('First Principal Component')
            plt.ylabel('Second Principal Component')
            plt.savefig('visual' + str(i) + '.png', bbox_inches='tight')
            plt.close()

        inputData["Labels"] = y_hat
        output.append(inputData)
    return output

## Uncomment below to check functionality if file paths are given as command line arguements
## Example: python websiteScriptKmeans.py ../data/mycovirus_genbank_all_refseq_nucleotide_unique.fasta ../data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta

# fastaList = sys.argv[1:]
# x = websiteScriptKmeans(fastaList)
# for df in x:
#     print(df.head())