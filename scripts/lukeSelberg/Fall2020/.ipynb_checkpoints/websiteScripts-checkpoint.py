import pandas as pd
from Bio import SeqIO
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import MeanShift
from sklearn import preprocessing 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def parseFasta(data):
    d = {fasta.id : str(fasta.seq) for fasta in SeqIO.parse(data, "fasta")}
    pd.DataFrame([d])

    s = pd.Series(d, name='Sequence')
    s.index.name = 'ID'
    s.reset_index()
    return pd.DataFrame(s)

def kmerXTable(s, a, b):
    tfid_vector = TfidfVectorizer(analyzer='char', ngram_range=(a,b))
    s_hat = tfid_vector.fit_transform(s.Sequence)
    kmerNames = tfid_vector.get_feature_names()
    kmers = s_hat.toarray()
    return pd.DataFrame(kmers,columns=kmerNames, index = s.index)
    
def kmeans(fasta, klength_min, klength_max, rNum, cNum):
    inputData = parseFasta(fasta)
#     temp = virus01.append(inputData)
#     temp = temp.drop_duplicates(keep="last")
        
    inputData["Sequence"] = inputData["Sequence"].apply(lambda x: x.replace("-", ""))
    kmerXTableInput = kmerXTable(inputData, klength_min, klength_max)
        
        
    km = KMeans(random_state = rNum, n_clusters = cNum)
    km.fit(kmerXTableInput) 
    y_hat = km.predict(kmerXTableInput)
        
    return y_hat, kmerXTableInput
        
def kmeans_semiSupervised(fasta, klength_min, klength_max, rNum, y_hat):
    inputData = parseFasta(fasta)
    inputData["Sequence"] = inputData["Sequence"].apply(lambda x: x.replace("-", ""))
    kmerXTableInput = kmerXTable(inputData, klength_min, klength_max)
    
    PCAembedding = PCA(n_components=2)
    NkmerXTableInput = preprocessing.normalize(kmerXTableInput)
    PCAembedding_low = PCAembedding.fit_transform(NkmerXTableInput)
    
    ms = MeanShift()
    ms.fit(PCAembedding_low)
    cluster_centers = ms.cluster_centers_

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        kmms = KMeans(init = cluster_centers, n_clusters = len(cluster_centers))
        kmms_labels = kmms.fit_predict(PCAembedding_low)

    # convert all clusters into two clusters
    kmerXTableInput["pLabels"] = kmms_labels
    kmerXTableInput["aLabels"] = y_hat
    newLabels_clusters_1 = kmerXTableInput[kmerXTableInput["aLabels"] == 1]["pLabels"].tolist()
    newLabels_clusters_0 = kmerXTableInput[kmerXTableInput["aLabels"] == 0]["pLabels"].tolist()
    newLabels = []

    for label in kmms_labels:
        if (label in newLabels_clusters_1) & (label not in newLabels_clusters_0):
            newLabels.append(1)
        else:
            newLabels.append(0)
            
    return newLabels, kmerXTableInput.drop(columns=["pLabels", "aLabels"])
    
def PCA2d(kTable, y_hat, filename):
    embedding = PCA(n_components= 2)
    nkTable = preprocessing.normalize(kTable)
    embedding.fit(nkTable)
    show = pd.DataFrame(embedding.transform(nkTable))
    # show kmeans clustering
    show["labels"] = y_hat
    ax = show[show["labels"]==1].plot.scatter(x=0, y=1, style="o", color="red", s=2)
    show[show["labels"]==0].plot.scatter(x=0, y=1, style="o", color="blue", s=2, ax=ax)
    show[(show["labels"]!= 0) & (show["labels"]!= 1) ].plot.scatter(x=0, y=1, style="o", color="gray", s=2, ax=ax)
    red = mpatches.Patch(color='red', label='Label 1')
    blue = mpatches.Patch(color='blue', label='Label 0')
    grey = mpatches.Patch(color='grey', label='Label -1')
    plt.legend(handles=[red, blue, grey])
    plt.title('PCA Visualization')
    plt.xlabel('First Principal Component')
    plt.ylabel('Second Principal Component')
    plt.savefig('nonNotebookFiles/' + filename + '.png', bbox_inches='tight')
    plt.close()
            
def tSNE2d(kTable, y_hat, filename):
    tSNEembedding = TSNE(n_components= 2, random_state = 0)
    tSNEembedding_low = tSNEembedding.fit_transform(kTable)
    show = pd.DataFrame(tSNEembedding_low)
    # show kmeans clustering
    show["labels"] = y_hat
    ax = show[show["labels"]==1].plot.scatter(x=0, y=1, style="o", color="red", s=2)
    show[show["labels"]==0].plot.scatter(x=0, y=1, style="o", color="blue", s=2, ax=ax)
    show[(show["labels"]!= 0) & (show["labels"]!= 1) ].plot.scatter(x=0, y=1, style="o", color="gray", s=2, ax=ax)
    red = mpatches.Patch(color='red', label='Label 1')
    blue = mpatches.Patch(color='blue', label='Label 0')
    grey = mpatches.Patch(color='grey', label='Label -1')
    plt.legend(handles=[red, blue, grey])
    plt.title('tSNE Visualization\n')
    plt.xlabel('First Component')
    plt.ylabel('Second Component')
    plt.savefig('nonNotebookFiles/' + filename + '.png', bbox_inches='tight')
    plt.close()