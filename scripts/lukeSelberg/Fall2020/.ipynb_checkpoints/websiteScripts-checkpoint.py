import pandas as pd
from Bio import SeqIO
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
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
    
def kmeans(fasta, klength, rNum, cNum):
    inputData = parseFasta(fasta)
    temp = virus01.append(inputData)
    temp = temp.drop_duplicates(keep="last")
        
    temp["Sequence"] = temp["Sequence"].apply(lambda x: x.replace("-", ""))
    kmerXTableInput = kmerXTable(temp, klength, klength)
        
        
    km = KMeans(random_state = rNum, n_clusters = cNum)
    km.fit(kmerXTableInput) 
    y_hat = km.predict(kmerXTableInput)
        
    return y_hat, kmerXTableInput
    
def PCA2d(kTable, y_hat, filename):
    embedding = PCA()
    embedding.fit(kTable)
    show = pd.DataFrame(embedding.transform(kTable))
    # show kmeans clustering
    show["labels"] = y_hat
    ax = show[show["labels"]==1].plot.scatter(x=0, y=1, style="o", color="red", s=2)
    show[show["labels"]==0].plot.scatter(x=0, y=1, style="o", color="blue", s=2, ax=ax)
    red = mpatches.Patch(color='red', label='Fungus Killers')
    blue = mpatches.Patch(color='blue', label='Fungus Non-Killers')
    plt.legend(handles=[red, blue])
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
    red = mpatches.Patch(color='red', label='Fungus Killers')
    blue = mpatches.Patch(color='blue', label='Fungus Non-Killers')
    plt.legend(handles=[red, blue])
    plt.title('tSNE Visualization\n')
    plt.xlabel('First Component')
    plt.ylabel('Second Component')
    plt.savefig('nonNotebookFiles/' + filename + '.png', bbox_inches='tight')
    plt.close()
    
# def kmeans_semiSupervised():
    
PATH1 = "../data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta"
PATH01 = "../data/mycovirus_genbank_all_refseq_nucleotide_unique.fasta"

virus1 = parseFasta(PATH1)
virus01 = parseFasta(PATH01)
virus01 = virus01.append(virus1)
virus01 = virus01.drop_duplicates(keep="last")