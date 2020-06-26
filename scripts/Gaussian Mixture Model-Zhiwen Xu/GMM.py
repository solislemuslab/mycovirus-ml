
# import packages
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import numpy as np

from pandas import Series, DataFrame
import Bio
from Bio import SeqIO,AlignIO
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.mixture import GaussianMixture as GMM

# methods

# parseFasta(data) credit to Luke
def parseFasta(data):
    d = {fasta.id : str(fasta.seq) for fasta in SeqIO.parse(data, "fasta")}
    pd.DataFrame([d])
    s = pd.Series(d, name='Sequence')
    s.index.name = 'ID'
    s.reset_index()
    return pd.DataFrame(s)

def get_kmer_table(path1,path2,k_min,k_max):
    genes,gene_len = read_fasta(path1,path2)
    count_vect = CountVectorizer(analyzer='char', ngram_range=(k_min, k_max))
    X = count_vect.fit_transform(genes)
    chars = count_vect.get_feature_names()
    kmers = X.toarray()
    kmer_freq = []
    for i in range(len(genes)):
        kmer_freq.append(kmers[i] / gene_len[i])
    input = pd.DataFrame(kmer_freq, columns=chars)
    return input

def get_gene_sequences(filename):
    genes = []
    for record in SeqIO.parse(filename, "fasta"):
        genes.append(str(record.seq))
    return genes

# genes: a list of gene sequences, which can directly be generated from get_gene_sequences().
def get_gene_len(genes):
    gene_len = []

    for i in range(len(genes)):
        gene_len.append(len(genes[i]))
    return gene_len

def read_fasta(path1,path2):
    virus1 = parseFasta(path1)
    # put confirmed virus killers at bottom, and removed the duplicates already in the data
    virus01 = parseFasta(path1)
    virus01 = virus01.append(virus1)
    virus01 = virus01.drop_duplicates(keep="last")
    genes = list(virus01['Sequence'])
    genes_0 = get_gene_sequences(path1)
    genes_1 = get_gene_sequences(path2)
    gene_len_0 = get_gene_len(genes_0)
    gene_len_1 = get_gene_len(genes_1)
    all_genes = genes
    all_gene_len = gene_len_0 + gene_len_1
    return all_genes,all_gene_len

# change the following parameters to user inputs
path1 = "label0.fasta"
path2 = "label1.fasta"
k_min = 2
k_max = 3
num_type = 2

kmer_table = get_kmer_table(path1,path2,k_min,k_max)
gmm = GMM(n_components=num_type).fit(kmer_table)
labels = gmm.predict(kmer_table)

# predicted results are stored in labels