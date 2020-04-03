import pandas as pd
from Bio import SeqIO
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# the websiteScriptKmeans() method creates a png image and returns a dataframe with names, sequences, and labels


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

# Replace path1 and path01 with your file paths
def websiteScriptKmeans():
    # replace these paths for your computer
    path1 = "../data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta"
    path01 = "../data/mycovirus_genbank_all_refseq_nucleotide_unique.fasta"

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
    km4 = KMeans(random_state = 42, n_clusters = 2)
    km4.fit(kmer7Table01[cols])
    # generate array of labels
    y_hat4 = km4.predict(kmer7Table01[cols])

    #PCA
    embedding2 = PCA()
    embedding2.fit(kmer7Table01[cols])
    show2 = pd.DataFrame(embedding2.transform(kmer7Table01[cols]))
    # show kmeans clustering
    show2.plot.scatter(x=0, y=1, style="o", c=y_hat4, cmap = "viridis", s=2)
    plt.title('PCA Visualization for KMeans Clustering')
    plt.xlabel('First Principal Component')
    plt.ylabel('Second Principal Component')
    plt.savefig('visual.png', bbox_inches='tight')

    virus01["Labels"] = y_hat4
    return virus01

# Make a script to call websiteScriptKmeans(). This will save a png image and return a dataframe with names, sequences, and labels.

# To test this script uncomment the line below and run this script:
# print(websiteScriptKmeans().head())