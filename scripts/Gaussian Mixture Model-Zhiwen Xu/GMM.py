# import packages
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pdz
import numpy as np
import pandas as pd

from pandas import Series, DataFrame
import Bio
from Bio import SeqIO,AlignIO

import sklearn
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.mixture import GaussianMixture as GMM

from sklearn.manifold import TSNE
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import time


# methods

def parseFasta(data):
    d = {fasta.id: str(fasta.seq) for fasta in SeqIO.parse(data, "fasta")}
    pd.DataFrame([d])
    s = pd.Series(d, name='Sequence')
    s.index.name = 'ID'
    s.reset_index()
    return pd.DataFrame(s)


def get_kmer_table(path, k_min, k_max):
    genes, gene_len = read_fasta(path)
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


def read_fasta(paths):
    all_genes = []
    all_gene_len = []
    
    for path in paths:
        virus = parseFasta(path)
        virus = virus.drop_duplicates(keep="last")
        genes = list(virus['Sequence'])
        genes_seq = get_gene_sequences(path)
        gene_len = get_gene_len(genes_seq)
        all_genes = all_genes + genes_seq
        all_gene_len = all_gene_len + gene_len
    return all_genes,all_gene_len


def get_predictions(path, k_min, k_max, num_class, cov_type, seed):
    kmer_table = get_kmer_table(path, k_min, k_max)
    gmm = GMM(n_components=num_class, covariance_type=cov_type, random_state=seed).fit(kmer_table)
    predictions = gmm.predict(kmer_table)
    return predictions


def get_predictions_semi(path,k_min,k_max,num_class,cov_type,seed,labels):
    targets = []
    kmer_table = get_kmer_table(path, k_min, k_max)
    finalDf = pd.concat([kmer_table, pd.Series(labels)], axis = 1)
    gmm = GMM(n_components=num_class,covariance_type=cov_type,random_state = seed)
    for i in range(num_class):
        if (i in list(finalDf.Labels)):
            targets.append(i)
    if (len(targets)==num_class):
        gmm.means_init = np.array([kmer_table[finalDf.Labels == i].mean(axis=0) for i in targets])
    gmm.fit(kmer_table)
    predictions = gmm.predict(kmer_table)
    return predictions


def cal_accuracy(labels, predictions):
    err = 0
    total_len = len(labels)
    for i in range(len(labels)):
        if (labels[i] == -1):
            total_len = total_len - 1
            continue
        if (labels[i] != predictions[i]):
            err += 1

    return 1 - err / (total_len)


def sammon(x, n, display=2, maxhalves=20, maxiter=500, tolfun=1e-9):
    import numpy as np
    from scipy.spatial.distance import cdist

    D = cdist(x, x)  # distance matrix
    init = 'pca'  # initialization from pca

    if np.count_nonzero(np.diagonal(D)) > 0:
        raise ValueError("The diagonal of the dissimilarity matrix must be zero")

    # Remaining initialisation
    N = x.shape[0]
    scale = 0.5 / D.sum()
    D = D + np.eye(N)

    if np.count_nonzero(D <= 0) > 0:
        raise ValueError("Off-diagonal dissimilarities must be strictly positive")

    Dinv = 1 / D
    [UU, DD, _] = np.linalg.svd(x)
    y = UU[:, :n] * DD[:n]

    one = np.ones([N, n])
    d = cdist(y, y) + np.eye(N)
    dinv = 1. / d
    delta = D - d
    E = ((delta ** 2) * Dinv).sum()

    for i in range(maxiter):

        # Compute gradient, Hessian and search direction (note it is actually
        # 1/4 of the gradient and Hessian, but the step size is just the ratio
        # of the gradient and the diagonal of the Hessian so it doesn't
        # matter).
        delta = dinv - Dinv
        deltaone = np.dot(delta, one)
        g = np.dot(delta, y) - (y * deltaone)
        dinv3 = dinv ** 3
        y2 = y ** 2
        H = np.dot(dinv3, y2) - deltaone - np.dot(2, y) * np.dot(dinv3, y) + y2 * np.dot(dinv3, one)
        s = -g.flatten(order='F') / np.abs(H.flatten(order='F'))
        y_old = y

        # Use step-halving procedure to ensure progress is made
        for j in range(maxhalves):
            s_reshape = np.reshape(s, (-1, n), order='F')
            y = y_old + s_reshape
            d = cdist(y, y) + np.eye(N)
            dinv = 1 / d
            delta = D - d
            E_new = ((delta ** 2) * Dinv).sum()
            if E_new < E:
                break
            else:
                s = 0.5 * s

        # Bomb out if too many halving steps are required
        if j == maxhalves - 1:
            print('Warning: maxhalves exceeded. Sammon mapping may not converge...')

        # Evaluate termination criterion
        if abs((E - E_new) / E) < tolfun:
            if display:
                print('TolFun exceeded: Optimisation terminated')
            break

        # Report progress
        E = E_new
        if display > 1:
            print('epoch = %d : E = %12.10f' % (i + 1, E * scale))

    if i == maxiter - 1:
        print('Warning: maxiter exceeded. Sammon mapping may not have converged...')

    # Fiddle stress to match the original Sammon paper
    E = E * scale

    return [y, E]


def sammon_plot(y, plot_labels, path):
    x_axis = []
    y_axis = []

    for i in range(len(y)):
        x_axis.append(y[i][0])
        y_axis.append(y[i][1])

    sns.scatterplot(x_axis, y_axis, hue=plot_labels, legend='full')
    # images.append(path)
    plt.savefig(path)
    plt.close()
    # plt.show()


def cal_accuracy(labels, predictions):
    err = 0
    total_len = len(labels)
    for i in range(len(labels)):
        if (labels[i] == -1):
            total_len = total_len - 1
            continue
        if (labels[i] != predictions[i]):
            err += 1

    return 1 - err / (total_len)


def PCA_plot(x, y, n_dim, path):
    x = StandardScaler().fit_transform(x)
    pca = PCA(n_components=n_dim)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data=principalComponents, columns=['principal component 1', 'principal component 2'])
    finalDf = pd.concat([principalDf, pd.Series(y)], axis=1)
    finalDf.columns = ['principal component 1', 'principal component 2', 'target']

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel('Principal Component 1', fontsize=15)
    ax.set_ylabel('Principal Component 2', fontsize=15)
    ax.set_title('2 component PCA', fontsize=20)
    targets = [0, 1]
    colors = ['r', 'g']
    for target, color in zip(targets, colors):
        indicesToKeep = finalDf['target'] == target
        ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
                   , finalDf.loc[indicesToKeep, 'principal component 2']
                   , c=color)
    ax.legend(targets)
    # images.append(path)
    fig.savefig(path)
    plt.close(fig)
    # plt.show()

def tsne_plot(x,y,path,title):
    tsne = TSNE()
    X_embedded = tsne.fit_transform(x)
    sns.scatterplot(X_embedded[:,0], X_embedded[:,1], hue=y, legend='full')
    #images.append(path)
    plt.xlabel('Dimension 1')
    plt.ylabel('Dimension 2')
    plt.title(title)
    plt.savefig(path)
    plt.close()

def model_selection(path,labels,num_class):
    best_accu = 0
    best_prediction = []
    cov_type = ['full','diag','tied','spherical']
    k_min = [2,3,4]
    k_max = [2,3,4,5]
    for cov in cov_type:
        for k1 in k_min:
            for k2 in k_max:
                if (k2 >= k1):
                    prediction = get_predictions_semi(path,k1,k2,num_class,cov,0,labels)
                    accu = cal_accuracy(labels,prediction)
                    if accu > best_accu:
                        best_accu = accu
                        best_kmin = k1
                        best_kmax = k2
                        best_cov = cov
                        best_prediction = prediction
    print('Best model has the following parameters:')
    print('minimum length of kmer: ', best_kmin)
    print('maximum length of kmer: ', best_kmax)
    print('covariance type: ', best_cov)
    print('It has an accuracy regard to known labels of ',best_accu)
    return best_kmin,best_kmax,best_cov,best_prediction

