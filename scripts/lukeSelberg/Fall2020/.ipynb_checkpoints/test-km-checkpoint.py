import unittest

from websiteScripts import *

def kmeans_semiSupervised(fasta, y_hat, klength_min = 6, klength_max = 6, rNum = 50):
    inputData = parseFasta(fasta)
    inputData["Sequence"] = inputData["Sequence"].apply(lambda x: x.replace("-", ""))
    kmerXTableInput = kmerXTable(inputData, klength_min, klength_max)
    
    PCAembedding = PCA(n_components=10)
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
        if (newLabels_clusters_1.count(label) > newLabels_clusters_0.count(label)):
            newLabels.append(1)
        else:
            newLabels.append(0)
            
    return newLabels, kmerXTableInput.drop(columns=["pLabels", "aLabels"])

PATH1 = "../data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta"
PATH01 = "../data/mycovirus_genbank_all_refseq_nucleotide_unique.fasta"

labels_true = [-1 for x in range(350)]

y_hat_us_true = kmeans(PATH01, 2, 6, 7, 50)
y_hat_s_true = kmeans_semiSupervised(PATH01, labels_true, 6, 7, 50)

class TestKM(unittest.TestCase):
    def test_km_us(self):
        data = PATH01
        result = kmeans(PATH01, 2, 6, 7, 50)
        for i in range(0, len(result)):
            self.assertEqual(result[0][i], y_hat_us_true[0][i])
        
    def test_km_s(self):
        data = PATH01
        result = kmeans_semiSupervised(PATH01, labels_true, 6, 7, 50)
        for i in range(0, len(result)):
            self.assertEqual(result[0][i], y_hat_s_true[0][i])
        
        
if __name__ == '__main__':
    unittest.main()