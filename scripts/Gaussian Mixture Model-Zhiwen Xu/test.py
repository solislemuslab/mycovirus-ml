import GMM
import numpy as np

# change the following parameters to user inputs
path1 = "datasets/bat_flu.fa"
path2 = "datasets/penguin_flu.fa"
k_min = 2
k_max = 3
num_class = 2
cov_type = 'full'
predictions = GMM.get_predictions(path1,path2,k_min,k_max,num_class,cov_type)

df = GMM.get_kmer_table(["datasets/bat_flu.fa","datasets/penguin_flu.fa"],k_min,k_max)
df2 = df.drop_duplicates(keep="last")
kept_viruses = df2.index.values
bat_len = len(GMM.get_gene_sequences("datasets/bat_flu.fa"))
penguin_len = len(GMM.get_gene_sequences("datasets/penguin_flu.fa"))
zeros = [0]*bat_len
labels1 = np.append(zeros, [1]*penguin_len, axis=None)
plot_labels = labels1[kept_viruses]
X = df2.as_matrix(columns = df2.columns)
[y, E] = GMM.sammon(X,2)
GMM.sammon_plot(y,plot_labels)
GMM.PCA_plot(df,labels1,num_class)
GMM.tsne_plot(df,labels1)

GMM.model_selection(["datasets/bat_flu.fa", "datasets/penguin_flu.fa"], labels1, num_class)