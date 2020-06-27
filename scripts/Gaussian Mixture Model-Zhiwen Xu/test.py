import GMM

# change the following parameters to user inputs
path1 = "label0.fasta"
path2 = "label1.fasta"
k_min = 2
k_max = 3
num_class = 2
cov_type = 'full'
predictions = GMM.get_predictions(path1,path2,k_min,k_max,num_class,cov_type)