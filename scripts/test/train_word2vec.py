import numpy as np
#import argparse
import gensim

class Virus2Vec(object):
    def __init__(self, k, hop, filename, dim=32, train_window=10):
        self.k = k
        self.hop = hop
        self.filename = filename
        self.dim = dim
        self.train_window = train_window
        self.model = None
        self.load_sequences()

    def load_sequences(self):
        genomes = {}
        with open(self.filename, "r") as f:
            for line in f:
                line = line.replace('\n', '')
                if line.startswith(">"):
                    curr = line
                    genomes[curr] = ''
                    continue
                genomes[curr] = genomes[curr] + line
        self.all_sequences = [genomes[a] for a in genomes]

    def train_model(self):
        all_seq_5 = []
        for s in self.all_sequences:
            temp = []
            i = 0
            while i+self.k < len(s):
                temp.append(s[i:i+self.k])
                i = i + self.hop
            all_seq_5.append(temp)
        self.model = gensim.models.Word2Vec(all_seq_5, min_count=1,
            size=self.dim, workers=5, window=self.train_window, 
            iter=300, sg=1)

if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument('train', metavar='T', type=int,
 #                    help='0 for load 1 for train')
    # args = parser.parse_args()
    # if args.train:
    #   embed = Virus2Vec(k=6, filename='sequences.fasta')
    #   embed.save("virus2vec.model")
    # else 
    #   embed = Word2Vec.load("virus2vec.model")

    virus2vec = Virus2Vec(k=6, hop=3, filename='sequences.fasta')
    virus2vec.train_model()
    virus2vec.model.save("virus2vec.model")