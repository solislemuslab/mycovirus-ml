#import argparse
import gensim
import numpy as np
import tensorflow as tf
from tensorflow.keras.layers import Input, Layer, Dense, Activation, Flatten, Reshape, Dropout, BatchNormalization, Conv2D, Conv2DTranspose, UpSampling2D, LeakyReLU, LSTM, RepeatVector, TimeDistributed, Bidirectional
from tensorflow.keras import models
from tensorflow.keras.optimizers import Adam

class AE(object):
    def __init__(self, model_code, sub_sequence_len=32):
        self.sub_sequence_len = sub_sequence_len
        self.model_code = model_code

        self.seq_vectors = None
        self.ae = None
        self.encoder = None
        self.decoder = None
        self.x = None
        #self.y = None
        self.miss_count = 0

    def train(self, epochs, bs=32):
        # model: 0 for CNN-AE, 1 for LSTM-AE
        if self.model_code == 0:
            self.ae.fit(self.x[:,:,:,np.newaxis],self.x[:,:,:,np.newaxis],batch_size=bs,epochs=epochs, shuffle=True, verbose=1)
        elif self.model_code == 1:
            self.ae.fit(self.x,self.x,batch_size=bs,epochs=epochs,shuffle=True, verbose=1)
        else:
            raise ValueError('invalid model code')

    def save_model(self, savename):
        self.ae.reset_metrics()
        self.ae.save(savename)

    def compile_model(self):
        # model: 0 for CNN-AE, 1 for LSTM-AE
        if self.model_code == 0:
            self.compile_CNN((32, self.sub_sequence_len, 1))
        elif self.model_code == 1:
            self.compile_LSTM((32, self.sub_sequence_len))
        else:
            raise ValueError('invalid model code')

    def compile_CNN(self, input_shape):
        # Total params: 195,585
        def conv2d(x, filters, shape=(3, 3), **kwargs) :
            x = Conv2D(filters, shape, strides=(2, 2),padding='same')(x)
            x = BatchNormalization()(x)
            x = LeakyReLU()(x)
            x = Dropout(0.3)(x)
            return x
        i = Input(shape=input_shape)
        x = i
        x = conv2d(x, 16)
        x = conv2d(x, 32)
        x = conv2d(x, 64)
        x = conv2d(x, 128)
        x = Flatten()(x)
        latent_shape = x.shape[1:]
        self.encoder = models.Model(inputs=i, outputs=x)
        self.encoder.compile(optimizer='adam', loss='mean_squared_error')

        def deconv2d(x, filters, shape=(3, 3)):
            x = Conv2DTranspose(filters, shape, padding='same',strides=(2, 2))(x)
            x = BatchNormalization()(x)
            x = LeakyReLU()(x)
            x = Dropout(0.3)(x)
            return x
        i = Input(shape=latent_shape)
        x = i
        x = Reshape((2, 2, 128))(x)
        x = deconv2d(x, 64)
        x = deconv2d(x, 32)
        x = deconv2d(x, 16)
        x = Conv2DTranspose(1, (3, 3), padding='same', activation=None,strides=(2, 2))(x)
        self.decoder = models.Model(inputs=i, outputs=x)
        self.decoder.compile(optimizer='adam', loss='mean_squared_error')

        x = Input(shape=input_shape)
        latent = self.encoder(x)
        x_h = self.decoder(latent)
        self.ae = models.Model(inputs=x, outputs=x_h)
        self.ae.compile(optimizer='adam', loss='mean_squared_error')
        self.ae.summary()

    def compile_LSTM(self, input_shape):
        # Total params: 366,112
        i = Input(shape=input_shape)
        x = i
        x = Bidirectional(LSTM(128, activation='relu'))(x)
        print(x.shape)
        latent_shape = x.shape[1:]
        self.encoder = models.Model(inputs=i, outputs=x)
        self.encoder.compile(optimizer='adam', loss='mean_squared_error')
        
        i = Input(shape=latent_shape)
        x = i
        x = RepeatVector(32)(x)
        print(x.shape)
        x = LSTM(128, activation='relu', return_sequences=True)(x)
        print(x.shape)
        x = TimeDistributed(Dense(32, activation=None))(x)
        print(x.shape)
        self.decoder = models.Model(inputs=i, outputs=x)
        self.decoder.compile(optimizer='adam', loss='mean_squared_error')

        x = Input(shape=input_shape)
        latent = self.encoder(x)
        x_h = self.decoder(latent)
        self.ae = models.Model(inputs=x, outputs=x_h)
        self.ae.compile(optimizer='adam', loss='mean_squared_error')
        self.ae.summary()

    def load_embedding_data(self, Word2Vec_file):
        model = gensim.models.Word2Vec.load(Word2Vec_file)
        self.seq_vectors = model.wv
        del model
        # Possible use of seq_vectors:
        # get_vector(word)
        # get_keras_embedding(train_embeddings=False)
        # distances(word_or_vector, other_words=())
        # distance(w1, w2)
        print('loading sequences and embeddings...')
        embed_0 = self.get_training_data(fasta='../data/mycovirus_genbank_all_refseq_nucleotide_unique.fasta', label=0)
        embed_1 = self.get_training_data(fasta='../data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta', label=1)
        self.x = np.concatenate((np.array(self.embed_full(embed_0)), np.array(self.embed_full(embed_1))), axis=0)
        print('Input shape is: ')
        print(self.x.shape)

    def embed_full(self, seqs, dim=32):
        def embed(seq, dim=32):
            embeddings = np.zeros((dim, len(seq)))
            for i, s in enumerate(seq):
                try:
                    embeddings[:, i] = self.seq_vectors.get_vector(s)
                except KeyError:
                    embeddings[:, i] = 0.2*np.random.randn((32)) # random noise if nonexist
                    self.miss_count = self.miss_count + 1 # see how many not in vocabulary
            return embeddings

        embeddings = np.zeros((len(seqs), dim, len(seqs[0])))
        for i, s in enumerate(seqs):
            embeddings[i, :, :] = embed(s, dim)
        print('There are ' + str(self.miss_count) + ' sequences not in word2vec vocabulary')
        return embeddings

    def get_training_data(self, fasta, label, length=32):
        kmer_vectors = self.get_kmers(fasta=fasta)
        training = []
        for s in kmer_vectors:
            i = 0
            while i+length < len(s):
                training.append(s[i:i+length])
                i = i + length

        return training

    def get_kmers(self, fasta, k=6, hop=3):
        sequences_dict = self.parse_fasta(file=fasta)

        sequences = [sequences_dict[a] for a in sequences_dict]

        all_seq = []
        for s in sequences:
            temp = []
            i = 0
            while i+k < len(s):
                temp.append(s[i:i+k])
                i = i + hop
            all_seq.append(temp)
        return all_seq

    def parse_fasta(self, file):
        genomes = {}
        with open(file, "r") as f:
            for line in f:
                line = line.replace('\n', '')
                if line.startswith(">"):
                    curr = line
                    genomes[curr] = ''
                    continue
                genomes[curr] = genomes[curr] + line
        return genomes

if __name__ == "__main__":
    model_code = int(input('model code (0 for CNN, 1 for RNN)'))

    savename = input('model save name:')
    nepochs = int(input('how many epochs?'))
    
    ae = AE(model_code=model_code)
    ae.load_embedding_data(Word2Vec_file='virus2vec.model')
    ae.compile_model()
    ae.train(bs=32, epochs=nepochs)
    ae.save_model(savename=savename)


