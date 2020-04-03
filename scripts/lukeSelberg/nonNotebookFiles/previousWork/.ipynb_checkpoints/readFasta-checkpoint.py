from Bio import SeqIO
import pandas as pd

#function to read in fasta files
def readFasta(file):
    # read in sequence and id separately
    fasta_sequences = SeqIO.parse(file,'fasta')
    df_1 = pd.DataFrame(fasta_sequences)

    df_1["ID"] = [fasta.id for fasta in SeqIO.parse(file, "fasta")]

    # place id column at front of dataframe
    cols = list(df_1.columns)
    cols = [cols[-1]] + cols[:-1]
    df_1 = df_1[cols]
    df_1.set_index('ID', inplace = True)
    
    return df_1