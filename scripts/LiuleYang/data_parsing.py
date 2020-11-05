from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import randint
cannot_kill_DNA_records = list(SeqIO.parse(open("dataset/mycovirus_genbank_all_refseq_nucleotide_unique(0).fasta"),'fasta'))
kill_DNA_records = list(SeqIO.parse("dataset/Sclerotinia_biocontrol_mycovirus_nucleotide(1).fasta",'fasta'))

print("Original counting for cannot: ", len(cannot_kill_DNA_records))
print("Original count for can: ", len(kill_DNA_records))


# counting:

# R	A or G
# Y	C or T
# S	G or C
# W	A or T
# K	G or T
# M	A or C
# B	C or G or T
# D	A or G or T
# H	A or C or T
# V	A or C or G
# N all

total_count = 0;

j = 0;
while(j <len(cannot_kill_DNA_records) ):
    i = 0;
    single_count = 1;
    while(i < len(cannot_kill_DNA_records[j])):
        if(cannot_kill_DNA_records[j][i] == 'R' or cannot_kill_DNA_records[j][i] == 'Y' or cannot_kill_DNA_records[j][i] == 'S' or cannot_kill_DNA_records[j][i] == 'W'
                or cannot_kill_DNA_records[j][i] == 'K' or cannot_kill_DNA_records[j][i] == 'M'):
            single_count = single_count * 2;
        elif cannot_kill_DNA_records[j][i] == 'B' or cannot_kill_DNA_records[j][i] == 'D' or cannot_kill_DNA_records[j][i] == 'H' or cannot_kill_DNA_records[j][i] == 'V':
            single_count = single_count * 3;
        elif cannot_kill_DNA_records[j][i] == 'N' :
            single_count = single_count * 4;
        i += 1;
    total_count += single_count;
    j += 1;
print("Cannot kill total count: ", total_count);

kill_totoal = 0;
q = 0
while(q <len(kill_DNA_records) ):
    p = 0;
    kill_single_count = 1;
    while(p < len(kill_DNA_records[q])):
        if(kill_DNA_records[q][p] == 'R' or kill_DNA_records[q][p] == 'Y' or kill_DNA_records[q][p] == 'S' or kill_DNA_records[q][p] == 'W'
                or kill_DNA_records[q][p] == 'K' or kill_DNA_records[q][p] == 'M'):
            kill_single_count = kill_single_count * 2;
        elif kill_DNA_records[q][p] == 'B' or kill_DNA_records[q][p] == 'D' or kill_DNA_records[q][p] == 'H' or kill_DNA_records[q][p] == 'V':
            kill_single_count = kill_single_count * 3;
        elif kill_DNA_records[q][p] == 'N':
            kill_single_count = kill_single_count * 4;
        p += 1;
    kill_totoal += kill_single_count;
    q += 1;
print("kill total count: ", kill_totoal);

# We have now 700 sequences in taotal

# Parsing:
# R	A or G
# Y	C or T
# S	G or C
# W	A or T
# K	G or T
# M	A or C
# B	C or G or T
# D	A or G or T
# H	A or C or T
# V	A or C or G
# N all
f= open("agumentation_cannot.txt","w+")
j = 0;
while(j <len(cannot_kill_DNA_records) ):
    i = 0;
    sequences = []
    while(i < len(cannot_kill_DNA_records[j])):
        if(cannot_kill_DNA_records[j][i] == 'R'):
            p = 0
            if(len(sequences) == 0):
                sequences.append(['A'])
                sequences.append(['G'])
            else:
                while(p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('A')
                    p += 1
                while (p < len(sequences)):
                    sequences[p].append('G')
                    p += 1
        elif(cannot_kill_DNA_records[j][i] == 'Y'):
            p = 0
            if (len(sequences) == 0):
                sequences.append(['C'])
                sequences.append(['T'])
            else:
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('C')
                    p += 1
                while (p < len(sequences)):
                    sequences[p].append('T')
                    p += 1
        elif (cannot_kill_DNA_records[j][i] == 'S'):
            p = 0
            if (len(sequences) == 0):
                sequences.append(['G'])
                sequences.append(['C'])
            else:
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('G')
                    p += 1
                while (p < len(sequences)):
                    sequences[p].append('C')
                    p += 1
        elif (cannot_kill_DNA_records[j][i] == 'W'):
            p = 0
            if (len(sequences) == 0):
                sequences.append(['A'])
                sequences.append(['T'])
            else:
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('A')
                    p += 1
                while (p < len(sequences)):
                    sequences[p].append('T')
                    p += 1
        elif (cannot_kill_DNA_records[j][i] == 'K'):
            p = 0
            if (len(sequences) == 0):
                sequences.append(['G'])
                sequences.append(['T'])
            else:
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('G')
                    p += 1
                while (p < len(sequences)):
                    sequences[p].append('T')
                    p += 1
        elif (cannot_kill_DNA_records[j][i] == 'M'):
            p = 0
            if (len(sequences) == 0):
                sequences.append(['A'])
                sequences.append(['C'])
            else:
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('A')
                    p += 1
                while (p < len(sequences)):
                    sequences[p].append('C')
                    p += 1
        elif (cannot_kill_DNA_records[j][i] == 'B'):
            p = 0
            if (len(sequences) == 0):
                sequences.append(['C'])
                sequences.append(['G'])
                sequences.append(['T'])
            else:
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('C')
                    p += 1
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('G')
                    p += 1
                while (p < len(sequences)):
                    sequences[p].append('T')
                    p += 1
        elif (cannot_kill_DNA_records[j][i] == 'D'):
            p = 0
            if (len(sequences) == 0):
                sequences.append(['A'])
                sequences.append(['G'])
                sequences.append(['T'])
            else:
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('A')
                    p += 1
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('G')
                    p += 1
                while (p < len(sequences)):
                    sequences[p].append('T')
                    p += 1
        elif (cannot_kill_DNA_records[j][i] == 'H'):
            p = 0
            if (len(sequences) == 0):
                sequences.append(['A'])
                sequences.append(['C'])
                sequences.append(['T'])
            else:
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('A')
                    p += 1
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('C')
                    p += 1
                while (p < len(sequences)):
                    sequences[p].append('T')
                    p += 1
        elif (cannot_kill_DNA_records[j][i] == 'V'):
            p = 0
            if (len(sequences) == 0):
                sequences.append(['A'])
                sequences.append(['C'])
                sequences.append(['G'])
            else:
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('A')
                    p += 1
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('C')
                    p += 1
                while (p < len(sequences)):
                    sequences[p].append('G')
                    p += 1
        elif (cannot_kill_DNA_records[j][i] == 'N'):
            p = 0
            if (len(sequences) == 0):
                sequences.append(['A'])
                sequences.append(['T'])
                sequences.append(['C'])
                sequences.append(['G'])
            else:
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('A')
                    p += 1
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('T')
                    p += 1
                while (p < len(sequences)):
                    sequences.append(sequences[p])
                    sequences[p].append('C')
                    p += 1
                while (p < len(sequences)):
                    sequences[p].append('G')
                    p += 1
        else:
            p = 0
            if (len(sequences) == 0):
                sequences.append([cannot_kill_DNA_records[j][i]])
            else:
                while (p < len(sequences)):
                    sequences[p].append(cannot_kill_DNA_records[j][i])
                    p += 1
        i += 1;
    print = 0
    while(print < len(sequences)):
        letter = 0
        while(letter < len(sequences[print])):
            f.write(sequences[print][letter])
            letter += 1
        f.write('\n')
        print +=1
    j += 1;



