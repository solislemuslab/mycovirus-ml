# Data

In UW Box from Helena Jaramillo Mesa.

```shell
(master) $ pwd
/Users/Clauberry/Dropbox/Documents/solislemus-lab/student-projects/mycovirus-ml/data
(master) $ ls
Sclerotinia_biocontrol_mycovirus_nucleotide.fasta
Sclerotinia_biocontrol_mycovirus_protein.fasta
mycovirus_genbank_all_refseq_nucleotide_unique.fasta
mycovirus_genbank_all_refseq_protein_unique_proteindb.fasta
mycovirus_list_viralzone.md
```

Number of taxa:
```shell
(master) $ grep ">" Sclerotinia_biocontrol_mycovirus_protein.fasta | wc -l
      12
(master) $ grep ">" Sclerotinia_biocontrol_mycovirus_nucleotide.fasta | wc -l
       7
(master) $ grep ">" mycovirus_genbank_all_refseq_nucleotide_unique.fasta | wc -l 
     350
(master) $ grep ">" mycovirus_genbank_all_refseq_protein_unique_proteindb.fasta | wc -l
     465
```

- My email to Helena (11/25):

I have several questions. 
If I understand correctly, in "Sclerotinia_biocontrol_mycovirus_nucleotide.fasta" (and "Sclerotinia_biocontrol_mycovirus_protein.fasta") you have the whole genome alignments (or protein alignments) for 7 (or 12) taxa that have biocontrol potential. 
While in "mycovirus_genbank_all_refseq_nucleotide_unique.fasta" ("mycovirus_genbank_all_refseq_protein_unique_proteindb.fasta") has 350 (or 465) taxa with whole genome alignments (or protein alignments) for viruses that we do not know if they are used as biocontrol or not.

Ideally, we need one fasta file with all (biocontrol, non-biocontrol), but I agree with what you say that the alignment will be meaningless as these viruses are very different.
Would it be the case also for the 110 viruses from the curated database (that they are difficult to be aligned)?
Because we could just use the closest to the sclerotinia viruses to have a semi-meaningful alignment. I would prefer not to just use one conserved gene because the algorithm will find for patterns in whatever alignment we give. If we only give this one gene, we are basically saying that the signal for biocontrol is in this gene only. 

So, the decision is to choose the right viruses so that we can have one fasta file with biocontrol and non-biocontrol viruses with whole genomes alignments (or protein alignments) so that alignment is not totally garbage, but also we have options outside sclerotinia. I know that this is difficult. If you think this is impossible, then we can meet and talk about alternatives.

- Helena's response (11/25):
The viruses that were reported as biocontrols for Sclerotinia were 7, 1 of which has no available sequences. 5 of these are non-segmented viruses (1 fasta sequence for the whole genome), and 1 has 2 segments (so 2 fasta sequences as genome). The file "Sclerotinia_biocontrol_mycovirus_nucleotide.fasta" has therefore, 7 sequences in total. They are not aligned, because I was waiting to see what you thought of what I asked you about the alignment. The file "Sclerotinia_biocontrol_mycovirus_protein.fasta" has the same 6 virus, but now it contains the protein sequences as they can be found in genbank (as a polyprotein for some of them, or as individual proteins for the others). These are also, not aligned. 

The files "mycovirus_genbank_all_refseq_nucleotide_unique.fasta" and mycovirus_genbank_all_refseq_protein_unique_proteindb.fasta" contain all the mycoviruses in genbank (including the Sclerotinia controls that I copied into the other files). None of the files is yet aligned.

The ~110 from ViralZone will indeed pose the same problem, but I don't think it is impossible, we just need to find an alignment that will work for what you want. An alignment can be made with all of them, I just think it will be too much noise to be useful. If you want, we can meet so that I can show you what I am talking about when I say that the viruses are too different, and we can decide what to use for the alignment. 


## After meeting 12/3

We will focus on the two files:
- `mycovirus_genbank_all_refseq_nucleotide_unique.fasta` (which have label 0)
- `Sclerotinia_biocontrol_mycovirus_nucleotide.fasta` (which have label 1)

These files have 350 and 7 viral sequences:
```shell
(master) $ grep ">" Sclerotinia_biocontrol_mycovirus_nucleotide.fasta | wc -l
       7
(master) $ grep ">" mycovirus_genbank_all_refseq_nucleotide_unique.fasta | wc -l 
     350
```

## Toy dataset

- Data input (genomes) of Staph bacteria: `core_gene_alignment-narsa.aln` 
    - rows=bacterial strains (individuals)
    - columns=nucleotide sites (features)
- Labels for Staph bacteria (0=susceptible to antibiotic, 1=resistant to antibiotic): `responses-staph.csv`

## New aligned data

Helena sent us the `Alignments.zip` file (1/16), and this is her email:

I just uploaded to the shared Box a folder named "Alignments" that contains the files with aligned sequences for the Hypoviridae family. This family contains most of the Sclerotinia mycoviruses and is the family with the largest number of viruses associated with hypovirulence in fungi. It has right now 19 members, of which only 16 has complete polyprotein sequences available. Those 16 sequences are the ones I used for the alignments. 

Of those 16 sequences, 7 have been shown to cause hypovirulence in fungi. In the file "Hypoviridae.csv" you will find the information for all of this including the name, the accession number, the abbreviation, if they have biocontrol characteristics and if there is a complete sequence. For the biocontrol column I added a Y for Yes (for the 7 that cause hypovirulence) or a N for No (the 9 that haven't been studied as hypovirulent). The last column is the Complete Sequence. If they have a C is because they are complete including the UTR regions, or if they have a P is because they only have the complete polyprotein (No UTRs). 

For the alignments I used the polyprotein. You will find a file for all of them (the 16 I used) aligned both in nucleotide ("hypoviridae_aligned_polyprotein_nucleotide_all.fas") or in amino acids ("hypoviridae_aligned_polyprotein_aminoacid_all.fas"). I also added the aligned filed that contains just the biocontrol ones (7), in both nucleotides ("hypoviridae_aligned_polyprotein_nucleotide_biocontrol.fas") and amino acids ("hypoviridae_aligned_polyprotein_aminoacid_biocontrol.fas"), and the last 2 files are the aligned sequences but for the non-biocontrol ones (9) in both nucleotides ("hypoviridae_aligned_polyprotein_nucleotide_nobiocontrol.fas") and amino acids ("hypoviridae_aligned_polyprotein_aminoacid_nobiocontrol.fas"). 

I also think that if the number of sequences is not sufficient, we could add other families (Endornaviridae and Narnaviridae) that also have 1 segment of ssRNA as genome and are non encapsidated. They might differ a lot, but at least they only have 1 RNA segment, which we could align and see if they have genetic similarities of some sort.

## Description of the data (meeting 1/31)
- We have two main groups of files (unaligned/aligned) and within each group, there are two subgroups (protein coded or nucleotide coded)
     - Unaligned
          - Protein coded: two files 1) virus that kill fungus, 2) virus that do not kill fungus
          - Nucleotide coded: two files 1) virus that kill fungus, 2) virus that do not kill fungus
     - Aligned
          - Protein coded: two files 1) virus that kill fungus, 2) virus that do not kill fungus
          - Nucleotide coded: two files 1) virus that kill fungus, 2) virus that do not kill fungus

Next steps:
- More complete picture of data based on dimension (how many viruses in each file, how many columns), missingness, how many kill fungus vs how many do not kill fungus, differences between protein and nucleotide coded files, how many columns are identical across virus?
- Conglomerate data: we want to file 1 genetic file for aligned (with fungus killers and non-killers in the same file) and 1 genetic file for unaligned (with fungus killers and non-killers in the same file). In addition, we want to have 1 "labels file" for aligned that simply list a 0/1 vector (0=non-killer, 1=killer) for each virus in the genetic file. Same for unaligned (1 "labels file")
- Think about ways in which could do the alignment automatically (as step 1 in the machine-learning process): the idea is to minimize discrepancies across columns. Maybe this is too hard to do (impossible?), so we should also explore data augmentation techniques

# Analyses

Next steps: We want to fit statistical/machine-learning models to accomplish two tasks:
- feature selection: identify genes associated with mycovirus fungus-killing capability
- prediction: for a given new genome, can we predict whether it will be fungus killer or not

Methods:
- regression (we need to explore penalized regression because we have more features than individuals)
- random forest
- neural networks
- ...

## Main difficulties of the project
- We do not have aligned sequences (can we align ourselves?)
- Input data is letters/categories: ACGT, cannot be treated as numbers 1234


## Previous work for toy dataset

Claudia had fit naive neural networks and random forest in Julia. All scripts in `scripts/previous-work`:
- `notebook.md`: explains pre-processing of the data
- `*.jl`: julia scripts, described in `notebook.md`


Next steps:
- wait for alignments
- use existing data files as "alignments" to start machine-learning and regression pipelines
- play with toy example