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

My email to Helena:

I have several questions. 
If I understand correctly, in "Sclerotinia_biocontrol_mycovirus_nucleotide.fasta" (and "Sclerotinia_biocontrol_mycovirus_protein.fasta") you have the whole genome alignments (or protein alignments) for 7 (or 12) taxa that have biocontrol potential. 
While in "mycovirus_genbank_all_refseq_nucleotide_unique.fasta" ("mycovirus_genbank_all_refseq_protein_unique_proteindb.fasta") has 350 (or 465) taxa with whole genome alignments (or protein alignments) for viruses that we do not know if they are used as biocontrol or not.

Ideally, we need one fasta file with all (biocontrol, non-biocontrol), but I agree with what you say that the alignment will be meaningless as these viruses are very different.
Would it be the case also for the 110 viruses from the curated database (that they are difficult to be aligned)?
Because we could just use the closest to the sclerotinia viruses to have a semi-meaningful alignment. I would prefer not to just use one conserved gene because the algorithm will find for patterns in whatever alignment we give. If we only give this one gene, we are basically saying that the signal for biocontrol is in this gene only. 

So, the decision is to choose the right viruses so that we can have one fasta file with biocontrol and non-biocontrol viruses with whole genomes alignments (or protein alignments) so that alignment is not totally garbage, but also we have options outside sclerotinia. I know that this is difficult. If you think this is impossible, then we can meet and talk about alternatives.