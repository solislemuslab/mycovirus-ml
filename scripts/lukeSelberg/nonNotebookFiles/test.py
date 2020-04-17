import websiteScriptKmeans as wsk

A = ["../data/hypoviridae_aligned_polyprotein_nucleotide_biocontrol.fas", "../data/hypoviridae_aligned_polyprotein_nucleotide_nobiocontrol.fas", "../data/hypoviridae_aligned_polyprotein_nucleotide_all.fas"]

UA = ["../data/mycovirus_genbank_all_refseq_nucleotide_unique.fasta",  "../data/Sclerotinia_biocontrol_mycovirus_nucleotide.fasta"]

output = wsk.websiteScriptKmeans(A)

