# Julia script to concatenate fasta files
# copied from chlamydia/scripts, but modified
# Claudia August 2019

datafolder = "../data/pseudomonas/sequences/"
outfile = string(datafolder,"pseudomonas.fasta")
datafiles = string(datafolder,"CDS-OGs/")
files = String[]

for f in filter(x -> endswith(x, ".fasta"), readdir(datafiles))
    push!(files,f)
end

println("found $(length(files)) fasta files")

taxa = String[]
for f in files
    io1 = open(string(datafiles,f), "r")
    lines1 = readlines(io1);
    for l in lines1
        if occursin('>',l)
            s = split(split(l,'|')[1],'>')[2]
            push!(taxa,s)
        end
    end
end
taxa = unique(taxa)
n = length(taxa)
seqs = fill("",n)

for f in files
    @show f
    lines = readlines(string(datafiles,f));
    tmptaxa = String[]
    tmpseqs = fill("",n)
    i = [0]
    for l in lines
        if occursin('>',l)
            s = split(split(l,'|')[1],'>')[2]
            push!(tmptaxa,s)
            i[1] += 1
        else
            tmpseqs[i[1]] = string(tmpseqs[i[1]], l)
        end
    end
    tmpseqs = tmpseqs[1:length(tmptaxa)];
    nseq = unique(length.(tmpseqs))
    length(nseq) == 1 || error("sequences on different length")
    nseq = nseq[1]
    @show nseq
    for j in 1:length(taxa)
        ind = findall(isequal(taxa[j]), tmptaxa)
        length(ind) > 1 && error("taxa found multiple times")
        if length(ind) == 0 ##taxa not in this file
            seqs[j] = string(seqs[j], join(fill('-',nseq)))
        else
            ind = ind[1]
            seqs[j] = string(seqs[j], tmpseqs[ind])
        end
    end
end

output = "concatenated.fasta"
out = open(string(datafolder,output),"w")

for j in 1:n
    write(out, string(">",taxa[j],"\n",seqs[j],"\n"))
end
close(out)
