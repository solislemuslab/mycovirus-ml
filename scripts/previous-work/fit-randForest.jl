## From https://github.com/bensadeghi/DecisionTree.jl
## Claudia August 2019

using Flux: normalise
using DecisionTree
using Statistics
using DataFrames, CSV
using JLD
using Random
include("functions.jl")

## ----------------------------------
## - Output file: this is the same ever,
## as new rows are just appended
## - universal parameters
## ----------------------------------
outputfile = "../results/model-accuracy.csv"
seed = 0841107
modelnum = 1  ## see list of options in the modelfile

## -----------------------------
## Which dataset and model?
## -----------------------------
datafolder = "../data/staph/"
modelfile = "model-staph-randForest.jl"
response = "y" ##only matters for pseudomona
maxentropy = true ## for sequences: use sites with max entropy only
nvar = 700 ## arbitrarily chosen most entropy-rich variants

datafolder = "../data/pseudomonas/images/"
modelfile = "model-pseudo-images-randForest.jl"
response = "carb" ## carb or toby (very unbalanced)
maxentropy = false

datafolder = "../data/pseudomonas/sequences/"
modelfile = "model-pseudo-seq-randForest.jl"
response = "carb" ## carb or toby (very unbalanced)
maxentropy = true
nvar = 700 ## arbitrarily chosen most entropy-rich variants

## -----------------------------
## Reading input JLD files
## -----------------------------

if modelfile == "model-pseudomonas-images.jl"
    include("data-preprocess-pseudomonas-image.jl")
    labels = load(string(datafolder,"labels-",response,".jld"), response) ## here we chose the response: carb or toby
elseif modelfile == "model-staph.jl"
    features = load(string(datafolder,"features.jld"), "features")
    labels = load(string(datafolder,"labels.jld"), "labels")
elseif modelfile == "model-pseudomonas-sequences.jl"
    features = load(string(datafolder,"features.jld"), "features")
    labels = load(string(datafolder,"labels-",response,".jld"), response) ## here we chose the response: carb or toby
end

## we want to choose the nvar with greater entropy (nvar chosen arbitrarily)
if maxentropy
    using StatsBase: entropy, proportionmap
    features = extractMaxEntropyFeatures(features, nvar)
end

## -----------------------------
## Normalizing features matrix
## -----------------------------
normed_features = normalise(features, dims=2)

## checking for invariant features:
normed_features = removeInvariantFeatures(normed_features)

## for DecisionTree we need nxp:
normed_features = normed_features'
normed_features = normed_features[:,:]

labels   = string.(labels)

## -----------------------------------
## Fitting random forest
## -----------------------------------
include(modelfile)
model = whichModel(modelnum)
rng = MersenneTwister(seed)

println("Starting training...")
tstart = time_ns() # in nanoseconds
accuracy = nfoldCV_forest(labels, features,
                          model["n_folds"],
                          model["n_subfeatures"],
                          model["n_trees"],
                          model["partial_sampling"],
                          model["max_depth"],
                          model["min_samples_leaf"],
                          model["min_samples_split"],
                          model["min_purity_increase"],
                          rng=rng)
tend = time_ns() # in nanoseconds
telapsed = round(convert(Int64, tend-tstart) * 1e-9, digits=2) # in seconds
println("\nTime elapsed: $telapsed seconds")
println("\nAccuracy: $(mean(accuracy))")

## Saving results to file (appending)
nn = maxentropy ? nvar : size(normed_features,2)
df = DataFrame(modelnum=modelnum, modelfile=modelfile, response=response,
               time=telapsed, accuracy = mean(accuracy),
               maxentropy=maxentropy, nvar= nn)

if isfile(outputfile)
    CSV.write(outputfile,df,append=true)
else
    CSV.write(outputfile,df)
end

# Sanity check.
@assert mean(accuracy) > 0.8

