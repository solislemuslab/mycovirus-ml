## From https://github.com/FluxML/model-zoo/blob/master/other/iris/iris.jl
## Claudia July 2019

using Flux
using Flux: crossentropy, normalise, onecold, onehotbatch
using Statistics
using DataFrames, CSV
using JLD
using Random
using LinearAlgebra
include("functions.jl")

## ----------------------------------
## - Output file: this is the same ever,
## as new rows are just appended
## - universal parameters
## ----------------------------------
outputfile = "../results/model-accuracy.csv"
seed = 1608949
modelnum = 6  ## see list of options in the modelfile
propTrain = 2/3 ## proportion of training data

## -----------------------------
## Which dataset and model?
## -----------------------------
datafolder = "../data/staph/"
modelfile = "model-staph.jl"
response = "y" ##only matters for pseudomona
maxentropy = true ## for sequences: use sites with max entropy only
nvar = 700 ## arbitrarily chosen most entropy-rich variants

datafolder = "../data/pseudomonas/images/"
modelfile = "model-pseudomonas-images.jl"
response = "carb" ## carb or toby (very unbalanced)
maxentropy = false

datafolder = "../data/pseudomonas/sequences/"
modelfile = "model-pseudomonas-sequences.jl"
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

klasses = sort(unique(labels))
onehot_labels = onehotbatch(labels, klasses)
## staph: 2×124 Flux.OneHotMatrix{Array{Flux.OneHotVector,1}}
## pseudomonas images: 2×327 Flux.OneHotMatrix{Array{Flux.OneHotVector,1}}


## ------------------------------------
## Split into training and test sets:
##    propTrain for training, 1-propTrain for test.
## ------------------------------------
n = size(features,2)
p = size(normed_features,1)

rng = MersenneTwister(seed)
train_indices = randsubseq(rng, collect(1:n), propTrain)
test_indices = setdiff(collect(1:n), train_indices)

X_train = normed_features[:, train_indices]
y_train = onehot_labels[:, train_indices]

X_test = normed_features[:, test_indices]
y_test = onehot_labels[:, test_indices]


## ----------------
## Neural network
## ----------------
include(modelfile)
model,loss,optimiser = whichModel(modelnum,p,2)

# Create iterator to train model over 110 epochs.
data_iterator = Iterators.repeated((X_train, y_train), 110)

println("Starting training...")
tstart = time_ns() # in nanoseconds
Flux.train!(loss, params(model), data_iterator, optimiser)
tend = time_ns() # in nanoseconds
telapsed = round(convert(Int64, tend-tstart) * 1e-9, digits=2) # in seconds
println("\nTime elapsed: $telapsed seconds")

# Evaluate trained model against test set.
accuracy(x, y) = mean(onecold(model(x)) .== onecold(y))
accuracy_score = accuracy(X_test, y_test)
println("\nAccuracy: $accuracy_score")


## Saving results to file (appending)
nn = maxentropy ? nvar : size(normed_features,1)
df = DataFrame(modelnum=modelnum, modelfile=modelfile, response=response,
               time=telapsed, accuracy = accuracy_score,
               maxentropy=maxentropy, nvar= nn)

if isfile(outputfile)
    CSV.write(outputfile,df,append=true)
else
    CSV.write(outputfile,df)
end

# Sanity check.
@assert accuracy_score > 0.8

