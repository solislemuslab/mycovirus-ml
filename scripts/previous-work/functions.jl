## Auxiliary functions for fit-nn.jl
## Claudia August 2019


## function to remove invariant features in
## normed_features
## returns the new normed_features
## (for some reason I cannot change in place)
function removeInvariantFeatures(normed_features)
    ind = findall(x->!isnan(x),sum(normed_features,dims=2))
    normed_features = normed_features[map(x->x[1],ind),:]
    ss = sum(normed_features,dims=2)
    sum(isnan.(ss)) > 0 && error("constant features in the data")
    return normed_features
end


## extract the nvar most entropy-rich features
function extractMaxEntropyFeatures(features, nvar)
    ee = fill(0.0, size(features,1))
    for i in 1:size(features,1)
        p = proportionmap(features[i,:])
        ee[i] = entropy(values(p))
    end
    dd = DataFrame(entropy=ee, id = 1:length(ee))
    sort!(dd, :entropy, rev=true)
    features = features[dd[:id][1:nvar],:]
    return features
end
