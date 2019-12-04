## Random Forest for pseudomonas sequence data
## needs to be read inside fit-randForest.jl
## Claudia August 2019


## -------------------------------------------
## function to keep track of all models tested
## -------------------------------------------
function whichModel(num::Integer)
    if num == 1
        model = Dict("n_subfeatures" => -1,
                     "n_folds" => 2, ##should be 1 to make comparable to NN, but error
                     "n_trees" => 10,
                     "partial_sampling" => 0.7,
                     "max_depth" => -1,
                     "min_samples_leaf" => 5,
                     "min_samples_split" => 2,
                     "min_purity_increase" => 0.0)

        return model
    end
end
