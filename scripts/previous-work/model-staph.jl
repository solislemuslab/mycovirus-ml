## Neural network model for staph data
## needs to be read inside fit-nn.jl
## Claudia July 2019


## -------------------------------------------
## function to keep track of all models tested
## -------------------------------------------
function whichModel(num::Integer, insize, outsize)
    if num == 1
        model = Chain(
                      Dense(insize, div(insize,1000)),
                      Dense(div(insize,1000), outsize),
                      softmax
                     )

        loss(x, y) = crossentropy(model(x), y)

        # Gradient descent optimiser with learning rate 0.5.
        optimiser = Descent(0.005) ##need small number or loss is infinite

        return model,loss,optimiser

    elseif num == 2
        model = Chain(
                      Dense(insize, div(insize,100)),
                      Dense(div(insize,100), outsize),
                      softmax
                     )

        loss2(x, y) = crossentropy(model(x), y)

        # Gradient descent optimiser with learning rate 0.5.
        optimiser = Descent(0.005) ##need small number or loss is infinite

        return model,loss2,optimiser
    elseif num == 3
        model = Chain(
                      Dense(insize, 1),
                      Dense(1, outsize),
                      softmax
                     )

        loss3(x, y) = crossentropy(model(x), y)

        # Gradient descent optimiser with learning rate 0.5.
        optimiser = Descent(0.005) ##need small number or loss is infinite

        return model,loss3,optimiser
    elseif num == 4
        model = Chain(
                      Dense(insize, outsize),
                      softmax
                     )

        loss4(x, y) = crossentropy(model(x), y)

        # Gradient descent optimiser with learning rate 0.5.
        optimiser = Descent(0.005) ##need small number or loss is infinite

        return model,loss4,optimiser

    elseif num == 5
        model = Chain(
            Dense(insize,256),
            Dense(256,256, relu),
            BatchNorm(256, relu),
            Dropout(0.7),
            Dense(256, outsize, σ),
            softmax
            )

        loss5(x, y) = crossentropy(model(x), y) + sum(norm, params(model))

        # Gradient descent optimiser with learning rate 0.5.
        optimiser = Descent(0.005) ##need small number or loss is infinite

        return model,loss5,optimiser

    elseif num == 6
        model = Chain(
            Dense(insize,256),
            Dense(256,256, relu),
            BatchNorm(256, relu),
            Dropout(0.5),
            Dense(256, 256, relu),
            BatchNorm(256, relu),
            Dropout(0.5),
            Dense(256,478, relu),
            BatchNorm(478, relu),
            Dropout(0.5),
            Dense(478, outsize, σ),
            softmax
            )

        loss6(x, y) = crossentropy(model(x), y) + sum(norm, params(model))

        # Gradient descent optimiser with learning rate 0.5.
        optimiser = Descent(0.005) ##need small number or loss is infinite

        return model,loss6,optimiser
    end
end


## model = Chain(
##     Conv((2, 2), p=>div(p,100), relu),
##     x -> maxpool(x, (2,2)),
##     Conv((2, 2), div(p,100)=>div(p,1000), relu),
##     x -> maxpool(x, (2,2)),
##     Dense(div(p,1000), 2),
##     softmax
## )



##loss(x, y) = Flux.mse(model(x), y)
##loss(x, y) = crossentropy(model(x), y) ## fixit: change to account for unbalancedness
##loss(x, y) = crossentropy(softmax(model(x)), y)


## model = Chain(
##               Dense(p,div(p,100)),
##               Dense(div(p,100),2),
##               softmax
##               )

##               function loss(x, y)
##                   return sum((model(x) .- y).^2)
##               end

