"""
partitionTrainTest(data, y_column, train_perc=0.7)

Splits data (rows × columns) into train/test subsets.

# Arguments
- data: array or matrix
- y_column: integer, the column index of the target variable
- train_perc: fraction of data used for training (default 0.7)

# Returns
- X_train, y_train, X_test, y_test
"""
function partitionTrainTest(data, y_column, train_perc)
    n = size(data, 1)
    idx = shuffle(1:n)
    n_train = floor(Int, train_perc * n)
    train_idx = view(idx, 1:n_train)
    test_idx = view(idx, n_train+1:n)

    X_train = Matrix(data[train_idx, Not(y_column)])
    y_train = Array(data[train_idx, y_column])
    X_test = Matrix(data[test_idx, Not(y_column)])
    y_test = Array(data[test_idx, y_column])

    return X_train, y_train, X_test, y_test
end





"""
    r2_score(y_true, y_pred) -> Float64

Compute the coefficient of determination R² between the true target values and predicted values.

R² is a statistical measure that indicates how well predictions approximate the actual data.
It is defined as:

    R² = 1 - (SS_res / SS_tot)

where:
- `SS_res` is the sum of squared residuals (errors),
- `SS_tot` is the total sum of squares (variance of the data).

An R² of 1 indicates perfect prediction; 0 indicates that the model does no better than predicting the mean.

# Arguments
- `y_true::Vector{<:Real}`: Actual target values.
- `y_pred::Vector{<:Real}`: Predicted target values.

# Returns
- `Float64`: The R² score.
"""
function r2_score(y_true, y_pred)
    ss_res = sum((y_true .- y_pred).^2)
    ss_tot = sum((y_true .- mean(y_true)).^2)
    return 1 - ss_res / ss_tot
end




"""
    random_forest_julia(X_train, y_train, X_test, y_test) -> Float64

Train and evaluate a random forest regressor using the `DecisionTree.jl` package.

# Arguments
- `X_train::Matrix{Float64}`: Feature matrix for training (features × samples).
- `y_train::Vector{Float64}`: Target vector for training.
- `X_test::Matrix{Float64}`: Feature matrix for testing (features × samples).


# n_subfeatures: number of features to consider at random per split (default: -1, sqrt(# features))
# n_trees: number of trees to train (default: 10)
# partial_sampling: fraction of samples to train each tree on (default: 0.7)
# max_depth: maximum depth of the decision trees (default: no maximum)
# min_samples_leaf: the minimum number of samples each leaf needs to have (default: 5)
# min_samples_split: the minimum number of samples in needed for a split (default: 2)
# min_purity_increase: minimum purity needed for a split (default: 0.0)
# keyword rng: the random number generator or seed to use (default Random.GLOBAL_RNG)
# multi-threaded forests must be seeded with an `Int`


# Returns
- `Surrogate{struct}`: model, prediction, loss (empty), r2 (empty).
"""
function random_forest_julia(X_train, y_train, X_test)
    m = build_forest(y_train, X_train)
    ŷ = apply_forest(m, X_test)
    return Surrogate(m, ŷ)
end




"""
    random_forest_sklearn(X_train, y_train, X_test, y_test) -> Float64

Train and evaluate a random forest regressor using the `ScikitLearn.jl` interface to Python's scikit-learn.

# Arguments
- `X_train::Matrix{Float64}`: Feature matrix for training (samples × features).
- `y_train::Vector{Float64}`: Target vector for training.
- `X_test::Matrix{Float64}`: Feature matrix for testing (samples × features).

# Returns
- `Surrogate{struct}`: model, prediction, loss (empty), r2 (empty).
"""

function random_forest_sklearn(X_train, y_train, X_test)
    m = RandomForestRegressor()
    ScikitLearn.fit!(m, X_train, y_train)
    ŷ = ScikitLearn.predict(m, X_test)
    return Surrogate(m, ŷ)
end




"""
    decision_tree_julia(X_train, y_train, X_test, y_test) -> Float64

Train and evaluate a single decision tree regressor using the `DecisionTree.jl` package.

# Arguments
- `X_train::Matrix{Float64}`: Feature matrix for training (features × samples).
- `y_train::Vector{Float64}`: Target vector for training.
- `X_test::Matrix{Float64}`: Feature matrix for testing (features × samples).

# Returns
- `Surrogate{struct}`: model, prediction, loss (empty), r2 (empty).

"""

function decision_tree_julia(X_train, y_train, X_test)
    m = build_tree(y_train, X_train)
    ŷ = apply_tree(m, X_test)
    return Surrogate(m, ŷ)
end




"""
    decision_tree_sklearn(X_train, y_train, X_test, y_test)

Trains and evaluates a decision tree regressor using ScikitLearn.jl.

# Arguments
- `X_train::Matrix`: Training feature matrix.
- `y_train::Vector`: Training target vector.
- `X_test::Matrix`: Test feature matrix.

# Returns
- `Surrogate{struct}`: model, prediction, loss (empty), r2 (empty).
"""

function decision_tree_sklearn(X_train, y_train, X_test)
    m = DecisionTreeRegressor()
    ScikitLearn.fit!(m, X_train, y_train)
    ŷ = ScikitLearn.predict(m, X_test)
    return Surrogate(m, ŷ)
end






"""
    linear_regression_sklearn(X_train, y_train, X_test, y_test)

Trains a linear regression model using ScikitLearn.jl and predicts outputs on the test set.

# Arguments
- `X_train::Matrix`: Training feature matrix (samples × features).
- `y_train::Vector`: Training target vector.
- `X_test::Matrix`: Test feature matrix.

# Returns
- `Surrogate{struct}`: model, prediction, loss (empty), r2 (empty).
"""

function linear_regression_sklearn(X_train, y_train, X_test)
    m = LinearRegression()
    ScikitLearn.fit!(m, X_train, y_train)
    ŷ = ScikitLearn.predict(m, X_test)
    return Surrogate(m, ŷ)
end






"""
    linear_regression_julia(X_train, y_train, X_test, y_test)

Fits a multiple linear regression model using `GLM.jl` on the training data and predicts target values for the test data.

# Arguments
- `X_train::Matrix{Float64}`: Feature matrix for training (samples × features).
- `y_train::Vector{Float64}`: Target vector for training.
- `X_test::Matrix{Float64}`: Feature matrix for testing.

# Returns
- `Surrogate{struct}`: model, prediction, loss (empty), r2 (empty).

"""


function linear_regression_julia(X_train, y_train, X_test)
    # Create DataFrame from feature matrix
    df = DataFrame(X_train, Symbol.("x", 1:size(X_train, 2)))
    features = names(df)
    # Add target variable
    df.y = y_train
    features_sym = Symbol.(features)  # convert each string to symbol
    formula = Term(:y) ~ sum(Term.(features_sym))
    m = GLM.lm(formula, df)
    ŷ = GLM.predict(m,  DataFrame(X_test, Symbol.("x", 1:size(X_test, 2))))
    return Surrogate(m, ŷ)
end





"""
    neural_network_model_flux(X_train, y_train, X_test, y_test; hidden_layer=1000, epochs=1000)

Train a feedforward neural network using Flux.jl for a regression task.
example taken from: https://apxml.com/courses/julia-for-machine-learning/chapter-6-julia-deep-learning-flux-jl/flux-jl-training-neural-networks

# Arguments
- `X_train::Matrix`: Training features matrix (samples × features).
- `y_train::Vector`: Target values for training data.
- `X_test::Matrix`: Test features matrix (samples × features).
- `y_test::Vector`: Target values for test data.
- `hidden_layer::Int=1000`: Number of neurons in the hidden layer (default: 1000).
- `epochs::Int=1000`: Number of training epochs (default: 1000).

# Returns
- `Surrogate{struct}`: model, prediction, loss, r2.

"""


function neural_network_model_flux(X_train, y_train, X_test, y_test; hidden_layer=1000, epochs=1000)

    x = permutedims(X_train)                    # shape: features × samples
    y = reshape(y_train, 1, :)                  # shape: 1 × samples

    m = Chain(
    Dense(size(X_train, 2), hidden_layer, relu),
    Dense(hidden_layer, 1)
)

    # track parameters
    θ = Flux.params(m)
    # select an optimizer
    α = 0.001 
    opt = ADAM(α)

    # train the model
    losses = []
    r2_list =  []

    loss(x, y) = Flux.Losses.mae(m(x), y)
    loss_fn(m, x, y) = Flux.mae(m(x), y) 

    train_loader = Flux.DataLoader((x, y), batchsize=32, shuffle=true)
    opt_state = Flux.setup(opt, m)
    
    for epoch in 1:epochs
        # train the model
        for (x_batch, y_batch) in train_loader
            current_loss, grads = Flux.withgradient(m) do m_in_grad
                loss_fn(m_in_grad, x_batch, y_batch)
            end
            Flux.update!(opt_state, m, grads[1])
        end

        # print report
        ŷ = m(permutedims(X_test))
        r2 = P2H_CapacityExpansion.r2_score(y_test, permutedims(ŷ))
        println("Epoch = $epoch : Training loss = $(loss(x, y)), R2 score : $(r2)")
        push!(r2_list, r2)
        push!(losses, loss(x, y))
    end 

    return Surrogate(m, ŷ)
end





"""
    simple_neural_network_sklearn(X_train, y_train, X_test, y_test; hidden_layer=1000, max_iter=1000)

Train a multi-layer perceptron regression model using the ScikitLearn.jl wrapper for Python's `MLPRegressor`.

# Arguments
- `X_train::Matrix`: Training features (samples × features).
- `y_train::Vector`: Training target values.
- `X_test::Matrix`: Test features.
- `y_test::Vector`: Test target values.
- `hidden_layer::Int=1000`: Number of neurons in the hidden layer (default: 1000).
- `max_iter::Int=1000`: Maximum number of training iterations (default: 1000).

# Returns
- `Surrogate{struct}`: model, prediction, loss, r2.

# Description
"""


function simple_neural_network_sklearn(X_train, y_train, X_test; hidden_layer=500, max_iter=1000) 

    m = MLPRegressor(hidden_layer_sizes=(size(X_train, 2), hidden_layer), max_iter=max_iter, random_state=42, verbose=true)

    ScikitLearn.fit!(m, X_train, y_train)

    # Predict
    ŷ = ScikitLearn.predict(m, X_test)

    return Surrogate(m, ŷ)
end



"""
    gaussian_process(X_train, y_train, X_test)

Train a Gaussian Process regression model using a Rational Quadratic kernel and predict on test data.

# Arguments
- `X_train::Matrix`: Training features (samples × features).
- `y_train::Vector`: Target values for training data.
- `X_test::Matrix`: Test features to predict on.

# Returns
- `Surrogate{struct}`: model, prediction, loss, r2.

"""


function gaussian_process(X_train, y_train, X_test)

    mZero = MeanZero()  
    kern =  RQ(0.0, 0.0, 0.0)
    
    gp = ScikitLearn.fit!(GPE(mean=mZero,kernel=kern, logNoise=-1.0), X_train, y_train)
    ŷ = ScikitLearn.predict(gp, X_test)

    return Surrogate(gp, ŷ)
end

"""
    svr_sklearn(X_train, y_train, X_test)

Train a Support Vector Regression (SVR) model using ScikitLearn.jl and make predictions on test data.

# Arguments
- `X_train::Matrix`: Feature matrix for training data (samples × features).
- `y_train::Vector`: Target values for training data.
- `X_test::Matrix`: Feature matrix for test data.

# Returns
- `Surrogate`: A custom struct containing:
    - `model`: The trained SVR model.
    - `prediction`: Predicted values on the test set.
"""

function svr_sklearn(X_train, y_train, X_test)

    # Create and fit the model
    m = SVR(kernel="rbf", C=1.0, epsilon=0.1)
    ScikitLearn.fit!(m, X_train, y_train)

    ŷ = ScikitLearn.predict(m, X_test)

    return Surrogate(m, ŷ)
end

