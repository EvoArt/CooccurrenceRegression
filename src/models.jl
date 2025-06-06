
@model function regression(Xs,k, n_obs,N,m, n_pred = size(Xs,1);priors =[4.0,1.0,1.0])
   
    # Priors
    α ~ Normal(0.0,priors[1])
    β ~ filldist(Normal(0.0,priors[2]),n_pred)
    λ ~ filldist(Normal(0.0,priors[3]),n_obs)
    lp = 0.0
    for j in 1:n_obs-1
        mB = m[j]
        for i in j+1:n_obs
            mA = m[i]
            if (0<mA<N) & (0<mB<N)
                μ = α+λ[i] +λ[j]
                for pred in 1:n_pred
                    μ += β[pred]*Xs[pred,i,j] 
                end
                lp+=  FNCH_logpdf(mA,N-mA,mB,exp(μ),k[i,j])
            end
        end
    end
    Turing.@addlogprob! lp
end

# User facing methods

# vector of X arrays
function cooccurrence_regression(Xs::Vector,Y,n_iter=1000;alpha_sd =4.0,beta_sd = 1.0,lambda_sd =1.0,n_adapts=100,delta=0.65)
    n_sample, n_obs = size(Y)
    k = zeros(Int64, n_obs,n_obs)
    for i in 2:n_obs
        for j in 1:i-1
            k[i,j] = sum(Y[:,i] .* Y[:,j])
        end
    end

    n_X = length(Xs)
    X_dims = size(X)[1]
    X_3D = Array{Float64}(undef,n_X,X_dims...)
    for i in 1:n_X
        X_3D[i,:,:] .= Xs[i]
    end

    params = X_3D,k, n_obs,n_sample, vec(sum(Y, dims = 1))

    inference_model = regression(params...;priors =[alpha_sd,beta_sd,lambda_sd])
    Turing.sample(inference_model, NUTS(n_adapts,delta; adtype=ADTypes.AutoMooncake(config=nothing)),n_iter)
end

# single X array
function cooccurrence_regression(X,Y,n_iter=1000;alpha_sd =4.0,beta_sd = 1.0,lambda_sd =1.0,n_adapts=100,delta=0.65) 
    cooccurrence_regression([X],Y,n_iter;alpha_sd,beta_sd,lambda_sd,n_adapts,delta)
end