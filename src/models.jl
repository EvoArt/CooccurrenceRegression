
@model function regression(Xs,k, n_obs,N,m, n_X = size(Xs,1);priors =[4.0,1.0,1.0])
   
    # Priors
    α ~ Normal(0.0,priors[1])
    β ~ filldist(Normal(0.0,priors[2]),n_X)
    λ ~ filldist(Normal(0.0,priors[3]),n_obs)
    lp = 0.0
    for j in 1:n_obs-1
        mB = m[j]
        for i in j+1:n_obs
            mA = m[i]
            if (0<mA<N) & (0<mB<N)
                μ = α+λ[i] +λ[j]
                for pred in 1:n_X
                    μ += β[pred]*Xs[pred,i,j] 
                end
                lp+=  FNCH_logpdf(mA,N-mA,mB,exp(μ),k[i,j])
            end
        end
    end
    Turing.@addlogprob! lp
end

function sub_regression(Xs,k,n_obs,N,m,n_X,β,λ,α)
    lp = 0.0
    for j in 1:n_obs-1
        mB = m[j]
        for i in j+1:n_obs
            mA = m[i]
            if (0<mA<N) & (0<mB<N)
                μ = α+λ[i] +λ[j]
                for pred in 1:n_X
                    μ += β[pred]*Xs[pred,i,j] 
                end
                lp+=  FNCH_logpdf(mA,N-mA,mB,exp(μ),k[i,j])
            end
        end
    end
    return lp
end


@model function h_regression(P,N_obs,n_p,n_X;priors =[4.0,1.0,1.0,1.0,2.0,0.1])
   
    # Priors
    α ~ Normal(0.0,priors[1])
    λ ~ filldist(Normal(0.0,priors[2]),N_obs)
    βμ ~ filldist(Normal(0.0,priors[3]),n_X)
    βσ ~ filldist(Gamma(priors[4],priors[5]),n_X)
    β ~ filldist(Normal(0.0,1.0),n_X,n_p)
    λinds = 0:0
    for i in eachindex(P)
        p=P[i]
        λinds = (λinds[end]+1):(λinds[end]+p[3])
        Turing.@addlogprob!  sub_regression(p...,n_X,(β[:,i] .* βσ) .+ βμ,λ[λinds],α)
    end
end

# User facing methods

# vector of X arrays
function cooccurrence_regression(Xs::Vector,Y,n_iter::Real=1000;drop_lambda=false,alpha_sd =4.0,beta_sd = 1.0,lambda_sd =1.0,n_adapts=100,delta=0.65)
    n_sample, n_obs = size(Y)
    k = zeros(Int64, n_obs,n_obs)
    for i in 2:n_obs
        for j in 1:i-1
            k[i,j] = sum(Y[:,i] .* Y[:,j])
        end
    end

    n_X = length(Xs)
    X_dims = size(Xs[1])
    X_3D = Array{Float64}(undef,n_X,X_dims...)
    for i in 1:n_X
        X_3D[i,:,:] .= Xs[i]
    end

    params = X_3D,k, n_obs,n_sample, vec(sum(Y, dims = 1))

    inference_model = regression(params...;priors =[alpha_sd,beta_sd,lambda_sd])
    chn = Turing.sample(inference_model, NUTS(n_adapts,delta; adtype=ADTypes.AutoMooncake(config=nothing)),n_iter)
    
    drop_lambda ? chn[["α",["β[$(i)]" for i in 1:n_X]...]] : chn
end

# single X array
function cooccurrence_regression(X,Y,n_iter::Real=1000;drop_lambda=false,alpha_sd =4.0,beta_sd = 1.0,lambda_sd =1.0,n_adapts=100,delta=0.65) 
    cooccurrence_regression([X],Y,n_iter;drop_lambda,alpha_sd,beta_sd,lambda_sd,n_adapts,delta)
end



function cooccurrence_regression(Xs::Vector,Y,g::Vector,n_iter::Real=1000;drop_lambda=false,alpha_sd =4.0,lambda_sd =1.0,
                                mu_sd =1.0,sigma_params=(2.0,0.1),n_adapts=100,delta=0.65)
    
    n_X = length(Xs)
    X_dims = size(Xs[1])
    X_3D = Array{Float64}(undef,n_X,X_dims...)
    for i in 1:n_X
        X_3D[i,:,:] .= Xs[i]
    end
    P= []
    for group in unique(sort(g))
        inds = g .== group
        y = Y[:,inds]
        x_3d = X_3D[:,inds,inds]
        n_sample, n_obs = size(y)
        k = zeros(Int64, n_obs,n_obs)
        for i in 2:n_obs
            for j in 1:i-1
                k[i,j] = sum(y[:,i] .* y[:,j])
            end
        end
        params = x_3d,k, n_obs,n_sample, vec(sum(y, dims = 1))
        push!(P,params)
    end
    P = Tuple.(P)

    inference_model = h_regression(P,size(Y,2),length(P),n_X; priors =[alpha_sd,lambda_sd,mu_sd,sigma_params...])
    chn = Turing.sample(inference_model, NUTS(n_adapts,delta; adtype=ADTypes.AutoMooncake(config=nothing)),n_iter)
    
    drop_lambda ? chn[["α",["β[$(i)]" for i in 1:n_X]...]] : chn
end
