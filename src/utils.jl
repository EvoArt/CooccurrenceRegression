

# code adapted from https://github.com/JuliaStats/Distributions.jl/blob/master/src/univariate/discrete/noncentralhypergeometric.jl 
# to work with reverse-mode automatic differentiation software
function _FNCH_mode(ns,nf,n,ω)
    A = ω - 1
    B = n - nf - (ns + n + 2)*ω
    C = (ns + 1)*(n + 1)*ω
    return -2C / (B - sqrt(B^2 - 4A*C))
end

FNCH_mode(ns,nf,n,ω) = floor(Int, _FNCH_mode(ns,nf,n,ω))

function FNCH_pdf(ns,nf,n,ω, k)

    ω, _ = promote(ω, float(k))
    # if nonsense value, reject sample
    -100000000000<ω<100000000000 ? nothing : return zero(ω)
    l = max(0, n - nf)
    u = min(ns, n)
    η = _FNCH_mode(ns,nf,n,ω)
    isfinite(η) ? η=floor(Int, η) : return zero(ω)
    s = one(ω)
    fᵢ = one(ω)
    fₖ = one(ω)
    for i in (η + 1):u
        rᵢ = (ns - i + 1)*ω/(i*(nf - n  + i))*(n - i + 1)
        fᵢ *= rᵢ
        # break if terms no longer contribute to s
        sfᵢ = s + fᵢ
        if sfᵢ == s && i > k
            break
        end
        s = sfᵢ

        if i == k
            fₖ = fᵢ
        end
    end
    fᵢ = one(ω)
    for i in (η - 1):-1:l
        rᵢ₊ = (ns - i)*ω/((i + 1)*(nf - n  + i + 1))*(n - i)
        fᵢ /= rᵢ₊

        # break if terms no longer contribute to s
        sfᵢ = s + fᵢ
        if sfᵢ == s && i < k
            break
        end
        s = sfᵢ

        if i == k
            fₖ = fᵢ
        end
    end

    return fₖ/s
end

FNCH_logpdf(ns,nf,n,ω, k) = log(FNCH_pdf(ns,nf,n,ω, k))


