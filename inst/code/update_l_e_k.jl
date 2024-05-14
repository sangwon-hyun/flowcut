function l_e_loglikelihood(vec, l, e, al, bl, ae, be, T)
    """
    Calculate log-likelihood for lⱼ and eⱼ for j=1,2,3
    """
    ll = logpdf(Gamma(al, bl), l) + logpdf(Gamma(ae, be), e) + log(l) + log(e) # jacobian 
    vec = vec .+ 10e-7 .* ones(T)
    K = get_dist_matrix(1:T, 1:T, l, e) 
    K = stablizeMatrix(K)

    ll += logpdf(MvNormal(zeros(T), K), vec)[1]
    return ll 
end 

function update_l_e_k(cur, hyper, dat, k)
    """
    Jointly update l1 and e1 
    """
    T = dat["T"]
    mu = cur["mu"*string(k)]

    al = hyper["al"*string(k)]
    bl = hyper["bl"*string(k)]
    ae = hyper["ae"*string(k)]
    be = hyper["be"*string(k)]

    lcur = cur["l"*string(k)]
    ecur = cur["e"*string(k)]

    lpro = exp(log(lcur) + rand(Normal(0,0.1), 1)[1])
    epro = exp(log(ecur) + rand(Normal(0,0.1), 1)[1])

    llcur = l_e_loglikelihood(mu, lcur, ecur, al, bl, ae, be, T)
    llpro = l_e_loglikelihood(mu, lpro, epro, al, bl, ae, be, T)

    if log(rand(Uniform(0,1),1)[1]) < (llpro - llcur)
        return [lpro, epro] 
    else
        return [lcur, ecur] 
    end
    
end
