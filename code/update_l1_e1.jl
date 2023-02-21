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

function update_l1_e1(cur, hyper, dat)
    """
    Jointly update l1 and e1 
    """
    T = dat["T"]
    mu1 = cur["mu1"]

    al = hyper["al"]
    bl = hyper["bl"]
    ae = hyper["ae"]
    be = hyper["be"]

    l1cur = cur["l1"]
    e1cur = cur["e1"]

    l1pro = exp(log(l1cur) + rand(Normal(0,1), 1)[1])
    e1pro = exp(log(e1cur) + rand(Normal(0,1), 1)[1])

    llcur = l_e_loglikelihood(mu1, l1cur, e1cur, al, bl, ae, be, T)
    llpro = l_e_loglikelihood(mu1, l1pro, e1pro, al, bl, ae, be, T)

    if log(rand(Uniform(0,1),1)[1]) < (llpro - llcur)
        return [l1pro, e1pro] 
    else
        return [l1cur, e1cur] 
    end
    
end
