function update_s_p_k(cur, hyper, dat, k)
    """
    Jointly update s and psi
    """
    T = dat["T"]
    phi = cur["phi"*string(k)]

    as = hyper["as"*string(k)]
    bs = hyper["bs"*string(k)]
    ap = hyper["ap"*string(k)]
    bp = hyper["bp"*string(k)]

    s_cur = cur["s"*string(k)]
    p_cur = cur["p"*string(k)]

    s_pro = exp(log(s_cur) + rand(Normal(0,0.1),1)[1])
    p_pro = exp(log(p_cur) + rand(Normal(0,0.1),1)[1])

    llcur = l_e_loglikelihood(phi, s_cur, p_cur, as, bs, ap, bp, T)
    llpro = l_e_loglikelihood(phi, s_pro, p_pro, as, bs, ap, bp, T)

    if log(rand(Uniform(0,1),1)[1]) < (llpro - llcur)
        return [s_pro, p_pro] 
    else
        return [s_cur, p_cur] 
    end
    
end