function update_l2_e2(cur, hyper, dat)
    T = dat["T"]
    mu2 = cur["mu2"]

    al = hyper["al"]
    bl = hyper["bl"]
    ae = hyper["ae"]
    be = hyper["be"]

    l2cur = cur["l2"]
    e2cur = cur["e2"]

    l2pro = exp(log(l2cur) + rand(Normal(0,0.1),1)[1])
    e2pro = exp(log(e2cur) + rand(Normal(0,0.1),1)[1])

    llcur = l_e_loglikelihood(mu2, l2cur, e2cur, al, bl, ae, be, T)
    llpro = l_e_loglikelihood(mu2, l2pro, e2pro, al, bl, ae, be, T)

    if log(rand(Uniform(0,1),1)[1]) < (llpro - llcur)
        return [l2pro, e2pro] 
    else
        return [l2cur, e2cur] 
    end
    
end