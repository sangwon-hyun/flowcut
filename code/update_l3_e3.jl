function update_l3_e3(cur, hyper, dat)
    T = dat["T"]
    phi1 = cur["phi1"]

    al = hyper["al"]
    bl = hyper["bl"]
    ae = hyper["ae"]
    be = hyper["be"]

    l3cur = cur["l3"]
    e3cur = cur["e3"]

    l3pro = exp(log(l3cur) + rand(Normal(0,0.1),1)[1])
    e3pro = exp(log(e3cur) + rand(Normal(0,0.1),1)[1])

    llcur = l_e_loglikelihood(phi1, l3cur, e3cur, al, bl, ae, be, T)
    llpro = l_e_loglikelihood(phi1, l3pro, e3pro, al, bl, ae, be, T)

    if log(rand(Uniform(0,1),1)[1]) < (llpro - llcur)
        return [l3pro, e3pro] 
    else
        return [l3cur, e3cur] 
    end
    
end