function phi1_loglikelihood(phi1, Z, l3, e3, T, nt)
    """
    Calculate log-likelihood of phi1 vector
    """

    Sigma3 = get_dist_matrix(1:T, 1:T, l3, e3)
    ll = logpdf(MvNormal(zeros(T), Sigma3), phi1)
    for t in 1:T
        for i in 1:nt[t]
            if Z[t][i] == 1
                ll = ll + phi1[t] - log(exp(phi1[t]) + 1)
            else
                ll = ll - log(exp(phi1[t]) + 1)
            end
        end
    end

    return ll 
end 

function update_phi1(cur, dat)
    """
    Update phi1 vector 
    """

    T  = dat["T"]
    nt = dat["nt"]

    Z  = cur["Z"]
    l3 = cur["l3"]
    e3 = cur["e3"]
    
    phi1cur = cur["phi1"]
    phi1pro = cur["phi1"] + rand(MvNormal(zeros(T), Matrix(Diagonal(0.01 .* ones(T)))),1)[:,1]

    llcur = phi1_loglikelihood(phi1cur, Z, l3, e3, T, nt)
    llpro = phi1_loglikelihood(phi1pro, Z, l3, e3, T, nt)

    if log(rand(Uniform(0,1),1)[1]) < (llpro - llcur)
        return phi1pro 
    else
        return phi1cur 
    end

end 