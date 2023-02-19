function phi1_loglikelihood(phi1, Z, l3, e3, T, nt)

    Sigma3 = get_dist_matrix(1:T, l3, e3)
    ll = logpdf(MvNormal(zeros(T), Sigma3), phi1)[1]
    for t in 1:T
        Nt = length(findall(Z[t] .== 1))
        ll += Nt * phi1[t] - nt[t]*log(exp(phi1[t]) + 1)
    end

    return ll 
end 

function update_phi1(cur, dat)

    T = dat["T"]
    nt = dat["nt"]
    Z = cur["Z"]
    l3 = cur["l3"]
    e3 = cur["e3"]
    
    
    phi1cur = cur["phi1"]
    phi1pro = cur["phi1"] + vec(rand(MvNormal(zeros(T), 0.1 .* Matrix(Diagonal(ones(T)))),1))

    llcur = phi1_loglikelihood(phi1cur, Z, l3, e3, T, nt)
    llpro = phi1_loglikelihood(phi1cur, Z, l3, e3, T, nt)

    if log(rand(Uniform(0,1),1)[1]) < (llpro - llcur)
        return phi1pro 
    else
        return phi1cur 
    end

end 