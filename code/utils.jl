function get_dist_matrix(mu, l, e)
    T = length(mu) 
    M = zeros(T,T) 
    for r in 1:T
        for c in r:T 
            M[r, c] = M[c, r] = e * exp(- (mu[r] - mu[c])^2 / (2 * l))
        end 
    end 
    return M 
end 


function svd2inv(M)

	X = svd(M)
	Minv = X.Vt' * Diagonal(1 ./ X.S) * X.U'
	Minv = (Minv + Minv')/2

	return Minv 
end