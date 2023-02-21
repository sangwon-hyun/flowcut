function get_dist_matrix(vec1, vec2, l, e)
    """
    Calculate GP covariance matrix 
    """

    N1 = length(vec1) 
    N2 = length(vec2) 

    M = zeros(N1,N2) 
    for r in 1:N1
        for c in 1:N2 
            M[r, c] = e * exp(- (vec1[r] - vec2[c])^2 / (2 * l))
        end 
    end  

    return M 
end 


function stablizeMatrix(K)
    """
    Used for fixing non-pos-def matrix caused by numerical issues 
    """
    if isposdef(K)
        return K
    end 

    while !isposdef(K)
        K = K + Diagonal(10e-7 .* ones(size(K)[1]))
    end
    return K
end


function svd2inv(M)
    """
    Matrix inversion using singular value decomposition
    """
    
    N = size(M)[1]

    if !isposdef(M)
        M = M + Matrix(Diagonal(10e-7 .* ones(N))) 
    end 
	X = svd(M)
	Minv = X.Vt' * Diagonal(1 ./ X.S) * X.U'
	Minv = (Minv + Minv')/2

	return Minv 
end