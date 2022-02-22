# Given the cost tensor C, index set gamma, marginals mu_j and regularization parameter epsilon this function computes
# the transport map M using multimarginal Sinkhorn iterations.

function compute_transport_plan_sinkhorn(C, gamma, mu, epsilon)

    J = ndims(C)
    K = exp.(-C / epsilon)
    #@show K
    N = size(C)[1] # Assuming all dimensions are of equal size.
    u = ones(J, N)

    function compute_U()
        U = ones(size(C))
        for i in CartesianIndices(U)
            U[i] = 1
            for j in gamma
                U[i] = U[i] * u[j, i[j]]
            end
        end
        #@show U
        return U
    end

    count = 1
    while count < 1000
        #for j in gamma
            U = compute_U()
            u[1, :] = u[1, :] .* mu[1, :] ./ reshape(sum(K .* U, dims=2), (3, 1))
            u[2, :] = u[2, :] .* mu[2, :] ./ reshape(sum(K .* U, dims=1), (3, 1))
            #display(K .* U)
        #end
        count = count + 1
    end

    return K .* compute_U()

end
C = [0 1 2; 1 0 1; 2 1 0]
C = C ./ maximum(C)
M = compute_transport_plan_sinkhorn(C, [1,2], [1 2 4; 4 2 1], 0.01)