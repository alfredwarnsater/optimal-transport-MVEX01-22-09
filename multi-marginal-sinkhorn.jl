# Given the cost tensor C, index set gamma, marginals mu, regularization parameter epsilon and convergence tolerance tol
# this function computes the transport map M using multimarginal Sinkhorn iterations. 

function compute_transport_plan_sinkhorn(C, gamma, mu, epsilon, tol)

    J = ndims(C)
    K = exp.(-C / epsilon)
    N = size(C)[1]  # Assuming all dimensions are of equal size.
    u = ones(J, N)

    function compute_U()
        U = ones(size(C))
        for i in CartesianIndices(U)
            U[i] = 1
            for j in gamma
                U[i] = U[i] * u[j, i[j]]
            end
        end
        return U
    end

    u_prev = Inf*ones(size(u))
    while any(abs.(u - u_prev) .> tol)  # Not sure if this is the best way to check for convergence, seems expensive.
        u_prev = copy(u)
        for j in gamma
            all_but_j = filter((x) -> (x != j), 1:J)
            u[j, :] = u[j, :] .* mu[j, :] ./ reshape(sum(K .* compute_U(), dims=all_but_j), (N, 1))
        end
    end
    return K .* compute_U()

end
C = [0 1 2; 1 0 1; 2 1 0]
C = C ./ maximum(C)
M = compute_transport_plan_sinkhorn(C, [1,2], [1 2 4; 4 2 1], 0.01, 0.1)