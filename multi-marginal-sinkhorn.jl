# Given the cost tensor C, index set gamma, marginals mu, regularization parameter epsilon and convergence tolerance tol
# this function computes the transport map M using Sinkhorn iterations. 

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

using Random, Distributions
include("transport-map-marginal-plot.jl")

# Pretty example

#= a = pdf.(Normal(10, 3), 0:1:20)
b = pdf.(Normal(5, 1.5), 0:1:20)/3 + pdf.(Normal(15, 1.5), 0:1:20)/3
#a = a .* 1e5
b = b .* sum(a) ./ sum(b)
a = transpose(a)
b = transpose(b) =#

a = [1 2 3]
b = [3 2 1]

n = length(a)
C = zeros(n, n)
for i in 1:n
    for j in 1:n
        C[i, j] = abs(i - j)
    end
end 

#C = [0 1 2; 1 0 1; 2 1 0]
#C = C ./ 1e5 #
M = compute_transport_plan_sinkhorn(C, [1, 2], [a; b], 0.01, 1e-5)

#M

transport_map_marginal_plot(a, b, M)
