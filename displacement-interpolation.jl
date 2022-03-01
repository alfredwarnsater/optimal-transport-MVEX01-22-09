using Plots
using Plots.PlotMeasures
using QuadGK
using LinearAlgebra

function cost_matrix(A, B, n)
    sigma, _ = quadgk(s -> exp(A*(1-s))*B*transpose(B)*exp(transpose(A)*(1-s)), 0, 1)
    sigma_inv = inv(sigma);
    cost(x_0, x_1) = transpose(x_1-exp(A)*x_0)*sigma_inv*(x_1-exp(A)*x_0)

    grid = range(0, 1, n)
    pts = hcat(repeat(grid, inner = n), repeat(grid, outer = n))
    C = [cost([pts[i, 1]; pts[i, 2]], [pts[j, 1]; pts[j, 2]]) for i in 1:(n*n), j in 1:(n*n)]

    return C
end

function compute_interpolation(mu_1, mu_T, m, A, B, epsilon, tol)
    T = m + 2
    n = size(mu_1, 1)
    C = cost_matrix(A, B, n)
    
    mu_1 = reshape(mu_1, :, 1)
    mu_T = reshape(mu_T, :, 1) 
    
    K = exp.(-C / epsilon)
    u = ones(T, n*n)

    function P(tau)
        function phi_hat(tau)
            ans = 1
            for i in reverse(1:(tau-1))
                ans = ans * transpose(K) * (i != 1 ? Diagonal(u[i, :]) : u[i, :]) 
            end
            return ans
        end

        function phi(tau)
            ans = 1
            for i in (tau+1):T
                ans = ans * K * (i != T ? Diagonal(u[i, :]) : u[i, :])
            end
            return ans
        end

        return u[tau, :] .* phi_hat(tau) .* phi(tau)
    end

    u_prev = Inf * ones(size(u))
    while any(abs.(u - u_prev) .> tol)
        u_prev = copy(u)
        u[1, :] = u[1, :] .* mu_1 ./ P(1)
        u[T, :] = u[T, :] .* mu_T ./ P(T)
    end

    return [P(i) for i in 1:T]
end

function plot_results(data, clsc)
    plot_array = []
    for i in 1:size(data, 1)
        dist = reshape(data[i], isqrt(size(data[i], 1)), :)
        push!(plot_array, heatmap(dist, showaxis=false, ticks=false, legend=false, framestyle = :box,
            clims = (0, maximum(maximum.(data))), c  = clsc))
    end
    plot(plot_array..., size = (1000, 1000), link = :both)
end

function main(f, g, n, m, epsilon, tol)
    A = [0 0; 0 0]
    B = [1 0; 0 1]
    
    grid = range(0, 1, n)
    mu_1 = [f(x, y) for x in grid, y in grid]
    mu_T = [g(x, y) for x in grid, y in grid]
    mu_T = mu_T / sum(mu_T) * sum(mu_1)

    data = compute_interpolation(mu_1, mu_T, m, A, B, epsilon, tol)
    return data
end

x_0 = 0.5
y_0 = 0.5
r = 0.5
sigma_X = 0.1
sigma_Y = 0.1

f(x, y) = exp(-((x-x_0)^2/(2*sigma_X^2)+(y-y_0)^2/(2*sigma_Y^2)))

g(x, y) = exp(-((x-x_0+r)^2/(2*sigma_X^2)+(y-y_0+r)^2/(2*sigma_Y^2))) +
    + exp(-((x-x_0-r)^2/(2*sigma_X^2)+(y-y_0-r)^2/(2*sigma_Y^2)))

n = 50
m = 23
epsilon = 0.001
tol = 0.001

data = main(f, g, n, m, epsilon, tol)
plot_results(data, :jet1)




