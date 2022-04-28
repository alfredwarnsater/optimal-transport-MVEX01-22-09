
using QuadGK
using LinearAlgebra

include("displacement-interpolation-general.jl")
include("plotting.jl")

epsilon = 0.5
tol = 0.001
maze3 = [
    1 1 0 0 0
    0 1 0 0 0
    0 1 1 0 0
    0 0 1 1 1
    0 0 0 0 1
]

maze3 = [
    1 1 1 1 1 1
    1 1 1 1 1 1
    1 1 0 0 1 1
    1 1 0 0 1 1
    1 1 1 1 1 1
    1 1 1 1 1 1
]

maze = [
    1 1 0 0 0 0 0 0 0 0 0 
    0 1 1 1 1 1 1 1 1 1 0 
    0 1 0 1 0 0 0 1 0 1 0 
    0 1 0 1 0 1 0 1 0 1 0 
    0 1 0 0 0 1 0 0 0 1 0 
    0 1 1 1 1 1 0 1 1 1 0 
    0 0 0 0 0 0 0 1 0 1 0 
    0 1 1 1 1 1 0 1 0 1 0 
    0 0 0 0 0 1 0 1 0 0 0 
    0 1 1 1 1 1 1 1 1 1 1
    0 0 0 0 0 0 0 0 0 0 1
]

L = 25
n = size(maze, 1)

pts = hcat(repeat(1:n, inner = n), repeat(1:n, outer = n))

function cost(x_1, y_1, x_2, y_2)
    if x_1 == x_2 && y_1 == y_2
        return 0
    elseif abs(x_1-x_2) + abs(y_1-y_2) == 1
        return 1
    else
        return Inf
    end
end

C = [cost(pts[i, 1], pts[i, 2], pts[j, 1], pts[j, 2]) for i in 1:(n*n), j in 1:(n*n)]

mu = zeros(L, n*n)

mu_start = zeros(n*n)
mu_start[1] = 1
mu_end = zeros(n*n)
mu_end[n*n] = 1

mu[1, :] = mu_start
mu[L, :] = mu_end

for l = 2:(L-1)
    mu[l, :] = reshape(maze, 1, :)
end

types = fill('<', L)
types[1] = '='
types[L] = '='
mu = mu
data = compute_interpolation(C, mu, types, epsilon, tol)

maze = 

plot_results(data, maze)