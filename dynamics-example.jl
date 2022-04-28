using QuadGK
using LinearAlgebra
#using CairoMakie
using Plots
#include("plotting.jl")
include("displacement-interpolation-general.jl")

function cost_matrix(A, B, N)
    sigma, _ = quadgk(s -> exp(A*(1-s))*B*transpose(B)*exp(transpose(A)*(1-s)), 0, 1)
    sigma_inv = inv(sigma);
    cost(x_0, x_1) = transpose(x_1-exp(A)*x_0)*sigma_inv*(x_1-exp(A)*x_0)

    grid = range(0, 1, N)
    pts = hcat(repeat(grid, inner = N), repeat(grid, outer = N))
    C = [cost([pts[i, 1]; pts[i, 2]], [pts[j, 1]; pts[j, 2]]) for i in 1:(N*N), j in 1:(N*N)]

    return C
end

function plot_results(data)
    plot_array = []
    for i in 1:size(data, 1)
        dist = reshape(data[i], isqrt(size(data[i], 1)), :)
        push!(plot_array, heatmap(dist, showaxis=false, ticks=false, legend=false, framestyle = :box,
            clims = (0, maximum(maximum.(data))), c  = :jet))
    end
    plot(plot_array..., size = (1000, 1000), link = :both)
end

function gen_obstacle(N, L)

    h = 1 / N
    T = 100

    r = 0.2
    x_curve(t) = r*cos.(t)
    y_curve(t) = 0.5*r*sin.(t)
    x_start, x_end = 0, 1
    y_start, y_end = 0, 1
    
    inside(x) = all((1 .<= x .<= N))
    rotate(p, theta) = [p[1]*cos(theta)-p[2]*sin(theta), p[1]*sin(theta)+p[2]*cos(theta)]
    snap(x) = round.(Int32, fld.(x, h)) .+ 1

    function transform(pts, t_vec, r_angle)
        t_pts = zeros(size(pts))
        for i in range(1, T)
            t_pts[i, :] = rotate(pts[i, :], r_angle) + t_vec[:]
        end
        return t_pts
    end

    function gen_hollow_obs(pts)
        hollow_obs = zeros(N, N)
        pts = snap(pts)
        for i in range(1, T)
            if inside(pts[i, :])
                hollow_obs[pts[i, :]...] = 1
            end
        end
        return hollow_obs
    end

    function fill_cavity(hollow_obs, p)
        obs = copy(hollow_obs)
        function fcr(p)
            if !inside(p) || obs[p...] == 1
                return
            else
                obs[p...] = 1
            end
            fcr(p .+ [1, 0])
            fcr(p .+ [0, 1])
            fcr(p .+ [-1, 0])
            fcr(p .+ [0, -1])
        end
        fcr(p)
        return obs
    end

    ts = range(0, 2*pi, T)
    pts = [x_curve(ts) y_curve(ts)]
    xs = range(x_start, x_end, L)
    ys = range(y_start, y_end, L)
    rs = range(0, 2*pi, L)

    obstacle = zeros(L, N, N)

    for l in range(1, L) 
        c_pt = [xs[l], ys[l]]
        t_pts = transform(pts, c_pt, rs[l])
        hollow_obs = gen_hollow_obs(t_pts)
        obstacle[l, :, :] = fill_cavity(hollow_obs, snap(c_pt))
    end

    return obstacle

end

L = 49
N = 100
epsilon = 0.01
tol = 0.01

A = [
    0 0
    0 0
]
B = [
    1 0
    0 1
]

x_0 = 0.5
y_0 = 0.5
r = 0.5
sigma_X = 0.1
sigma_Y = 0.1

f(x, y) = exp(-((x-x_0)^2/(2*sigma_X^2)+(y-y_0)^2/(2*sigma_Y^2)))

g(x, y) = exp(-((x-x_0+r)^2/(2*sigma_X^2)+(y-y_0+r)^2/(2*sigma_Y^2))) +
    + exp(-((x-x_0-r)^2/(2*sigma_X^2)+(y-y_0-r)^2/(2*sigma_Y^2)))

C = cost_matrix(A, B, N)
grid = range(0, 1, N)
mu_start = [f(x, y) for x in grid, y in grid]
total_mass = sum(mu_start)
mu_end = [g(x, y) for x in grid, y in grid]
mu_end = mu_end / sum(mu_end) * total_mass

obstacle = gen_obstacle(N, L-2)
obstacle = (obstacle .- 1) * (-total_mass)

mu = zeros(L, N*N)
types = fill('-', L)
mu[1, :] = reshape(mu_start, 1, :)
types[1] = '='
mu[L, :] = reshape(mu_end, 1, :)
types[L] = '='

for l in 1:(L-2)
    mu[l+1, :] = reshape(obstacle[l, :, :], 1, :)
    types[l+1] = '<'
end

data = @time @eval compute_interpolation(C, mu, types, epsilon, tol)

plot_results(data)




