using Plots
using Plots.PlotMeasures
using JuMP
using GLPK
using Random, Distributions

include("transport-map-marginal-plot.jl")

sizeOfGrid = 20

a = pdf.(Normal(sizeOfGrid/2, 3/20*sizeOfGrid), 0:1:sizeOfGrid-1)
b = pdf.(Normal(sizeOfGrid/4, 1.5/20*sizeOfGrid), 0:1:sizeOfGrid-1)/3 + pdf.(Normal(sizeOfGrid*3/4, sizeOfGrid*1.5/20), 0:1:sizeOfGrid-1)/3

b = b .* sum(a) ./ sum(b)

n = length(a)

c(i, j) = abs(i - j)^2 / 32340.0

model = Model(GLPK.Optimizer)
@variable(model, 0 <= M[1:n, 1:n])
@objective(model, Min, sum(sum(c(i, j)*M[i, j] for j in 1:n) for i in 1:n))
@constraint(model, c1, M*ones(n) .== b)
@constraint(model, c2, M'*ones(n) .== a)
optimize!(model)

best_perm_mat = value.(M)

transport_map_marginal_plot(transpose(a), transpose(b), best_perm_mat)
