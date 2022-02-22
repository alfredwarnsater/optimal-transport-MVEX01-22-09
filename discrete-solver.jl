using Plots
using Plots.PlotMeasures
using JuMP
using GLPK
using Random, Distributions

include("transport-map-marginal-plot.jl")

a = pdf.(Normal(10, 3), 0:1:20)
b = pdf.(Normal(5, 1.5), 0:1:20)/3 + pdf.(Normal(15, 1.5), 0:1:20)/3

a = a / sum(a) * sum(b);

n = length(a)

c(i, j) = abs(i - j)

model = Model(GLPK.Optimizer)
@variable(model, 0 <= M[1:n, 1:n])
@objective(model, Min, sum(sum(c(i, j)*M[i, j] for j in 1:n) for i in 1:n))
@constraint(model, c1, M*ones(n) .== b)
@constraint(model, c2, M'*ones(n) .== a)
optimize!(model)
best_perm_mat = value.(M)

transport_map_marginal_plot(a, b, best_perm_mat)
