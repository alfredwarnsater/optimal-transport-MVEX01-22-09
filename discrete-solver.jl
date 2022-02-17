using Plots
using Plots.PlotMeasures
using JuMP
using GLPK
using Random, Distributions

gr()

a = pdf.(Normal(10, 3), 0:1:20)
b = pdf.(Normal(5, 1.5), 0:1:20)/3 + pdf.(Normal(15, 1.5), 0:1:20)/3

a = a / sum(a) * sum(b);

n = length(a)
m = maximum(a)

c(i, j) = abs(i - j)

model = Model(GLPK.Optimizer)
@variable(model, 0 <= M[1:n, 1:n])
@objective(model, Min, sum(sum(c(i, j)*M[i, j] for j in 1:n) for i in 1:n))
@constraint(model, c1, M*ones(n) .== b)
@constraint(model, c2, M'*ones(n) .== a)
optimize!(model)
best_perm_mat = value.(M)

transport_map = heatmap(best_perm_mat, color = :greys, showaxis=false, ticks=false, legend = false, framestyle = :box)
margin = -20;
marginal1 = bar(a, orientation = :vertical, bar_width=1, yflip=false, bottom_margin=margin*px, 
showaxis=false, ticks=false, legend = false, ylims = (0, m), color = "red", linecolor= "red")

marginal2 = bar(b, orientation = :horizontal, bar_width=1, xflip=true, right_margin=margin*px, 
showaxis=false, ticks=false, legend = false, xlims = (0, m), color = "blue", linecolor= "blue")

l = @layout[_ a; b c{0.7w, 0.7h}]

plot(marginal1, marginal2, transport_map, size = (600, 600), layout = l, link = :both)



# heights = [0.3, 0.7], widths = [0.3, 0.7]
# empty_white = plot(legend = false, grid = false, foreground_color_subplot=:white, margin = 0px)