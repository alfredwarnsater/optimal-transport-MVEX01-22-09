using JuMP
using GLPK
using Random, Distributions
using Plots
using Plots.PlotMeasures

# Denna funktion plottar resultatet till exemplet som genereras nedan.
function transport_map_marginal_plot(marginal_1, marginal_2, transport_map, filename)
    margin = -20;
    bw = 0.8
    m1 = bar(transpose(marginal_1), orientation = :vertical, bar_width=bw, yflip=false, bottom_margin=margin*Plots.px, 
        showaxis=false, ticks=false, legend = false,
        color = "red", linecolor= "red")
    m2 = bar(transpose(marginal_2), orientation = :horizontal, bar_width=bw, xflip=true, right_margin=margin*Plots.px, 
        showaxis=false, ticks=false, legend = false,
        color = "blue", linecolor= "blue")
    hm = Plots.heatmap(transport_map, color = :greys, showaxis=false, ticks=false, legend = false, framestyle = :box)

    l = @layout[_ a; b c{0.7w, 0.7h}]

    display(Plots.plot(m1, m2, hm, size = (600, 600), layout = l, link = :both))
    savefig(filename)
end

# Denna funktion genererar och löser ett exempel av bimarginell optimaltransport.
function bimarginal_transport(N)
    a = pdf.(Normal(N/2, 3/20*N), 0:1:N-1)
    b = pdf.(Normal(N/4, 1.5/20*N), 0:1:N-1)/3 + pdf.(Normal(N*3/4, N*1.5/20), 0:1:N-1)/3
    b = b .* sum(a) ./ sum(b)

    n = length(a)

    c(i, j) = abs(i - j)^2

    model = Model(GLPK.Optimizer)
    @variable(model, 0 <= M[1:n, 1:n])
    @objective(model, Min, sum(sum(c(i, j)*M[i, j] for j in 1:n) for i in 1:n))
    @constraint(model, c1, M*ones(n) .== b)
    @constraint(model, c2, M'*ones(n) .== a)
    optimize!(model)
    best_perm_mat = value.(M)

    transport_map_marginal_plot(transpose(a), transpose(b), best_perm_mat, "plots/bimarginal_transport.pdf")
end

bimarginal_transport(20)



