using JuMP
using GLPK
using Plots

# Denna funktion beräknar och plottar en lösning till matchningsproblemet på enhetskvadraten 
# med linjärprogrammering.
function matching_lp(n, filename)
    x0 = rand(Float64, n)
    y0 = rand(Float64, n)
    x1 = rand(Float64, n)
    y1 = rand(Float64, n)
    c(i, j) = (x0[i]-x1[j])^2 + (y0[i]-y1[j])^2

    model = Model(GLPK.Optimizer)
    @variable(model, 0 <= M[1:n, 1:n])
    @objective(model, Min, sum(sum(c(i, j)*M[i, j] for j in 1:n) for i in 1:n))
    @constraint(model, c1, M*ones(n) .== ones(n))
    @constraint(model, c2, M'*ones(n) .== ones(n))
    optimize!(model)
    best_perm_mat = value.(M)
    Plots.plot()
    for i = 1:n
        for j = 1:n
            if best_perm_mat[i, j] > 0
                Plots.plot!([x0[i], x1[j]], [y0[i], y1[j]], color = "black", lw=1.5)
            end
        end
    end
    Plots.scatter!(x0, y0, color = "blue", markersize = 5)
    display(Plots.scatter!(x1, y1, color = "red", size = (400, 400), markersize = 5,
        xlims = (0, 1), ylim = (0, 1), legend = false))
    savefig(filename)
end

matching_lp(11, "plots/matching.pdf")
