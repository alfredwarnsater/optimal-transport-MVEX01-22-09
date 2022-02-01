using JuMP
using GLPK
using Plots

function main(n)
    x0 = rand(Float64, n)
    y0 = rand(Float64, n)
    x1 = rand(Float64, n)
    y1 = rand(Float64, n)
    c(i, j) = (x0[i]-x1[j])^2 + (y0[i]-y1[j])^2

    model = Model(GLPK.Optimizer)
    @variable(model, 0 <= M[1:n, 1:n] <= 1)
    @objective(model, Min, sum(sum(c(i, j)*M[i, j] for j in 1:n) for i in 1:n))
    @constraint(model, c1, M*ones(n) .== ones(n))
    @constraint(model, c2, M'*ones(n) .== ones(n))
    optimize!(model)
    best_perm_mat = value.(M)

    scatter(x0, y0, color = "blue")
    for i = 1:n
        for j = 1:n
            if best_perm_mat[i, j] > 0
                plot!([x0[i], x1[j]], [y0[i], y1[j]], color = "black")
            end
        end
    end
    display(scatter!(x1, y1, color = "red", size = (400, 400), xlims = (0, 1), ylim = (0, 1), legend = false))
end

main(30)