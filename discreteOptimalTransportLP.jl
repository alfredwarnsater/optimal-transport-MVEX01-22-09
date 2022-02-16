using JuMP
using GLPK
using Plots

function mainModel(n)
    global x0 = rand(Float64, n)
    global y0 = rand(Float64, n)
    global x1 = rand(Float64, n)
    global y1 = rand(Float64, n)
    global mu_0 = zeros(n)
    global mu_1 = zeros(n)

    mu_0 = [0,0,0,2,3]  #Blue nodes
    mu_1 = [1,1,1,1,1]  #Red nodes

    if ones(n)'*mu_0 != ones(n)'*mu_1
         throw(DomainError("mu_0 is not euqal to mu_1"))
    end


    c(i, j) = (x0[i]-x1[j])^2 + (y0[i]-y1[j])^2

    model = Model(GLPK.Optimizer)
    @variable(model, M[1:n, 1:n] >= 0)
    @objective(model, Min, sum(sum(c(i, j)*M[i, j] for j in 1:n) for i in 1:n))
    @constraint(model, M*ones(n) .== mu_0)
    @constraint(model, M'*ones(n) .== mu_1)
    optimize!(model)
    return value.(M)
end

function mainPlot(best_perm_mat)
    scatter(x0, y0, color = "blue")
    for i = 1:agents
        for j = 1:agents
            #If there is a way between node i and j, draw a black line.
            if best_perm_mat[i, j] > 0
                plot!([x0[i], x1[j]], [y0[i], y1[j]], color = "black")
            end
        end
    end
    display(scatter!(x1, y1, color = "red", size = (400, 400), xlims = (0, 1), ylim = (0, 1), legend = false))
end

#Here is where the code is executed!
println("\n\nv-v-v-v-v-v-v-v-v-v-v-v-v-v-v-v    New run    v-v-v-v-v-v-v-v-v-v-v-v-v-v-v-v-v-v-v-v")
agents = 5 #Number of agents in the field...
@time best_perm_mat = mainModel(agents)
@time mainPlot(best_perm_mat)

if (agents < 15)
    for i = 1:agents
        print("  ",i,"  ")
    end
    print("\n")
    for i = 1:agents
        println(best_perm_mat[i,:]," = ", sum(best_perm_mat[i,:]))
    end
    for i = 1:agents
        print("  || ")
    end
    print("\n")
    for i = 1:agents
        print(" ",sum(best_perm_mat[:,i])," ")
    end
    print("\n")
end
println("\nMu_0 = ",mu_0)
println("\nMu_1 = ",mu_1)
