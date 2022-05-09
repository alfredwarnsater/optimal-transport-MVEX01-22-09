using Plots
using Permutations

function main(n)
    x0 = rand(Float64, n)
    y0 = rand(Float64, n)
    x1 = rand(Float64, n)
    y1 = rand(Float64, n)

    best_perm = Permutation(n)
    best_sum = Base.Inf64

    for p in PermGen(n)
        sum = calculate_distance_sum(x0, y0, x1, y1, p)
        if sum < best_sum
            best_sum = sum
            best_perm = p
        end
    end

    scatter(x0, y0, color = "blue")
    for i = 1:length(best_perm)
        j = best_perm[i]
        plot!([x0[i], x1[j]], [y0[i], y1[j]], color = "black")
    end
    display(scatter!(x1, y1, color = "red", size = (400, 400), xlims = (0, 1), ylim = (0, 1), legend = false))
end

function calculate_distance_sum(x0, y0, x1, y1, p)
    s = 0
    for i = 1:length(p)
        s = s + (x0[i]-x1[p[i]])^2 + (y0[i]-y1[p[i]])^2;
    end
    return s
end

@time main(11)
