using CairoMakie

function plot_results(data, obstacle, filename)
    n_cols = 8
    n_rows = ceil(Int32, size(data, 1) / n_cols)
    replace!(x -> isapprox(x, 0) ? -1 : 0, obstacle)
    replace!(x -> (x > 0) ? 0 : x, obstacle)
    size_cm = 4 .* (n_cols, n_rows)
    size_pt = 28.3465 .* size_cm
    fig = Figure(resolution = size_pt, fontsize = 12)
    count = 1
    for row in 1:n_rows
        for col in 1:n_cols
            dist = reshape(data[count, :], isqrt(size(data[count, :], 1)), :)
            dist = dist .+ (count == 1 || count == n_rows * n_cols ? 0 : obstacle[count-1, :, :])
            ax, _ = heatmap(fig[row, col][1, 1], dist, colorrange=(0, maximum(data)), lowclip=:grey70, colormap=:jet)
            rowgap!(fig.layout, 10)
            colgap!(fig.layout, 10)
            hidedecorations!(ax)
            count = count + 1
        end
    end
    display(fig)
    save(filename, fig, pt_per_unit = 1)
end