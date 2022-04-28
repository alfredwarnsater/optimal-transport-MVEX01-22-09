using CairoMakie

function plot_results(data, maze)
    replace!(maze, 0 => -1)
    replace!(maze, 1 => 0)
    fig = Figure()
    count = 1
    for row in 1:5
        for col in 1:5
            dist = reshape(data[count], isqrt(size(data[count], 1)), :)
            dist = dist + maze
            ax, hm = heatmap(fig[row, col][1, 1], dist, colorrange = (0, 1), lowclip = :blue)
            #Colorbar(fig[1, 1][1, 2], hm)
            hidedecorations!(ax)
            count = count + 1
        end
    end


#=     ax, hm = heatmap(fig[1, 2][1, 1], xs, ys, zs, colormap = :grays,
        colorrange = (0, 1), highclip = :red, lowclip = :blue)
    Colorbar(fig[1, 2][1, 2], hm)
    
    ax, hm = contourf(fig[2, 1][1, 1], xs, ys, zs,
        levels = -1:0.25:1, colormap = :heat)
    Colorbar(fig[2, 1][1, 2], hm, ticks = -1:0.25:1)
    
    ax, hm = contourf(fig[2, 2][1, 1], xs, ys, zs,
        colormap = :Spectral, levels = [-1, -0.5, -0.25, 0, 0.25, 0.5, 1])
    Colorbar(fig[2, 2][1, 2], hm, ticks = -1:0.25:1) =#
    
    fig

        #plot_array = []
#=     for i in 1:size(data, 1)
        dist = reshape(data[i], isqrt(size(data[i], 1)), :)
        push!(plot_array, heatmap(dist, showaxis=false, ticks=false, legend=false, framestyle = :box,
            clims = (0, maximum(maximum.(data))), c = :jet))
    end
    plot(plot_array..., size = (1000, 1000), link = :both) =#
end

#plot_results([1])