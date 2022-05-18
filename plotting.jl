using CairoMakie

# Denna funktion plottar agenternas fördelning och eventuella hinder.
function plot_results(data, obstacle, filename, is_maze)
    L = size(data, 1)
    n_cols = 8
    n_rows = ceil(Int32, L / n_cols)
    replace!(x -> isapprox(x, 0) ? -1 : 0, obstacle)
    replace!(x -> (x > 0) ? 0 : x, obstacle)
    size_cm = 4 .* (n_cols, n_rows)
    size_pt = 28.3465 .* size_cm
    global fig = Figure(resolution = size_pt, fontsize = 12)
    count = 1
    done = false
    for row in 1:n_rows
        for col in 1:n_cols
            if count > L
                done = true
                break
            end
            dist = reshape(data[count, :], isqrt(size(data[count, :], 1)), :)
            if !is_maze
                dist = dist .+ (count == 1 || count == L ? 0 : obstacle[count-1, :, :])
                ax, _ = CairoMakie.heatmap(fig[row, col][1, 1], dist, colormap=:jet,
                    colorrange=(0, maximum(data)), lowclip=:grey80, highclip=:red)
            else
                dist = dist .+ obstacle[count, :, :]
                ax, _ = CairoMakie.heatmap(fig[row, col][1, 1], dist, colormap=:jet,
                    colorrange=(0, 0.9), lowclip=:grey80, highclip=:red)
            end
            rowgap!(fig.layout, 10)
            colgap!(fig.layout, 10)
            hidedecorations!(ax)
            count = count + 1
        end
        if done
            break
        end
    end
    display(fig)
    #Hämtar arbetskatalogen
    WD = dirname(@__FILE__)
    save(string(WD,"\\",filename), fig, pt_per_unit = 1)
    #Tar bort .pdf i namnet på filen
    foldername = replace(filename,".pdf" => "")
    pathToFolder = string(WD,"\\",foldername,"\\")
    for i in 1:(n_rows*n_cols)
        resize!(fig.scene.children[i],(100,100));
        absolutePathPNG = pathToFolder*string(i)*".png"
        absolutePathPDF = pathToFolder*string(i)*".pdf"
        save(absolutePathPNG, fig.scene.children[i], pt_per_unit = 1)
        save(absolutePathPDF, fig.scene.children[i], pt_per_unit = 1)
    end
end
