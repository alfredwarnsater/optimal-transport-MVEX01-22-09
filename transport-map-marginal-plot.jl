using Plots
using Plots.PlotMeasures

function transport_map_marginal_plot(marginal_1, marginal_2, transport_map)

    margin = -20;
    m1 = bar(transpose(marginal_1), orientation = :vertical, bar_width=1, yflip=false, bottom_margin=margin*px, 
        showaxis=false, ticks=false, legend = false,
        color = "red", linecolor= "red")
    m2 = bar(transpose(marginal_2), orientation = :horizontal, bar_width=1, xflip=true, right_margin=margin*px, 
        showaxis=false, ticks=false, legend = false,
        color = "blue", linecolor= "blue")
    hm = heatmap(transport_map, color = :greys, showaxis=false, ticks=false, legend = false, framestyle = :box)

    l = @layout[_ a; b c{0.7w, 0.7h}]

    plot(m1, m2, hm, size = (600, 600), layout = l, link = :both)

end