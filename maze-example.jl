include("displacement-interpolation.jl")
include("plotting.jl")

function cost(x_1, y_1, x_2, y_2)
    if x_1 == x_2 && y_1 == y_2
        return 0
    elseif abs(x_1-x_2) + abs(y_1-y_2) == 1
        return 1
    else
        return Inf
    end
end

# Denna funktion definierar exemplets parametrar, startar lösaren och plottar resultatet.
function maze_example(L, epsilon, tol)
    maze = [  # Labyrinten som agenterna ska utrymma.
        1 1 0 0 0 0 0 0 0 0 0
        0 1 1 1 1 1 1 1 1 1 0 
        0 1 0 1 0 0 0 1 0 1 0 
        0 1 0 1 0 1 0 1 0 1 0 
        0 1 0 0 0 1 0 0 0 1 0 
        0 1 1 1 1 1 0 1 1 1 0 
        0 0 0 0 0 0 0 1 0 1 0 
        0 1 1 1 1 1 0 1 0 1 0 
        0 0 0 0 0 1 0 1 0 0 0 
        0 1 1 1 1 1 1 1 1 1 1
        0 0 0 0 0 0 0 0 0 0 1
    ]
    N = size(maze, 1)
    obstacle = zeros(L, N, N)
    for l in 1:L
        obstacle[l, :, :] = maze
    end
    pts = hcat(repeat(1:N, inner = N), repeat(1:N, outer = N))
    println("Genererar kostnadsmatris...")
    C = [cost(pts[i, 1], pts[i, 2], pts[j, 1], pts[j, 2]) for i in 1:(N*N), j in 1:(N*N)]
    mu = zeros(L, N*N)
    types = fill('<', L)
    mu[1, 1] = 1  # Agenterna börjar i ena hörnet.
    types[1] = '='
    for l = 2:(L-1)
        mu[l, :] = reshape(maze, 1, :)
    end
    mu[L, N*N] = 1  # Agenterna slutar i det motsatta hörnet.
    types[L] = '='
    println("Beräknar interpolation...")
    data = compute_interpolation(C, mu, types, epsilon, tol)
    println("Plottar...")
    plot_results(data, obstacle, "plots/maze-example.pdf", true)
    return
end

maze_example(40, 0.25, 0.01)