using QuadGK
using ProgressMeter

include("displacement-interpolation.jl")
include("plotting.jl")

# Denna funktion beräknar kostnadsmatrisen som motsvarar dynamiken definierad av A och B.
function cost_matrix(grid_points, A, B)
    N = size(grid_points, 1)
    # Numerisk integrering.
    sigma, _ = quadgk(s -> exp(A*(1-s)) * B*transpose(B) * exp(transpose(A)*(1-s)), 0, 1)
    sigma_inv = inv(sigma)
    exp_A = exp(A)
    function cost(x_0, x_1)
        next!(p)
        tmp = x_1-exp_A*x_0
        return transpose(tmp)*sigma_inv*tmp
    end
    # Här beräknas kostnaderna mellan alla par av möjliga tillstånd.
    p = Progress(N^2,
                 dt=0.5,
                 desc="Genererar kostnadsmatrisen",
                 barglyphs=BarGlyphs("[=> ]"),
                 barlen=50)
    C = [cost(grid_points[i, :], grid_points[j, :]) for i in 1:N, j in 1:N]
    return C
end

# Denna funktion beräknar rutnätet utifrån de givna punkterna.
function gen_grid_points(points)
    N = length(points)-1
    grid = points[1:N] .+ 0.5/N  # Vi samplar i mitten av varje ruta.
    grid_points = hcat(repeat(grid, inner = N), repeat(grid, outer = N))
    return grid_points
end

# Denna funktion genererar hindret för agenterna.
function gen_obstacle(N, n_steps)

    # Parametriserad ellips.
    r = 0.2
    x_curve(t) = 2*r*cos.(t)
    y_curve(t) = r*sin.(t)
    x_start, x_end = 0, 1
    y_start, y_end = 0, 1

    h = 1 / N
    n_a_pts = 10 * N
    n_r_pts = 10 * N

    inside(x) = all((1 .<= x .<= N))
    rotate(p, theta) = [p[1]*cos(theta)-p[2]*sin(theta), p[1]*sin(theta)+p[2]*cos(theta)]
    snap(x) = round.(Int32, x / h) .+ 1

    function transform(pts, t_vec, r_angle)
        t_pts = zeros(size(pts))
        for i in 1:size(pts, 1)
            t_pts[i, :] = rotate(pts[i, :], r_angle) + t_vec[:]
        end
        return t_pts
    end

    function gen_obs(pts)
        obs = zeros(N, N)
        pts = snap(pts)
        for i in 1:size(pts, 1)
            if inside(pts[i, :])
                obs[pts[i, :]...] = 1
            end
        end
        return obs
    end

    a_pts = range(0, 2*pi, n_a_pts)
    r_pts = range(0, 1, n_r_pts)
    obs_pts = zeros(n_r_pts*n_a_pts, 2)
    for i in 1:n_r_pts
        obs_pts[(i-1).*n_a_pts .+ (1:n_a_pts), :] = r_pts[i] .*[x_curve(a_pts) y_curve(a_pts)]
    end

    xs = range(x_start, x_end, n_steps)
    ys = range(y_start, y_end, n_steps)
    rs = range(0, 2*pi, n_steps)
    obstacle = zeros(n_steps, N, N)

    for l in range(1, n_steps)
        t_pts = transform(obs_pts, [xs[l], ys[l]], rs[l])
        obstacle[l, :, :] = gen_obs(t_pts)
    end

    return obstacle
end

# Denna funktion definierar exemplets parametrar, startar lösaren och plottar resultatet.
function dynamics_example(N, L, epsilon, tol)
    # Definition av matriserna som definierar det dynamiska systemet.
    A = [
        0 0
        0 0
    ]
    B = [
        1 0
        0 1
    ]
    # Parametrar för startfördelningen.
    m = [0.5, 0.5]
    s = [0.2, 0.2]
    # Exempelfördelning. Baserad på täthetsfunktionen för en tvådimensionell normalfördelning.
    f(x, m, s) = exp(-((x[1]-m[1])^2/(2*s[1]^2) + (x[2]-m[2])^2/(2*s[2]^2)))
    pts = range(0, 1, N+1)  # Ekvidistanta punkter som definierar diskretiseringen.
    grid_pts = gen_grid_points(pts)
    println("Genererar kostnadsmatris...")
    C = cost_matrix(grid_pts, A, B)
    mu = zeros(L, N*N)
    types = fill('-', L)
    mu[1, :] = reshape([f(row, m, s) for row in eachrow(grid_pts)], 1, :)  # Startfördelningen.
    types[1] = '='
    total_mass = sum(mu[1, :])
    obstacle = gen_obstacle(N, L-2)  # Vi lägger till ett hinder för agenterna.
    obstacle = (obstacle .- 1) * (-total_mass)
    for l in 2:(L-1)  # Hindret modelleras med olikhetsbivillkor.
        mu[l, :] = reshape(obstacle[l-1, :, :], 1, :)
        types[l] = '<'
    end
    mu[L, :] = reshape(total_mass * ones(N, N) / N^2, 1, :)  # Slutfördelningen.
    types[L] = '='
    println("Beräknar interpolation...")
    data = compute_interpolation(C, mu, types, epsilon, tol)
    println("Plottar...")
    plot_results(data, obstacle, "plots\\dynamics-example.pdf", false)
    return
end

dynamics_example(10, 8*8, 0.01, 0.01)
