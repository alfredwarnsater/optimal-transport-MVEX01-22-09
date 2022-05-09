using QuadGK
using JuMP
using GLPK
using Plots
using ControlSystems
using Gurobi

c(x0, x1, A, B, sigma_inv) = transpose(x1-exp(A)*x0)*sigma_inv*(x1-exp(A)*x0)

function get_C(sigma_inv, ls, A, B, n)
    C = zeros(n, n, n, n)
    for i in 1:n
        for j in 1:n
            for k in 1:n
                for l in 1:n
                    C[i,j,k,l] = c([ls[i], ls[j]], [ls[k], ls[l]], A, B, sigma_inv)
                end
            end
        end
    end
    return C
end

function get_sigma_inv(A, B)
    I, _ = quadgk(s -> exp(A*(1-s))*B*transpose(B)*exp(transpose(A)*(1-s)), 0, 1)
    I = inv(I);
    return I
end

function get_middle_dist(fd, sd, A, B, n)

    @time begin
        ls = 0:1/(n-1):1
        sigma_inv = get_sigma_inv(A, B)
        C = get_C(sigma_inv, ls, A, B, n)
    end
    
    @time begin
        model = Model(Gurobi.Optimizer)
        @variable(model, 0 <= M1[1:n, 1:n, 1:n, 1:n])
        @variable(model, 0 <= M2[1:n, 1:n, 1:n, 1:n])
        @variable(model, 0 <= md[1:n, 1:n])
    end

    @time begin
        @objective(model, Min,
            sum(C .* (M1 +  M2))
            #sum(sum(sum(sum(C[i,j,k,l]*M1[i, j, k, l] for i in 1:n) for j in 1:n) for k in 1:n) for l in 1:n) +
            #sum(sum(sum(sum(C[i,j,k,l]*M2[i, j, k, l] for i in 1:n) for j in 1:n) for k in 1:n) for l in 1:n)
    )
    end

    @time begin
        @constraint(model, con1[i=1:n, j=1:n], sum(M1[i,j,k,l] for k=1:n, l=1:n) == fd[i, j])
        @constraint(model, con2[k=1:n, l=1:n], sum(M1[i,j,k,l] for i=1:n, j=1:n) == md[k, l])
        @constraint(model, con3[i=1:n, j=1:n], sum(M2[i,j,k,l] for k=1:n, l=1:n) == md[i, j])
        @constraint(model, con4[k=1:n, l=1:n], sum(M2[i,j,k,l] for i=1:n, j=1:n) == sd[k, l])
    end
    
    @time begin
        #set_optimizer_attribute(model, "msg_lev", GLPK.GLP_MSG_ALL)
        optimize!(model)
    end

    return value.(md)
end

function main()
    A = [
        0 0
        0 0
    ]
    B = [
        1 0
        0 1
    ]

    n = 50
    total_mass = 1000

    n_0 = 5
    mu0 = zeros(n, n)
    mu0[1:n_0, 1:n_0] .= total_mass / n_0^2

    mu2 = zeros(n, n)
    mu2[n, 1:n] .= total_mass / (2*n-1)
    mu2[1:(n-1), n] .= total_mass / (2*n-1)

#=     mu0 = zeros(n, n)
    mu0[1:n, 1] .= total_mass / n
    mu2 = zeros(n, n)
    mu2[1:n, n] .= total_mass / n =#
    
    mu1 = get_middle_dist(mu0, mu2, A, B, n);

    return mu1
end

data = main();