using LinearAlgebra

# Denna funktion löser interpolationsproblemet med hjälp av Sinkhorniterationer.
function compute_interpolation(C, mu, types, epsilon, tol)
    N = size(mu, 2)
    L = size(mu, 1)

    K = exp.(-C / epsilon)
    K_inv = inv(K)
    K_trans = transpose(K)
    u = ones(L, N)

    function gen_phi()
        phi = ones(L, N)
        phi[L-1, :] = K * u[L, :]
        for l in (L-2):(-1):1
            phi[l, :] =  K * Diagonal(u[l+1, :]) * phi[l+1, :]
        end
        return phi
    end

    count = 1
    max_diff = Inf
    # I denna loop används Sinkhorniterationerna för att beräkna lösningen.
    while max_diff > tol
        print("Iterationsrunda $(count): Maxdifferens: $(max_diff), " *
            "Vald tolerans: $(tol)\n")
        u_prev = copy(u)
        phi_hat = 1
        phi = gen_phi()
        for l in 1:L
            proj = u[l, :] .* phi_hat .* phi[l, :]
            tmp = mu[l, :] ./ proj[:]
            replace!(tmp, NaN => 1)
            u_tmp = u[l, :] .* tmp
            if types[l] == '<'  # Vi delar upp i fall beroende på typ av bivillkor.
                u[l, :] = min.(u_tmp, 1)
            elseif types[l] == '='
                u[l, :] = u_tmp
            elseif types[l] == '>'
                u[l, :] = max.(u_tmp, 1)
            end
            phi_hat = (l != 1 ? K_trans * Diagonal(u[l, :]) : K_trans * u[l, :]) * phi_hat
        end
        count = count + 1
        max_diff = maximum(abs.(u - u_prev))
    end
    
    phi_hat = 1
    phi = gen_phi()
    projs = zeros(L, N)
    for l in 1:L
        projs[l, :] = u[l, :] .* phi_hat .* phi[l, :]
        phi_hat = (l != 1 ? K_trans * Diagonal(u[l, :]) : K_trans * u[l, :]) * phi_hat
    end
    
    return projs
end



