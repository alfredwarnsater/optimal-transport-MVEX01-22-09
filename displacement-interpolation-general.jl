function compute_interpolation(C, mu, types, epsilon, tol)
    L = size(mu, 1)
    N = size(mu, 2)
    
    K = exp.(-C / epsilon)
    K_inv = inv(K)
    K_trans = transpose(K)
    u = ones(L, N)

    function P(l)
        function phi_hat(l)
            ans = 1
            for i in reverse(1:(l-1))
                ans = ans * transpose(K) * (i != 1 ? Diagonal(u[i, :]) : u[i, :]) 
            end
            return ans
        end

        function phi(l)
            ans = 1
            for i in (l+1):L
                ans = ans * K * (i != L ? Diagonal(u[i, :]) : u[i, :])
            end
            return ans
        end

        return u[l, :] .* phi_hat(l) .* phi(l)
    end

    function gen_phi()
        phi = ones(L, N)
        phi[L-1, :] = K * u[L, :]
        for l in (L-2):(-1):1
            phi[l, :] =  K * Diagonal(u[l+1, :]) * phi[l+1, :]
        end
        return phi
    end

    u_prev = Inf * ones(size(u))
    count = 1
    last_proj = zeros(L, N)
    while any(abs.(u - u_prev) .> tol)
        u_prev = copy(u)
        phi_hat = 1
        phi = gen_phi()
        for l in 1:L
            proj = u[l, :] .* phi_hat .* phi[l, :]
            last_proj[l, :] = proj
            tmp = mu[l, :] ./ proj[:]
            replace!(tmp, NaN => 1)
            it = u[l, :] .* tmp
            if types[l] == '<'
                u[l, :] = min.(it, 1)
            elseif types[l] == '='
                u[l, :] = it
            elseif types[l] == '>'
                u[l, :] = max.(it, 1)
            end
            phi_hat = (l != 1 ? K_trans * Diagonal(u[l, :]) : K_trans * u[l, :]) * phi_hat
        end
        count = count + 1
    end

    print("count: ")
    println(count-1)
    return [last_proj[i, :] for i in 1:L]
end



