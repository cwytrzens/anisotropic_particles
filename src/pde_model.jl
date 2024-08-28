
function laplace_nonlinear!(du, u, G, dV)
    
    # period BC
    get_u(i,j) = u[ mod(i, axes(u,1)), mod(j, axes(u,2)) ]

    @inbounds for i in axes(u,1), j in axes(u,2)

        uc = get_u(i,j)
        ul = get_u(i-1,j)
        ur = get_u(i+1,j)
        ub = get_u(i,j-1)
        ut = get_u(i,j+1)

        Gl = G(0.5*ul + 0.5*uc)
        Gr = G(0.5*ur + 0.5*uc)
        Gb = G(0.5*ub + 0.5*uc)
        Gt = G(0.5*ut + 0.5*uc)

        du[i,j] += (
            Gl * (ul - uc) / dV[1]^2 + 
            Gr * (ur - uc) / dV[1]^2 +
            Gb * (ub - uc) / dV[2]^2 +
            Gt * (ut - uc) / dV[2]^2
            )
    end
    return nothing
end




function init(p)
    (;nx, ny, Lx, Ly) = p

    xs = LinRange(0, p.Lx, p.nx)
    ys = LinRange(0, p.Ly, p.ny)
    sigma = 2.0
    u0 = @. 30 + exp(-(xs - 2*p.Lx/3)^2/sigma^2 - (ys' - p.Ly/2)^2/sigma^2)  + exp(-(xs - p.Lx/3)^2/sigma^2 - (ys' - p.Ly/2)^2/sigma^2) # .* ( 1.1 .+ rand(nx, ny) )
end



function rhs!(du, u, p, t)
    (; D_x, mu, cache) = p
    (; K_rho) = cache

    rho_domain = (cache.data.rhos[1], cache.data.rhos[end])

    du .= 0.0

    dV = ( p.Lx / p.nx, p.Ly / p.ny)
    G(rho) = if rho_domain[1] < rho < rho_domain[2]
        D_x + mu * K_rho(rho) * rho
    else        
        0.0
    end

    laplace_nonlinear!(du, u, G, dV)

    return nothing
end

