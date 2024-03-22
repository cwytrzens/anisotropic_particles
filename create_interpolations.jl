using DataInterpolations, Integrals
include("definitions.jl")

# define functions to compute S₂:

function Z_g(eta)
    g(theta, eta) = exp(eta * cos(theta)^2)

    prob = IntegralProblem(g, (0.0, 1.0*π), eta)
    return solve(prob, QuadGKJL()).u
end

g(eta, theta) = exp(eta * cos(theta)^2) / Z_g(eta)

function S2(eta, n) 
    f(theta, n) = n/(n-1) * (cos(theta)^2 - 1/n) * g(eta, theta) * sin(theta)^(n-2)

    prob = IntegralProblem(f, (0.0, 1.0*π), n)
    return solve(prob, QuadGKJL()).u
end


# S2(0.1, 10)

function create_interpolations(p; eta_range = 1000.0 * (1 - p.chi^2), n_interp = 2000)

    # interpolate eta(rho)
    alpha = 1 - p.chi^2

    etas = LinRange(0, eta_range, n_interp)[2:end]
    S2s = S2.(etas, 2)
    rhos = @. etas / (S2s * alpha)

    printstyled("Domain of ρ values: $(minimum(rhos)) to $(maximum(rhos))", color = :red)


    # create interpolation object 
    eta_rho = CubicSpline(etas, rhos)
    S2_eta = CubicSpline(S2s, etas)


    function K(rho, p, cache) 
        (; chi, sigma) = p
        (; eta_rho, S2_eta) = cache
        n = 2

        eta = eta_rho(rho)
        eta_prime = DataInterpolations.derivative(eta_rho, rho)

        return 1 - chi^2 / n - sigma * (n-1)/n * S2_eta(eta) * eta_prime 
    end

    Ks = [K(rho, p, (;eta_rho, S2_eta)) for rho in rhos]

    K_rho = CubicSpline(Ks, rhos) 

    return (; 
        eta_rho, S2_eta, K_rho, 
        data = (;
            etas, rhos, S2s, Ks    
        )
    )
end
