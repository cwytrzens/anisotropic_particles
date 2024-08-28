
# define functions to compute S₂:
function Z_g(eta)
    prob = IntegralProblem((θ, η) -> exp(η*cos(θ)^2), (0.0, 1.0*π), eta)
    return solve(prob, QuadGKJL()).u
end

g(eta, theta) = exp(eta * cos(theta)^2) / Z_g(eta)

function S2(eta, n) 
    function f(θ, p) 
        (n, η) = p 
        return n/(n-1) * (cos(θ)^2 - 1/n) * g(η, θ) * sin(θ)^(n-2)
    end

    prob = IntegralProblem(f, (0.0, 1.0*π), (n, eta))
    return max(eps(1.0), solve(prob, QuadGKJL()).u)
end


function create_interpolations(p; kwargs...)
    return create_interpolations(p.chi, p.D_u, p.lambda; kwargs...)
end 

function create_interpolations(chi, D_u, lambda; deg = 1.0, eta_range = 100.0 * (1 - chi^2), n_interp = 1000)
    
    alpha = chi^2 * lambda / D_u
    sigma = D_u / lambda

    # interpolate eta(rho)
    etas = LinRange(0.0, 1.0, n_interp) .^ deg .* eta_range
    S2s = S2.(etas, 2)
    rhos = @. etas / (S2s * alpha)

    printstyled("Domain of ρ values: $(minimum(rhos)) to $(maximum(rhos))\n", color = :red)


    # create interpolation object 
    eta_rho = AkimaInterpolation(etas, rhos)  # η(ρ)
    S2_eta = AkimaInterpolation(S2s, etas)    # S2(η)


    function K(rho, p, cache) 
        (; chi, sigma) = p

        (; eta_rho, S2_eta) = cache
        n = 2

        eta = eta_rho(rho)  # evaluate η(ρ)
        eta_prime = DataInterpolations.derivative(eta_rho, rho)  # η'(ρ)

        return 1 - chi^2 / n - sigma * (n-1)/n * S2_eta(eta) * eta_prime 
    end

    Ks = [K(rho, (; chi, sigma), (; eta_rho, S2_eta)) for rho in rhos]


    K_eta = AkimaInterpolation(Ks, etas) 
    K_rho = AkimaInterpolation(Ks, rhos) 

    return (; 
        eta_rho, S2_eta, K_rho, K_eta,
        data = (;etas, rhos, S2s, Ks)
    )
end


function doit()
    p = loadparameters("paper/inputs/base.toml")
    # p = (; chi = 0.996, D_u = 0.48, lambda = 0.1)
    intp = create_interpolations(p; eta_range = 40.0, n_interp = 2000, deg = 1.0)
    Makie.inline!(true)
    begin 
        fig = Figure()
        ax = Axis(fig[1,1], xlabel = "η", ylabel = "K(η(ρ))", title = "Positivity of mass density")

        lines!(intp.data.etas, intp.data.Ks)

        fig 
    end 
end
