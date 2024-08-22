using DataInterpolations, Integrals
include("definitions.jl")

# define functions to compute S₂:

function Z_g(eta)
    g_(theta, eta) = exp(eta * cos(theta)^2)

    prob = IntegralProblem(g_, (0.0, 1.0*π), eta)
    return solve(prob, QuadGKJL()).u
end

g(eta, theta) = exp(eta * cos(theta)^2) / Z_g(eta)

function S2(eta, n) 
    function f(theta, p) 
        (n, eta) = p 
        return n/(n-1) * (cos(theta)^2 - 1/n) * g(eta, theta) * sin(theta)^(n-2)
    end

    prob = IntegralProblem(f, (0.0, 1.0*π), (n, eta))
    return max(eps(1.0), solve(prob, QuadGKJL()).u)
end


# S2(0.1, 10)
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
    eta_rho = AkimaInterpolation(etas, rhos)
    S2_eta = AkimaInterpolation(S2s, etas)


    function K(rho, p, cache) 
        (; chi, sigma) = p

        (; eta_rho, S2_eta) = cache
        n = 2

        eta = eta_rho(rho)
        eta_prime = DataInterpolations.derivative(eta_rho, rho)

        return 1 - chi^2 / n - sigma * (n-1)/n * S2_eta(eta) * eta_prime 
    end

    Ks = [K(rho, (; chi, sigma), (; eta_rho, S2_eta)) for rho in rhos[1:end]]

    K_rho = AkimaInterpolation(Ks, rhos[1:end]) 

    return (; 
        eta_rho, S2_eta, K_rho, 
        data = (;
            etas = etas[1:end], rhos = rhos[1:end], S2s, Ks    
        ),
        min_K = minimum(Ks)
    )
end

using ProgressMeter

samples = []

run = Ref(true)
run[] = false


@async while run[] 

    # rescale lambda around 150 
    params = (;chi = rand(), D_u = tan(pi/2 * rand()), lambda = tan(pi/2 * rand()))
    min_K = create_interpolations( Tuple(params)...).min_K 

    println(length(samples), ": ", params, "    min_K: ", min_K)
    push!(samples, (; params..., min_K))

    if min_K < 0.0 
        run[] = false 
    end
end

create_interpolations(0.996, 0.48, 0.1; eta_range = 5.0, n_interp = 2000).min_K


begin 
 fig = Figure() 
 ax = Axis(fig[1,1])
k = 8
 data = [s for s in samples if k*0.1 < s.chi < (k+1)*0.1]
x = [s.D_u for s in data]
y = [s.lambda for s in data]
c = [s.min_K for s in data]
    scatter!(x,y,color = c)
    xlims!(ax, 0, 5)
    ylims!(ax, 0, 5)
 fig
end

scatter(Chi, Kmin)

function doit()
    p = loadparameters("inputs/mass_density/mass_density_example.toml")
    p = (; chi = 0.996, D_u = 0.48, lambda = 0.1)
    intp = create_interpolations(p; eta_range = 40.0, n_interp = 2000, deg = 1.0)
    Makie.inline!(true)
    begin 
        fig = Figure()
        ax = Axis(fig[1,1], xlabel = "ρ", ylabel = "K(η(ρ))", title = "Positivity of mass density")

        lines!(intp.data.etas, intp.data.Ks)

        fig 
    end 
end

doit()






begin 
    n_interp = 5000
    eta_range = 1.0

    chi = 0.996
    D_u = 0.48
    lambda = 0.1
    
    alpha = chi^2 * lambda / D_u
    sigma = D_u / lambda

    # interpolate eta(rho)
    etas = LinRange(0.0, eta_range, n_interp)[2:end]
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

    Ks = [K(rho, (; chi, sigma), (; eta_rho, S2_eta)) for rho in rhos]

end