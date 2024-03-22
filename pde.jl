using DataInterpolations, Integrals, OrdinaryDiffEq
using ProgressLogging
using ForwardDiff
using GLMakie
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


p = loadparameters("params.toml")

# interpolate eta(rho)

alpha = 1 - p.chi^2
eta_range = 1000.0 * alpha
n_interp = 2000

etas = LinRange(0, eta_range, n_interp)[2:end]
S2s = S2.(etas, 2)
rhos = @. etas / (S2s * alpha)

printstyled("Domain of ρ values: $(minimum(rhos)) to $(maximum(rhos))", color = :red)


# create interpolation object 
eta_rho = CubicSpline(etas, rhos)
S2_eta = CubicSpline(S2s, etas)

begin
    fig = Figure()

    inds = etas .< 40.0
    
    Axis(fig[1,1], title = "ρ(η)", xlabel = "η", ylabel = "ρ")
    lines!(etas[inds], rhos[inds])

    Axis(fig[1,2], title = "S₂(η)", xlabel = "η", ylabel = "S₂")
    lines!(etas[inds], S2s[inds])

    fig
end


# usage:
# eval:  eta_rho(1.0)
# deriv: DataInterpolations.derivative(eta_rho, 40.0)


function K(rho, p, cache) 
    (; chi, sigma) = p
    (; eta_rho, S2_eta) = cache
    n = 2

    eta = eta_rho(rho)
    eta_prime = DataInterpolations.derivative(eta_rho, rho)

    return 1 - chi^2 / n - sigma * (n-1)/n * S2_eta(eta) * eta_prime 
end

 
# interpolate from rho to K 
begin 
    cache = (;eta_rho, S2_eta)

    fig = Figure()
    ax = Axis(fig[1,1], title = "ρ ↦ K(η(ρ))", xlabel = "ρ", ylabel = "K")

    xs = LinRange(60.0, 1000.0, 400)
    Ks = [K(rho, p, cache) for rho in xs]

    #scatter!(xs, Ks .- minimum(Ks))
    scatter!(xs, Ks)  # rounding to Float32 seems to be a problem here
    ax.ylabel[] = "K, (shifted by -$(round(minimum(Ks), digits=2)))"

    fig
end



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











const Dim = 2

function init(p)
    (;nx, ny, Lx, Ly) = p

    xs = LinRange(0, p.Lx, p.nx)
    ys = LinRange(0, p.Ly, p.ny)

    u0 = minimum(rhos) .* ( 1.1 .+ rand(nx, ny) )
end

u0 = init(p)
tspan = (p.t_start, p.t_end)

function rhs!(du, u, p, t)
    (; D_x, mu, cache) = p

    rho_domain = (cache.eta_rho.t[1], cache.eta_rho.t[end])

    du .= 0.0

    dV = ( p.Lx / p.nx, p.Ly / p.ny)
    G(rho) = if rho_domain[1] < rho < rho_domain[2]
        D_x + mu * K(rho, p, cache) * rho
    else        
        0.0
    end

    laplace_nonlinear!(du, u, G, dV)

    return nothing
end

p_ode = (;p..., cache)
prob = ODEProblem(rhs!, u0, tspan, p_ode)

sol = solve(prob, ROCK2(), progress = true)


begin
    fig = Figure() 
    ax = Axis(fig[1,1])
    sl = Slider(fig[2,1], range = eachindex(sol))
    sol_obs = @lift sol[$(sl.value)]

    heatmap!(sol_obs, colorrange = (minimum(sol), maximum(sol)))

   fig
end

# u = u0
# du = similar(u0)
# @time rhs!(du, u, p_ode, 0.1)