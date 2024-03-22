using DataInterpolations, Integrals, OrdinaryDiffEq
using ProgressLogging
using ForwardDiff
using GLMakie

# define functions to compute S₂:

function Z_g(eta)
    g(theta, eta) = exp(eta * cos(theta)^2)

    prob = IntegralProblem(g, (0.0, 1.0*π), eta)
    return solve(prob, QuadGKJL()).u
end

g(eta, theta) = exp(eta * cos(theta)^2) / Z_g(eta)

function S_2(eta, n) 
    f(theta, n) = n/(n-1) * (cos(theta)^2 - 1/n) * g(eta, theta) * sin(theta)^(n-2)

    prob = IntegralProblem(f, (0.0, 1.0*π), n)
    return solve(prob, QuadGKJL()).u
end


# S_2(0.1, 10)


# interpolate eta(rho)

alpha = 1 - p.chi^2
eta_range = 200.0
n_interp = 2000

etas = LinRange(0, eta_range, n_interp)[2:end]
S2s = S_2.(etas, 2)
rhos = @. etas / (S2s * alpha)

printstyled("Minimal value of ρ: $(minimum(rhos))", color = :red)


# create interpolation object 
eta_rho = CubicSpline(etas, rhos)
S2_eta = CubicSpline(S2s, etas)

begin
    fig = Figure()
    
    Axis(fig[1,1], title = "ρ(η)", xlabel = "η", ylabel = "ρ")
    lines!(etas, rhos)
    xlims!(0, 40)

    Axis(fig[1,2], title = "S₂(η)", xlabel = "η", ylabel = "S₂")
    lines!(etas, S2s)
    xlims!(0, 40)

    fig
end


# usage:
# eval:  eta_rho(1.0)
# deriv: DataInterpolations.derivative(eta_rho, 40.0)


p = loadparameters("params.toml")

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
    ax = Axis(fig[1,1])

    xs = LinRange(60.0, 1000.0, 400)
    Ks = [K(rho, p, cache) for rho in xs]

    scatter!(xs, Ks .- minimum(Ks))
    #scatter!(xs, Ks)


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

        Gc = G(uc)
        Gl = 0.5*G(ul) + 0.5*Gc
        Gr = 0.5*G(ur) + 0.5*Gc
        Gb = 0.5*G(ub) + 0.5*Gc
        Gt = 0.5*G(ut) + 0.5*Gc


        du[i,j] += (
            Gl * (ul - uc) / dV[1]^2 + 
            Gr * (ur - uc) / dV[1]^2 +
            Gb * (ub - uc) / dV[2]^2 +
            Gt * (ut - uc) / dV[2]^2
            )
    end
    return nothing
end












p = loadparameters("params.toml")

const Dim = 2

function init(p)
    (;nx, ny, Lx, Ly) = p

    xs = LinRange(0, p.Lx, p.nx)
    ys = LinRange(0, p.Ly, p.ny)

    u0 = 80.0 .+ rand(nx, ny)
end

u0 = init(p)
tspan = (p.t_start, p.t_end)

function rhs!(du, u, p, t)
    (; D_x, mu, cache) = p

    du .= 0.0

    dV = ( p.Lx / p.nx, p.Ly / p.ny)
    G(rho) = D_x + mu * K(rho, p, cache) * rho

    laplace_nonlinear!(du, u, G, dV)

    return nothing
end

p_ode = (;p..., cache)
prob = ODEProblem(rhs!, u0, tspan, p_ode)

sol = solve(prob, Euler(), dt = 0.1, progress = true)

u = u0
du = similar(u0)
rhs!(du, u, p_ode, 0.1)