using Interpolations, Integrals, OrdinaryDiffEq
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

alpha = 0.5
eta_range = 100.0
n_interp = 1000

etas = LinRange(0, eta_range, n_interp)
S2s = S_2.(etas, 2)
rhos = @. etas / (S2s * alpha)

begin
    fig = Figure()
    Axis(fig[1,1])

    lines(etas, rhos)
    

end

K(u) = 1.0 
nonlinear_term(u) = K(u) * u



struct PeriodicBoundary{T}
    data::T
end

struct NeumannBoundary{T}
    data::T
end

struct DirichletBoundary{T}
    data::T
end


bcindex(l, ax, bnd::NeumannBoundary) = clamp(l, ax)
bcindex(l, ax, bnd::PeriodicBoundary) = mod(l, ax)
bcindex(l, ax, bnd::DirichletBoundary) = l in ax ? l : 0

bc(u::Float64, inds) = u
function bc(u::AbstractArray, inds)    
    I = map(ind -> bcindex(ind...), inds)
    i_bnd = findfirst(I .== 0)
    return isnothing(i_bnd) ? u[I...] : inds[i_bnd][3].data::Float64
end

function laplace!(du::AbstractArray{Float64,2}, u::AbstractArray{Float64,2}, D, dV; factor = 1.0, boundaries = (NeumannBoundary(nothing), NeumannBoundary(nothing)))
    bnds = boundaries
    
    get_D(i,j) = D * bc(factor, ((i,axes(factor,1), NeumannBoundary(nothing)), (j, axes(factor,2), NeumannBoundary(nothing))))
    get_u(i,j) = bc(u, ((i,axes(u,1),bnds[1]), (j,axes(u,2),bnds[2])))

    @inbounds for i in axes(u,1), j in axes(u,2)
        Dc = get_D(i,j)
        Dl = (0.5 * Dc + 0.5 * get_D(i-1,j))
        Dr = (0.5 * Dc + 0.5 * get_D(i+1,j))
        Db = (0.5 * Dc + 0.5 * get_D(i,j-1))
        Dt = (0.5 * Dc + 0.5 * get_D(i,j+1))

        uc = get_u(i,j)
        ul = get_u(i-1,j)
        ur = get_u(i+1,j)
        ub = get_u(i,j-1)
        ut = get_u(i,j+1)

        du[i,j] += (
            Dl * (ul - uc) / dV[1]^2 + 
            Dr * (ur - uc) / dV[1]^2 +
            Db * (ub - uc) / dV[2]^2 +
            Dt * (ut - uc) / dV[2]^2
            )
    end
    return nothing
end





const Dim = 2

p = (
    nx = 100,
    ny = 100,
    Lx = 1.0, 
    Ly = 1.0
)

function init(p)
    (;nx, ny, Lx, Ly) = p

    xs = LinRange(0, p.Lx, p.nx)
    ys = LinRange(0, p.Ly, p.ny)

    u0 = xs * ys'
end

u0 = init(p)