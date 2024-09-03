struct State
    X
    theta
    t
end

function State(sol, t)
    (; X, theta) = unpack(sol(t))
    return State(X, theta, t)
end


function pack(X, theta)
    return [reinterpret(Float64, X); theta]
end

pack(s) = pack(s.X, s.theta)

function unpack(u)
    n = div(length(u), 3) * 2

    X = reinterpret(SVec2, @view u[1:n])
    theta = @view u[n+1:end]

    return (; X, theta)
end




#= 
    Apply periodic boundary conditions to positions
=#
@inline function proj(p, x; dom=SA[p.Lx, p.Ly])
    if p.periodic
        mod.(x, dom)
    else
        x
    end
end

#= 
    Apply periodic boundary conditions to geodesics/relative positions
=#
@inline function wrap(p, dx; dom=SA[p.Lx, p.Ly], inv_dom=inv.(dom))
    if p.periodic
        return @. dx - dom * round(dx * inv_dom)
    else
        return dx
    end
end


#############################
#initial random values for X and theta
#############################
function initstate(p, rng=Random.default_rng())
    (; N, Lx, Ly) = p

    # create initial positions and angles
    X = [SVec2(rand(rng) * Lx, rand(rng) * Ly) for _ in 1:N]
    theta = [rand(rng, FloatT) * pi for _ in 1:N]

    if hasproperty(p, :init)
        if p.init == :Havva_1
            for i in eachindex(X)
                if X[i][1] < Lx / 2
                    theta[i] = 0.0
                else
                    theta[i] = pi / 2.0
                end
            end
        elseif p.init == :Antoine_1
            for i in eachindex(X)
                if sin(X[i][2] / Ly * p.freq * 2π) < 0
                    theta[i] = 0.0
                else
                    theta[i] = pi / 2.0
                end
            end
        elseif p.init == :bump_1
            center = SA[p.Lx, p.Ly] ./ 2

            for i in eachindex(X)
                X[i] = randn(SVec2) .* sqrt(p.init_spread) .+ center
            end
        end
    end

    s = (; X, theta, t=FloatT(0.0))

    if use_gpu
        s = map(x -> cu(x), p)
    end

    return State(s.X, s.theta, s.t)
end



function potential(x, p)
    (; l, d) = p

    R = @SVector[x[1], x[2]]
    alpha = x[3]
    beta = x[4]

    sa, ca = sincos(alpha)
    sb, cb = sincos(beta)

    gamma_i = (l^2 - d^2) * @SMatrix[ca^2 ca*sa; ca*sa sa^2] + d^2 * I
    gamma_j = (l^2 - d^2) * @SMatrix[cb^2 cb*sb; cb*sb sb^2] + d^2 * I
    Σ = gamma_i + gamma_j

    expo = hasproperty(p, :expo) ? p.expo : FloatT(0.5)
    # expo = FloatT(0.5)

    return FloatT(1 / (4π)) * det(Σ)^(expo) * exp(-dot(R, inv(Σ) * R))
end

@kernel function particle_interaction!(p, X, theta, grid, dX, dTheta, dom, inv_dom)
    i = @index(Global)
    Xi = X[i]
    dXi = zero(Xi)
    dThetai = FloatT(0.0)

    N = length(X)
    cutoff² = p.cutoff^2

    for j in neighbours(grid, Xi, p.cutoff)
        Xij = wrap(p, Xi - X[j]; dom, inv_dom)
        d² = dist²(Xij)

        if 0 < d² < cutoff²
            alpha = theta[i]
            beta = theta[j]

            z = SA[Xij[1], Xij[2], alpha, beta]
            dzij = ForwardDiff.gradient(z -> potential(z, p), z)

            dXi += -1 / N * SA[dzij[1], dzij[2]]
            dThetai += -1 / N * dzij[3]
        end
    end

    dX[i] = p.mu * dXi
    dTheta[i] = p.lambda * dThetai
end


function sde_rhs!(du, u, p, t)
    (; particle_interaction_kernel, dom, inv_dom, grid) = p

    X, theta = unpack(u)
    dX, dTheta = unpack(du)  # make changes in place! 

    particle_interaction_kernel(p, X, theta, grid, dX, dTheta, dom, inv_dom, ndrange=length(X))
end

function sde_noise!(du, u, p, t)
    n = div(length(u), 3)
    du[1:2*n] .= p.coeff_x
    du[2*n+1:end] .= p.coeff_theta
end

function updategrid!(u, t, integrator)
    p = integrator.p
    X, theta = unpack(u)
    updatecells!(p.grid, X)
end

function steadynematicmean(integrator, abstol, reltol, min_t)

    if !isnothing(min_t) && t < min_t
        return false
    end

    _, theta = unpack(integrator.u)
    now = nematic_mean(theta)

    _, thetaprev = unpack(integrator.uprev)
    prev = nematic_mean(thetaprev)

    err = norm(prev - now) / abs(integrator.t - integrator.tprev)
    aim = max(abstol, reltol * norm(prev))

    return err < aim
end

function construct_particle_ode(p)

end


function simulate(p, s=initstate(p);
    rng=Random.default_rng(),
    seed=0,
    preview=false,
    fig=nothing,
    progress=true,
    steadystate=false,
    abstol=1e-8,
    reltol=1e-6,
    repetitions=1,
    callbacks=())

    # convert to GPU types...
    p = map(x -> x isa Int ? IntT(x) : x isa Float64 ? FloatT(x) : x, p)

    # set random seed
    Random.seed!(seed)

    # create a callback for the grid
    grid = BoundedGrid(p.cutoff, SA[p.Lx, p.Ly], s.X, IndexVecT)
    particle_interaction_kernel = particle_interaction!(get_backend(s.X), 64)

    dom = SVec2(p.Lx, p.Ly)
    inv_dom = FloatT.(1 ./ dom)

    p_gpu = map(x -> (x isa IntT || x isa FloatT || x isa Bool) ? x : nothing, p)

    p_sde = (;
        p_gpu...,
        particle_interaction_kernel,
        dom, inv_dom,
        grid,
        # noise  
        coeff_x=sqrt(2 * p.D_x),
        coeff_theta=sqrt(2 * p.D_u)
    )

    tspan = (FloatT(0.0), p.t_end)
    cb_grid = FunctionCallingCallback(updategrid!, func_everystep=true)
    cb_terminate = TerminateSteadyState(1e-8, 1e-6, steadynematicmean)
    # cb_showtime = PeriodicCallback((i) -> println("t = ", i.t), 10.0)

    cb = CallbackSet(cb_grid, (steadystate ? (cb_terminate,) : ())..., callbacks...)

    u0 = pack(s.X, s.theta)
    prob = SDEProblem(sde_rhs!, sde_noise!, u0, tspan, p_sde, callback=cb)

    if repetitions == 1
        return solve(prob, RKMil())
    else
        return [solve(prob, RKMil(), u0=pack(initstate(p))) for _ in 1:repetitions]
    end
end

function simulate_ensemble(p, n_reps; kwargs...)
    sols = simulate(p; kwargs..., repetitions=n_reps)

    if !(sols isa Vector)  # for the special case of only one repetition...
        sols = [sols]
    end

    return sols
end