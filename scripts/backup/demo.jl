include("definitions.jl")
using StochasticDiffEq, ProgressLogging
using OrdinaryDiffEq



function tostate(sol, ind)
    return (; unpack(sol[ind])..., t = sol.t[ind])
end

function attime(sol, t)
    s = sol(t) 
    return (;unpack(s)..., t = t)
end


using OrdinaryDiffEq, DiffEqCallbacks




p = loadparameters("params.toml")
p = rescale(p)
# convert to GPU types...


p = map( x -> x isa Int ? IntT(x) : x isa Float64 ? FloatT(x) : x, p)


ts = [p.t_start]
s = initstate(p)


u0 = pack(s.X, s.theta)

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
    coeff_x = sqrt(2 * p.D_x), 
    coeff_theta = sqrt(2 * p.D_u)
)

tspan = (FloatT(0.0), p.t_end)
cb_grid = FunctionCallingCallback(updategrid!, func_everystep = true)
cb_terminate = TerminateSteadyState(1e-8, 1e-6, steadynematicmean)
cb_showtime = PeriodicCallback((i) -> println("t = ", i.t), 10.0)

# cb = CallbackSet(cb_grid)

s_obs = with_theme(theme_dark()) do 
    fig = Figure()
    s_obs = Observable((;unpack(u0)..., t = 0.0))
    init_plot(s_obs, p, fig[1,1])
    display(fig)
    return s_obs
end
cb_updateplot = PeriodicCallback( (i) -> s_obs[] = (; unpack(i.u)..., t = i.t), 100.0)
cb = CallbackSet(cb_grid, cb_terminate, cb_updateplot)

u0 = pack(s.X, s.theta)
prob = SDEProblem(sde_rhs!, sde_noise!, u0, tspan, p_sde, callback = cb)

println("Particle area per domain area = ", (pi * p.d * p.l * p.N) / (p.Lx * p.Ly))


@time sol = solve(prob, RKMil());

# 1.6


# 
# # set_theme!()
with_theme(theme_dark()) do
    fig = Figure()
    sl = Slider(fig[2,1], range = eachindex(sol))
    s_obs = @lift tostate(sol, $(sl.value))
    # s_obs = Observable(sol[1])  
    init_plot(s_obs, p, fig[1,1])
    fig
end


with_theme(theme_dark()) do 
    fig = Figure() 
    ax = Axis(fig[1,1])
    ts = LinRange(sol.t[1], sol.t[end], 1000)
    nn = norm.(nematic_mean.(attime(sol, t).theta for t in ts))

    lines!(ts, nn)
    fig
end