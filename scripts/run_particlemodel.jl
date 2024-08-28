include("../src/AnisotropicParticles.jl")


p = loadparameters("params.toml")
p = @set p.N = 30_000
p = @set p.t_end = 5000.0

Random.seed!(0)  # set random seed
s = initstate(p)
sol = simulate(p, s)

####
# Single plot
####
# fig = init_plot(State(sol, sol.t[end]), p)

#####
# Interactive Window 
#####
fig = Figure()
ts = LinRange(0, sol.t[end], 100)
sl = Slider(fig[2,1], range = ts, startvalue = ts[end])
s_obs = @lift State(sol, $(sl.value))
init_plot(s_obs, p, fig[1,1])
fig

####
# Make video 
####
fig = Figure()
s_obs = Observable(State(sol, 0.0))
init_plot(s_obs, p, fig[1,1])

ts = LinRange(0, sol.t[end], 100)
prog = Progress(length(ts), 1.0)
record(fig, "movie.mp4", ts) do t 
    s_obs[] = State(sol, t)
    next!(prog)
end

