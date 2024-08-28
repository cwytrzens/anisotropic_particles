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
sl = Slider(fig[2,1], range = LinRange(0, sol.t[end], 100), startvalue = sol.t[end])
s_obs = @lift State(sol, $(sl.value))
init_plot(s_obs, p, fig[1,1])
fig

####
# Make video 
####
# record(fig, "movie_antoine.mp4", eachindex(sol)[1:10:end]) do i 
#     s_obs[] = sol[i]  # update state
#     @show i
# end

