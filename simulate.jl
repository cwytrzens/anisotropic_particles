include("definitions.jl")


p = loadparameters("params.toml")
p = @set p.N = 50
p = @set p.Lx = 20.0
p = @set p.Ly = 20.0

Random.seed!(0)  # set random seed
s = init(p)
ts, sol = simulate(s, p)

#result_1 = deepcopy( (;ts, sol) )
#result_2 = deepcopy( (;ts, sol) )

sum( norm(x - y) for (x,y) in zip(sol[end-1].X, result_2.sol[end-1].X))

####
# Single plot
####
# s_obs = Observable(s)
# fig = init_plot(s_obs, p)

####
# Make video 
####
# record(fig, "movie.mkv", eachindex(sol)) do i 
#     s_obs[] = sol[i]  # update state
# end


#####
# Interactive Window 
#####
fig = Figure()
sl = Slider(fig[2,1], range = eachindex(sol))
s_obs = @lift sol[$(sl.value)]

init_plot(s_obs, p, fig)
fig


#####
# Interactive Window, compare
#####
fig = Figure()
sl = Slider(fig[2,1], range = eachindex(sol))
s_obs = @lift sol[$(sl.value)]

init_plot(s_obs, p, fig[1,1])

s_obs_2 = @lift result_1.sol[$(sl.value)]

init_plot(s_obs_2, p, fig[1,2])
fig