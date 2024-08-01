include("definitions.jl")


p = loadparameters("params.toml")
p = @set p.N = 30_000
p = @set p.t_end = 5000.0

Random.seed!(0)  # set random seed
s = init(p)
ts, sol = simulate(s, p)

# benchmarks:
# - 1.056791 seconds (2.98 M allocations: 131.389 MiB)
# - 1.028484 seconds (2.98 M allocations: 131.385 MiB)

# SHT: 
# - 0.488747 seconds (15.01 k allocations: 16.006 MiB)

#result_1 = deepcopy( (;ts, sol) )
#result_2 = deepcopy( (;ts, sol) )

# sum( norm(x - y) for (x,y) in zip(sol[end-1].X, result_2.sol[end-1].X))

####
# Single plot
####
# s_obs = Observable(s)
# fig = init_plot(s_obs, p)

#####
# Interactive Window 
#####
fig = Figure()
sl = Slider(fig[2,1], range = eachindex(sol))
s_obs = @lift sol[$(sl.value)]
# s_obs = Observable(sol[1])  
init_plot(s_obs, p, fig[1,1])
fig

####
# Make video 
####
record(fig, "movie_antoine.mp4", eachindex(sol)[1:10:end]) do i 
    s_obs[] = sol[i]  # update state
    @show i
end



#####
# Interactive Window, compare
#####
# fig = Figure()
# sl = Slider(fig[2,1], range = eachindex(sol))
# s_obs = @lift sol[$(sl.value)]

# init_plot(s_obs, p, fig[1,1])

# s_obs_2 = @lift result_1.sol[$(sl.value)]

# init_plot(s_obs_2, p, fig[1,2])
# fig