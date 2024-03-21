include("definitions.jl")



p = loadparameters("params.toml")
s = init(p)
ts, sol = simulate(s, p)




s_obs = Observable(s)
fig = init_plot(s_obs, p)


######################################
# Make video 
######################################
# record(fig, "movie.mkv", eachindex(sol)) do i 
#     s_obs[] = sol[i]  # update state
# end


#######################################
# Interactive Window 
#######################################

fig = Figure()
sl = Slider(fig[2,1], range = eachindex(sol))
s_obs = @lift sol[$(sl.value)]

init_plot(s_obs, p, fig)
fig