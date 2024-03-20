include("definitions.jl")


#############################
# parameters
#############################
p = let 
    l = 1.5
    d = 0.2 
    cutoff = l * 2
    
    (
        N = 5000,  # total number of particles 
        l = l,  # length of major axis of ellipse 
        d = d,  # length of minor axis of ellipse
        t_step = 0.1,
        t_save = 0.0,
        t_start = 0.0,
        t_end = 500.0,
        Lx = 50.0,
        Ly = 50.0,
        mu = 100.0,
        lambda = 1.0,
        D_x  = 0 * 0.01,
        D_u =  0 * 0.001,
        periodic = true,
        cutoff = cutoff
    )
end 

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