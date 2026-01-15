
p = loadparameters("params.toml")
cache = create_interpolations(p)

rho_min = minimum(cache.data.rhos)
u0 = init(p)
tspan = (p.t_start, p.t_end)
p_ode = (;p..., cache)
prob = ODEProblem(rhs!, u0, tspan, p_ode)

sol = solve(prob, ROCK2(), progress = true)


begin
    fig = Figure() 
    sl = Slider(fig[2,1], range = eachindex(sol))
    sol_obs = @lift sol[$(sl.value)]

    ax = Axis(fig[1,1], title = @lift string("t = ", sol.t[$(sl.value)]))

    heatmap!(sol_obs, colorrange = (minimum(sol), maximum(sol)))

   fig
end





# create video
# begin
#     fig = Figure() 
#     title_str = Observable("t = ")
#     ax = Axis(fig[1,1], title = title_str)
#     sol_obs = Observable(sol[end])

#     heatmap!(sol_obs, colorrange = (minimum(sol), maximum(sol)))

#     record(fig, "plots/video_alpha=$(round(1 - p.chi^2, digits=2)).mp4", LinRange(p.t_start, p.t_end, 360)) do t 
#         title_str[] = string("t = ", round(t, digits = 4))
#         sol_obs[] = sol(t)
#     end
# end