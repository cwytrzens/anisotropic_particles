# plotting function
function ellipse(X, theta, p)
    (;l, d) = p
    s, c = sincos(theta) 
    R = @SMatrix[ c -s; s c ]
    return Polygon(
        [Point2f(X + R * SVec2( l * cos(t + theta), d * sin(t + theta) )) for t in LinRange(0,2π,20)]
        )
end



############################################
# Create plot 
############################################
init_plot(s, p, fig_pos = Figure()[1,1]) = init_plot(Observable(s), p, fig_pos)

function init_plot(s::Observable, p, fig_pos = Figure()[1,1])
    ax = Axis(fig_pos, aspect = DataAspect(), title = @lift string("time = ", round($s.t, digits = 1)))

    plotparticles!(ax, s)
    # X = @lift map( x -> Point2f(proj(p,x)), $s.X) 
    # angles = @lift mod.($s.theta, π)

    # E = @lift [ ellipse($X[i], $angles[i], p) for i in 1:p.N ]
    # poly!(E, color = angles, alpha = 0.5, colorrange = (0.0, π), colormap = :cyclic_mygbm_30_95_c78_n256_s25) # :phase, :twilight,  :cyclic_mygbm_30_95_c78_n256_s25   

    current_figure() 
end

plotparticles!(p, s; kwargs...) = plotparticles!(current_axis(), p, s; kwargs...)
function plotparticles!(ax, p, s::Observable{State}; offset = (0,0), kwargs...)
    shift = SA[p.Lx, p.Ly] .* offset
    X = @lift map( x -> Point2f(proj(p,x) + shift), $s.X) 
    angles = @lift mod.($s.theta, π)
    E = @lift [ ellipse($X[i], $angles[i], p) for i in 1:p.N ]
    poly!(ax, E; color = angles, alpha = 0.5, colorrange = (0.0, π), colormap = :cyclic_mygbm_30_95_c78_n256_s25, kwargs...) # :phase, :twilight,  :cyclic_mygbm_30_95_c78_n256_s25   
end

function interact(sol, p, darktheme = true)

    with_theme(darktheme ? theme_dark() : theme_light()) do 
        fig = Figure(size = (1000, 1000))
        sl = Slider(fig[2,1], range = LinRange(0, sol.t[end], 1000), startvalue = sol.t[end])
        s_obs = @lift State(sol, $(sl.value))
        init_plot(s_obs, p, fig[1,1])
        fig
    end
end