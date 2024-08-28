using CairoMakie, GLMakie

theme_mods = Theme(
    palette = (color = Makie.wong_colors(), linestyle = [:solid, :dash, :dot],),
    Lines = (
        cycle = Cycle([:color, :linestyle], covary = true,),
    )
)

video_theme = merge(theme_dark(), theme_latexfonts())
manuscript_theme = merge(theme_minimal(), theme_latexfonts(), theme_mods)


function trajectoryplot(p, X, ts, data, config)
    with_theme(manuscript_theme) do 
        fig = Figure(size= (420, 320))
        ax = Axis(fig[1,1], 
                xlabel = L"t \text{ (time)}",
                ylabel = L"S2", 
                title = L"\text{Trajectories of order parameter } (D_u = 0)")
        
        for (i, x) in enumerate(X.range)
            band!(ax, ts, 
                mean(data[:,:,i], dims=2)[:,1] - std(data[:,:,i], dims=2)[:,1] ./ 2, 
                mean(data[:,:,i], dims=2)[:,1] + std(data[:,:,i], dims=2)[:,1] ./ 2, 
                color = Makie.wong_colors(0.2)[i],
                label = latexstring("2^{$(Int(log2(x)))}"))
        end
            

        for (i, x) in enumerate(X.range)
            lines!(ax, ts, mean(data[:,:,i], dims=2)[:,1], color = Cycled(i), linestyle = Cycled(i),
                linewidth = 1.5,label = latexstring("2^{$(Int(log2(x)))}"))
        end

        Legend(fig[1,2], ax, X.name, merge = true, framevisible = false)

        xlims!(ax, 0, ts[end])
        ylims!(ax, 0, 1.)
        
        fig
    end
    # save(joinpath(config.plotdir, "S2_vs_$(X.sym).png"), fig)
end


# function to quickly see how the simulation for a parameterfile looks like 
function interact(fn::String, changes::Tuple = ())
    GLMakie.activate!()
    p = loadparameters(fn)
    p = modify(p, changes)
    sol = simulate(p; seed = rand(1:1000))
    fig = interact(sol, p) 
    display(fig)
    CairoMakie.activate!()

    return fig, sol, p
end

function savevideo(fn, p;
                        changes::Tuple = (), 
                        sol = simulate(modify(p, changes)), 
                        frames = LinRange(0, sol.t[end], 30 * 10))

    s_obs = Observable(State(sol, 0.0))


    GLMakie.activate!()
    with_theme(merge(theme_dark(), theme_latexfonts())) do 
        fig = Figure()
        init_plot(s_obs, p, fig[1,1])
        display(fig)
        record(fig, fn, frames) do t 
            s_obs[] = State(sol, t)
        end            
    end
    CairoMakie.activate!()
end
