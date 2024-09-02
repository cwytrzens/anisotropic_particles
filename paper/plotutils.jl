using CairoMakie, GLMakie

theme_mods = Theme(
    fontsize = 16,
    palette = (color = Makie.wong_colors(), linestyle = [:solid, :dash, :dot],),
    Lines = (
        cycle = Cycle([:color, :linestyle], covary = true,),
    ),
    Heatmap = (
        colormap = :thermal,
    ),
    Axis = (
        ygridvisible = true,
    )
)

video_theme = merge(theme_dark(), theme_latexfonts())
manuscript_theme = merge(theme_minimal(), theme_latexfonts(), theme_mods)


function movingwindow(x, window = 5)
    return [mean(@view x[max(i-window, 1):min(i+window, length(x))]) for i in eachindex(x)]
end

function trajectoryplot(p, X, ts, data, config; 
                            xlabel = L"t \text{ (time)}", 
                            ylabel = L"\gamma_f \text{ (order parameter)}",
                            labels = :normal,
                            reversed = false,
                            title_extra = "",
                            title = L"\textbf{Trajectories of order parameter } %$(title_extra)")
                            
    CairoMakie.activate!()
    fig = with_theme(manuscript_theme) do 
        fig = Figure(size= (420, 320))

    label_format = if labels == :base2
        (x) -> L"2^{%$(round(Int, log2(x)))}"
    elseif labels == :base2_rounded
        (x) -> L"2^{%$(round(log2(x),digits=1))}"
    elseif labels == :base10
        (x) -> L"10^{%$(round(Int, log10(x)))}"
    else
        x -> L"%$(x)"
    end

    labels = [ x != 0 ? label_format(x) : L"0" for x in X.range]

        ax = Axis(fig[1,1]; xlabel, ylabel, title, ygridvisible = true, yticks = 0:0.2:1.0)

        inds = reversed ? reverse(eachindex(X.range)) : eachindex(X.range)

        for i in inds
            mu = mean(data[:,:,i], dims=2)[:,1]
            st = std(data[:,:,i], dims=2)[:,1]
            up = movingwindow(mu + st / 2)
            down = movingwindow(mu - st / 2)
            
            band!(ax, ts, up, down, 
                color = Makie.wong_colors(0.2)[i],
                label = labels[i])
        end
            

        for i in inds
            mu = mean(data[:,:,i], dims=2)[:,1]
            mu = movingwindow(mu)
            lines!(ax, ts, mu, 
                color = Cycled(i), linestyle = Cycled(i),
                linewidth = 1.5,label = labels[i])
        end

        Legend(fig[1,2], ax, X.name, merge = true, framevisible = false)

        xlims!(ax, 0, ts[end])
        ylims!(ax, 0, 1.)
        
        fig
    end
    save(joinpath(config.plotdir, "gamma_f_time_vs_$(X.sym).png"), fig)    
    save(joinpath(config.plotdir, "gamma_f_time_vs_$(X.sym).eps"), fig)
    fig
end


function plotheatmap(X, Y, data; xscale = identity, yscale = identity)
    CairoMakie.activate!()
    with_theme(manuscript_theme) do 
        fig = Figure(size = (420, 320))
        ax = Axis(fig[1,1]; xscale, yscale, 
                    xlabel = X.name, ylabel = Y.name, 
                    title = L"\gamma_f \text{  at } t = %$(p.t_end)")
        
        hm = heatmap!(X.range, Y.range, mean(data, dims = 1)[1,:,:], 
                    colorrange = (0, 1))
        Colorbar(fig[1,2], hm, label = L"\gamma_f")

        save(joinpath(plotdir, "gamma_f_$(X.sym)_$(Y.sym).png"), fig)
        save(joinpath(plotdir, "gamma_f_$(X.sym)_$(Y.sym).eps"), fig)
        fig
    end
end





# function to quickly see how the simulation for a parameterfile looks like 
function interact(fn::String, changes::Tuple = ())
    GLMakie.activate!()
    p = loadparameters(fn)
    p = updateparameters(p, changes)
    sol = simulate(p; seed = rand(1:1000))
    fig = AnisotropicParticles.interact(sol, p) 
    display(fig)
    CairoMakie.activate!()

    return fig, sol, p
end

function savevideo(fn, p;
                        changes::Tuple = (), 
                        p_mod = updateparameters(p, changes),
                        sol = simulate(p_mod), 
                        n_frames = 30 * 10,
                        frames = LinRange(0, sol.t[end], n_frames),
                        resolution = (1280, 1280))

    s_obs = Observable(State(sol, 0.0))


    GLMakie.activate!()
    with_theme(merge(theme_dark(), theme_latexfonts())) do 
        fig = Figure(size = resolution)
        AnisotropicParticles.init_plot(s_obs, p_mod, fig[1,1])
        display(fig)

        prog = Progress(length(frames), 1)
        record(fig, fn, frames) do t 
            s_obs[] = State(sol, t)
            next!(prog)
        end            
    end
    CairoMakie.activate!()

    return p_mod, sol
end
