using CairoMakie, GLMakie

video_theme = merge(theme_dark(), theme_latexfonts())
manuscript_theme = merge(theme_minimal(), theme_latexfonts()) 
defaultcolor(i, alpha = 1.0) = Makie.wong_colors(alpha)[i]





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
