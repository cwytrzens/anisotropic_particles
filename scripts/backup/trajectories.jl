include("../definitions.jl")
using GLMakie.Makie.LaTeXStrings: latexstring


fig, sol = interact("paper/inputs/base.toml")


# main parameter file 
p_base = loadparameters("paper/inputs/base.toml")
p_base = rescale(p_base)
sol = simulate(p_base; seed = 3)
# statistical number of repetitions 
n_reps = 8



X = ParameterRange(:D_x, L"D_x", 2.0 .^ (-11:2:-1) )
# Y = ParameterRange(:lambda, L"\lambda", LinRange(1, 400, 16))

data = Array{Any}(undef, length(X.range))

t_end = p_base.t_end
analyze(sol) = [nematic_mean(State(sol, t).theta) for t in LinRange(0,t_end,500)]

fn = "$(X.sym)_nematic_alignment_trajectory"


@progress for (i, x) in enumerate(X.range) 
    p_ = modify(p_base, (X.sym => x,))
    sols = simulate_ensemble(p_, n_reps)

    samples = analyze.(sols)
    data[i] = samples

    println(X.sym, ": ", x, "\toutput: ", mean(samples))
end


nn_trajs = stack( [norm.(d[i]) for d in data, i in 1:n_reps] )  # time x param x reps 
nn_mean = mean(nn_trajs, dims = 3)[:,:,1]
nn_std = std(nn_trajs, dims = 3)[:,:,1]

JLD2.@save "$(fn).jld2" nn_mean nn_std p_base t_end X 

using CairoMakie

CairoMakie.activate!()
cols = Makie.wong_colors()
cols_transp = Makie.wong_colors(0.2)

fn = "D_x_nematic_alignment_trajectory"
JLD2.@load "$(fn).jld2" nn_mean nn_std p_base t_end X 
with_theme(merge(theme_latexfonts(), theme_minimal())) do
    
    fig = Figure( size = (420, 320))

    ax = Axis(fig[1,1], 
        xlabel = L"t \text{ (time)}", 
        ylabel = L"\Vert \Omega \Vert \text{ (nematic alignment)}", 
        title = L"\textbf{Nematic alignment trajectory } (D_u = 0)",
        yticks = 0:0.2:1.0, ygridvisible = true)

    for (i, x) in reverse(collect(enumerate(X.range)))
        band!(ax, LinRange(0,t_end,500), nn_mean[:,i] - nn_std[:,i] ./ 2, nn_mean[:,i] + nn_std[:,i] ./ 2,
            color = cols_transp[i],
            label = latexstring("2^{$(Int(log2(x)))}"))
    end

    for (i, x) in reverse(collect(enumerate(X.range)))
        lines!(ax, LinRange(0,t_end,500), nn_mean[:,i],
            linewidth = 1.5, colormap = :tab10, color = cols[i],
            label = latexstring("2^{$(Int(log2(x)))}"))
    end

    Legend(fig[1,2], ax, L"D_x", merge = true, framevisible = false)

    xlims!(ax, 0, t_end)
    ylims!(ax, 0, 1.)

    save("$(fn).png", fig)
    save("$(fn).eps", fig)

    fig
end
















p_base =rescale(loadparameters("paper/inputs/figure_2_alignment.toml"))

X = ParameterRange(:D_u, L"D_u", [0.0; 2.0 .^ [-11, -10.5, -10, -9.5, -9]] )
# Y = ParameterRange(:lambda, L"\lambda", LinRange(1, 400, 16))

data = Array{Any}(undef, length(X.range))

t_end = p_base.t_end
analyze(sol) = [nematic_mean(State(sol, t).theta) for t in LinRange(0,t_end,500)]

fn = "$(X.sym)_nematic_alignment_trajectory"

# GLMakie.activate!()
# p_base = loadparameters("paper/figure_2_alignment.toml")
# p_base = rescale(p_base)
# sol = simulate(p_base; seed = 3)
# interact(sol, p_base) 

@progress for (i, x) in enumerate(X.range) 

    p_ = modify(p_base, (X.sym => x,))

    sols = simulate_ensemble(p_, n_reps)

    samples = analyze.(sols)

    data[i] = samples

    println(X.sym, ": ", x, "\toutput: ", mean(samples))
end


nn_trajs = stack( [norm.(d[i]) for d in data, i in 1:n_reps] )  # time x param x reps 
nn_mean = mean(nn_trajs, dims = 3)[:,:,1]
nn_std = std(nn_trajs, dims = 3)[:,:,1]

JLD2.@save "paper/data/$(fn).jld2" nn_mean nn_std p_base t_end X 

using CairoMakie

CairoMakie.activate!()
cols = Makie.wong_colors()
cols_transp = Makie.wong_colors(0.2)


JLD2.@load "paper/data/$(fn).jld2" nn_mean nn_std p_base t_end X 
with_theme(merge(theme_latexfonts(), theme_minimal())) do
    
    fig = Figure( size = (420, 320))

    labels = [ x != 0 ? latexstring("2^{$(round(log2(x),digits=1))}") : L"0" for x in X.range]

    ax = Axis(fig[1,1], 
        xlabel = L"t \text{ (time)}", 
        ylabel = L"\Vert \Omega \Vert \text{ (nematic alignment)}", 
        title = L"\textbf{Nematic alignment trajectory } (D_x = 2^{-3})",
        yticks = 0:0.2:1.0, ygridvisible = true)

    for (i, x) in collect(enumerate(X.range))
        band!(ax, LinRange(0,t_end,500), nn_mean[:,i] - nn_std[:,i] ./ 2, nn_mean[:,i] + nn_std[:,i] ./ 2,
            color = cols_transp[i],
            label = labels[i])
    end

    for (i, x) in collect(enumerate(X.range))
        lines!(ax, LinRange(0,t_end,500), nn_mean[:,i],
            linewidth = 1.5, color = cols[i],
            label = labels[i])
    end

    Legend(fig[1,2], ax, L"D_x", merge = true, framevisible = false)

    xlims!(ax, 0, t_end)
    ylims!(ax, 0, 1.)

    save("paper/plots/$(fn).png", fig)
    save("paper/plots/$(fn).eps", fig)

    fig
end











p_base =rescale(loadparameters("paper/inputs/base.toml"))
interact("paper/inputs/base.toml", (:d => 0.1,))
X = ParameterRange(:d, L"d", [0.1, 0.2, 0.4, 0.8] )
# Y = ParameterRange(:lambda, L"\lambda", LinRange(1, 400, 16))

data = Array{Any}(undef, length(X.range))

t_end = p_base.t_end
analyze(sol) = [nematic_mean(State(sol, t).theta) for t in LinRange(0,t_end,500)]

fn = "$(X.sym)_nematic_alignment_trajectory"

# GLMakie.activate!()
# p_base = loadparameters("paper/figure_2_alignment.toml")
# p_base = rescale(p_base)
# sol = simulate(p_base; seed = 3)
# interact(sol, p_base) 
n_reps = 4
@progress for (i, x) in enumerate(X.range) 

    p_ = modify(p_base, (X.sym => x,))

    sols = simulate_ensemble(p_, n_reps)

    samples = analyze.(sols)

    data[i] = samples

    println(X.sym, ": ", x, "\toutput: ", mean(samples))
end


nn_trajs = stack( [norm.(d[i]) for d in data, i in 1:n_reps] )  # time x param x reps 
nn_mean = mean(nn_trajs, dims = 3)[:,:,1]
nn_std = std(nn_trajs, dims = 3)[:,:,1]

JLD2.@save "paper/data/$(fn).jld2" nn_mean nn_std p_base t_end X 

using CairoMakie

CairoMakie.activate!()
cols = Makie.wong_colors()
cols_transp = Makie.wong_colors(0.2)


JLD2.@load "paper/data/$(fn).jld2" nn_mean nn_std p_base t_end X 
with_theme(merge(theme_latexfonts(), theme_minimal())) do
    
    fig = Figure( size = (420, 320))

    labels = [latexstring("$(round((1^2-x^2)/(1^2+x^2),digits=1))") for x in X.range]

    ax = Axis(fig[1,1], 
        xlabel = L"t \text{ (time)}", 
        ylabel = L"S_2 \text{ (order parameter)}", 
        title = L"\textbf{Trajectories of order parameter}",
        yticks = 0:0.2:1.0, ygridvisible = true)

    for (i, x) in collect(enumerate(X.range))
        band!(ax, LinRange(0,t_end,500), nn_mean[:,i] - nn_std[:,i] ./ 2, nn_mean[:,i] + nn_std[:,i] ./ 2,
            color = cols_transp[i],
            label = labels[i])
    end

    for (i, x) in collect(enumerate(X.range))
        lines!(ax, LinRange(0,t_end,500), nn_mean[:,i],
            linewidth = 1.5, color = cols[i],
            label = labels[i])
    end

    Legend(fig[1,2], ax, L"\chi", merge = true, framevisible = false)

    xlims!(ax, 0, t_end)
    ylims!(ax, 0, 1.)

    save("paper/plots/$(fn).png", fig)
    save("paper/plots/$(fn).eps", fig)

    fig
end





p_base =rescale(loadparameters("paper/inputs/base.toml"))
interact("paper/inputs/base.toml", (:d => 0.1,))
X = ParameterRange(:density, L"\overline{\rho}", [0.5, 1.0, 1.5, 2.0] )
# Y = ParameterRange(:lambda, L"\lambda", LinRange(1, 400, 16))

data = Array{Any}(undef, length(X.range))

t_end = p_base.t_end
analyze(sol) = [nematic_mean(State(sol, t).theta) for t in LinRange(0,t_end,500)]

fn = "$(X.sym)_nematic_alignment_trajectory"

# GLMakie.activate!()
# p_base = loadparameters("paper/figure_2_alignment.toml")
# p_base = rescale(p_base)
# sol = simulate(p_base; seed = 3)
# interact(sol, p_base) 
n_reps = 4
@progress for (i, x) in enumerate(X.range) 

    p_ = modify(p_base, (X.sym => x,))

    sols = simulate_ensemble(p_, n_reps)

    samples = analyze.(sols)

    data[i] = samples

    println(X.sym, ": ", x, "\toutput: ", mean(samples))
end


nn_trajs = stack( [norm.(d[i]) for d in data, i in 1:n_reps] )  # time x param x reps 
nn_mean = mean(nn_trajs, dims = 3)[:,:,1]
nn_std = std(nn_trajs, dims = 3)[:,:,1]

JLD2.@save "paper/data/$(fn).jld2" nn_mean nn_std p_base t_end X 

using CairoMakie

CairoMakie.activate!()
cols = Makie.wong_colors()
cols_transp = Makie.wong_colors(0.2)


JLD2.@load "paper/data/$(fn).jld2" nn_mean nn_std p_base t_end X 
with_theme(merge(theme_latexfonts(), theme_minimal())) do
    
    fig = Figure( size = (420, 320))

    labels = [latexstring("$(round((1^2-x^2)/(1^2+x^2),digits=1))") for x in X.range]

    ax = Axis(fig[1,1], 
        xlabel = L"t \text{ (time)}", 
        ylabel = L"\overline{\rho} \text{ (order parameter)}", 
        title = L"\textbf{Trajectories of order parameter}",
        yticks = 0:0.2:1.0, ygridvisible = true)

    for (i, x) in collect(enumerate(X.range))
        band!(ax, LinRange(0,t_end,500), nn_mean[:,i] - nn_std[:,i] ./ 2, nn_mean[:,i] + nn_std[:,i] ./ 2,
            color = cols_transp[i],
            label = labels[i])
    end

    for (i, x) in collect(enumerate(X.range))
        lines!(ax, LinRange(0,t_end,500), nn_mean[:,i],
            linewidth = 1.5, color = cols[i],
            label = labels[i])
    end

    Legend(fig[1,2], ax, L"\chi", merge = true, framevisible = false)

    xlims!(ax, 0, t_end)
    ylims!(ax, 0, 1.)

    save("paper/plots/$(fn).png", fig)
    save("paper/plots/$(fn).eps", fig)

    fig
end