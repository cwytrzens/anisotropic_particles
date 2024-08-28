include("../src/AnisotropicParticles.jl")
include("parameterscan.jl")
include("plotutils.jl")

# config 
fn_p = "paper/inputs/base.toml"
datadir = "paper/data"
plotdir = "paper/plots"
videodir = "paper/videos"

n_reps = 8
n_reps_heatmap = 6
p = loadparameters(fn_p)
ts = LinRange(0, p.t_end, 400)
redo = false

config = (;datadir, plotdir, videodir, redo)



# Figure 1a
X = ParameterRange(:D_x, L"D_x", 0.5 .^ [2,3,4,5,7,9])
fn_data = joinpath(datadir, "S2_$(X.sym).jld2")

data = if redo || !isfile(fn_data)
    data = parameterscan(p, X, 8, analyze_S2_traj(ts); changes = (:D_u => 0.0,))
    data = stack(data)  # timesteps × n_reps × n_cases

    JLD2.@save fn_data data p X ts
    data
else
    JLD2.@load fn_data data
    data 
end

trajectoryplot(p, X, ts, data, config)



# Figure 1b
X = ParameterRange(:D_u, L"D_u", 0.5 .^ [11,10,9.5,9])
fn_data = joinpath(datadir, "S2_$(X.sym).jld2")

data = if redo || !isfile(fn_data)
    data = parameterscan(p, X, 8, analyze_S2_traj(ts))
    data = stack(data)  # timesteps × n_reps × n_cases

    JLD2.@save fn_data data p X ts
    data
else
    JLD2.@load fn_data data
    data 
end

trajectoryplot(p, X, ts, data, config)


# Figure 1c
X = ParameterRange(:D_x, L"D_x", 0.5 .^ LinRange(0,5,16))
Y = ParameterRange(:D_u, L"D_u", 0.5 .^ LinRange(9,14,16))
fn_data = joinpath(datadir, "S2_$(X.sym)_$(Y.sym).jld2")

data = if redo || !isfile(fn_data)
    data = parameterscan(p, X, Y, n_reps_heatmap, analyze_S2_fixedtime(p.t_end); changes = ())
    data = stack(data)  # n_reps × n_cases × n_cases

    JLD2.@save fn_data data p X ts
    data
else
    JLD2.@load fn_data data
    data 
end



using JLD2
JLD2.@save joinpath(datadir, "S2_DxDu_2.jld2") data p X Y

begin
    CairoMakie.activate!()
    with_theme(manuscript_theme) do 
        fig = Figure(size = (420, 320))
        ax = Axis(fig[1,1], xscale = log2, yscale = log2, xlabel = X.name, ylabel = Y.name, title = "S2 at t = $(p.t_end)")

        
        hm = heatmap!(X.range, Y.range, mean(data, dims = 1)[1,:,:], colormap = :plasma, colorrange = (0, 1))
        Colorbar(fig[1,2], hm, label = "S2")

        save(joinpath(plotdir, "S2_DxDu.png"), fig)
        save(joinpath(plotdir, "S2_DxDu.eps"), fig)
        fig
    end
end