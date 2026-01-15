
include("generate_base.jl")

n_reps = 14
n_reps_heatmap = 6
p = loadparameters(fn_p)
ts = LinRange(0, p.t_end, 400)
redo = true

config = (;datadir, plotdir, videodir, redo)



# Figure 1a
X = ParameterRange(:D_x, L"D_x", 0.5 .^ [2,3,4,5,6,7])
fn_data = joinpath(datadir, "gamma_f_$(X.sym).jld2")

data = if redo || !isfile(fn_data)
    data = parameterscan(p, X, 8, analyze_S2_traj(ts); changes = (:D_u => 0.0,))
    data = stack(data)  # timesteps × n_reps × n_cases

    JLD2.@save fn_data data p X ts
    data
else
    JLD2.@load fn_data data
    data 
end

trajectoryplot(p, X, ts, data, config; labels = :base2, title_extra = "(D_u = 0)")



# Figure 1b
X = ParameterRange(:D_u, L"D_u", 0.5 .^ [Inf,11,10,9.75,9.5,9])
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

trajectoryplot(p, X, ts, data, config; labels = :base2_rounded, title_extra = "(D_x = 2^{-3})")


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

plotheatmap(X, Y, data; xscale = log2, yscale = log2)