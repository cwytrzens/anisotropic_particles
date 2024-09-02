
include("generate_base.jl")

n_reps = 8
n_reps_heatmap = 2
p = loadparameters(fn_p)
ts = LinRange(0, p.t_end, 400)
redo = true

config = (;datadir, plotdir, videodir, redo)





# Figure 2a
X = ParameterRange(:chi, L"\chi", [0.98, 0.9, 0.85,0.8, 0.5, 0.0] )
fn_data = joinpath(datadir, "gamma_f_$(X.sym).jld2")

data = if redo || !isfile(fn_data)
    data = parameterscan(p, X, n_reps, analyze_S2_traj(ts))
    data = stack(data)  # timesteps × n_reps × n_cases

    JLD2.@save fn_data data p X ts
    data
else
    JLD2.@load fn_data data
    data 
end

trajectoryplot(p, X, ts, data, config)


# Figure 2b
X = ParameterRange(:density, L"\overline{\rho}", [0.5, 0.75, 1.0, 1.5, 2.0] )
fn_data = joinpath(datadir, "gamma_f_$(X.sym).jld2")

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
X = ParameterRange(:chi, L"\chi", LinRange(0.8, 0.98, 16))
Y = ParameterRange(:density, L"\overline{\rho}", LinRange(0.5, 1.5, 16))
fn_data = joinpath(datadir, "gamma_f_$(X.sym)_$(Y.sym)_2.jld2")

data = if redo || !isfile(fn_data)
    data = parameterscan(p, X, Y, 3, analyze_S2_fixedtime(p.t_end); changes = ())
    data = stack(data)  # n_reps × n_cases × n_cases

    JLD2.@save fn_data data p X ts
    data
else
    JLD2.@load fn_data data
    data 
end


plotheatmap(X, Y, data; xscale = identity, yscale = identity)
