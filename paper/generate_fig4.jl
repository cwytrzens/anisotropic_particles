
include("generate_base.jl")

n_reps = 6
n_reps_heatmap = 2
p = loadparameters(fn_p)
ts = LinRange(0, p.t_end, 400)
redo = true

config = (;datadir, plotdir, videodir, redo)




# Figure 2a
X = ParameterRange(:xi, L"\xi", 0.0:0.1:1.0 )
fn_data = joinpath(datadir, "gamma_f_$(X.sym).jld2")

data = if redo || !isfile(fn_data)
    data = parameterscan(p, X, n_reps, analyze_S2_traj(ts); extras = (p, x) -> (:lambda => p.lambda / (1 - p.chi^2)^(2x),))
    data = stack(data)  # timesteps × n_reps × n_cases

    JLD2.@save fn_data data p X ts
    data
else
    JLD2.@load fn_data data
    data 
end

trajectoryplot(p, X, ts, data, config; inds = [1,3,5,7,9,10,11])