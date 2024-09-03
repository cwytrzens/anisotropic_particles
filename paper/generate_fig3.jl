
include("generate_base.jl")

n_reps = 3
n_reps_heatmap = 2
p = loadparameters(fn_p)
ts = LinRange(0, p.t_end, 400)
redo = true

config = (;datadir, plotdir, videodir, redo)





# Figure 2a
X = ParameterRange(:lambda, L"\lambda", [2^4, 2^5, 2^6, 2^6.5, 2^7, 2^7.5, 2^8] )
fn_data = joinpath(datadir, "gamma_f_$(X.sym)_3.jld2")

data = if redo || !isfile(fn_data)
    data = parameterscan(p, X, 5, analyze_S2_traj(ts); changes = (:mu => 2^12,))
    data = stack(data)  # timesteps × n_reps × n_cases

    JLD2.@save fn_data data p X ts
    data
else
    JLD2.@load fn_data data
    data 
end

trajectoryplot(p, X, ts, data, config; labels = :base2_rounded, title_extra = L"(\mu = 2^{12})", reversed = true)


# Figure 2b
X = ParameterRange(:mu, L"\mu", [0, 2^4, 2^12, 2^13, 2^14, 2^15, 2^16] )
fn_data = joinpath(datadir, "gamma_f_$(X.sym).jld2")


data = if redo || !isfile(fn_data)
    data = parameterscan(p, X, 3, analyze_S2_traj(ts))
    data = stack(data)  # timesteps × n_reps × n_cases

    JLD2.@save fn_data data p X ts
    data
else
    JLD2.@load fn_data data
    data 
end

trajectoryplot(p, X, ts, data, config; labels = :base2_rounded, title_extra = L"(\lambda = 150)")



# Figure 1c
X = ParameterRange(:lambda, L"\lambda", 2.0 .^ LinRange(6,8,16))
Y = ParameterRange(:mu, L"\mu", 2.0 .^ (11:0.5:15))
fn_data = joinpath(datadir, "gamma_f_$(X.sym)_$(Y.sym)_2.jld2")

data = if redo || !isfile(fn_data)
    data = parameterscan(p, X, Y, 8, analyze_S2_fixedtime(p.t_end); changes = ())
    data = stack(data)  # n_reps × n_cases × n_cases

    JLD2.@save fn_data data p X ts
    data
else
    JLD2.@load fn_data data
    data 
end


plotheatmap(X, Y, data; xscale = log2, yscale = log2)






# Figure 1c
X = ParameterRange(:lambda, L"\lambda", 2.0 .^ LinRange(6,8,4))
Y = ParameterRange(:mu, L"\mu", 2.0 .^ (11:1.0:15))
fn_data = joinpath(datadir, "gamma_f_$(X.sym)_$(Y.sym)_3.jld2")

data = if redo || !isfile(fn_data)
    data = parameterscan(p, X, Y, 1, analyze_S2_fixedtime(p.t_end); changes = ())
    data = stack(data)  # n_reps × n_cases × n_cases

    JLD2.@save fn_data data p X ts
    data
else
    JLD2.@load fn_data data
    data 
end


plotheatmap(X, Y, data; xscale = log2, yscale = log2)


