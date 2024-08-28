include("../src/AnisotropicParticles.jl")
include("parameterscan.jl")
include("plotutils.jl")

# config 
fn_p = "paper/inputs/base.toml"
datadir = "paper/data"
plotdir = "paper/plots"
videodir = "paper/videos"

n_reps = 8
p = loadparameters(fn_p)
ts = LinRange(0, p.t_end, 400)
redo = false

config = (;datadir, plotdir, videodir, redo)


X = ParameterRange(:D_x, L"D_x", 0.5 .^ [2,3,4,5,7,9])
data = parameterscan(p, X, n_reps, analyze_S2_traj(ts); changes = (:D_u => 0.0,))

data = stack(data)  # timesteps × n_reps × n_cases

trajectoryplot(p, X, ts, data, config)

X_vid = ParameterRange(X.sym, X.name, [1, 3])
savevideos(p, X_vid, config)

