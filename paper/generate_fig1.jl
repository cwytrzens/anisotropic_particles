include("../src/AnisotropicParticles.jl")
include("parameterscan.jl")
include("plotutils.jl")

# config 
fn_p = "paper/inputs/base.toml"
datadir = "paper/data"
plotdir = "paper/plots"
videodir = "paper/videos"

n_reps = 4
p = loadparameters(fn)
ts = LinRange(0, p.t_end, 400)
redo = false

config = (;datadir, plotdir, videodir, redo)

X = ParameterRange(:D_u, L"D_u", [1,2,3])
data = parameterscan(p, X, n_reps, analyze_S2_traj(ts))

trajectoryplot(p, X, data, config)

X_vid = ParameterRange(X.sym, X.name, [1, 3])
savevideos(p, X_vid, config)

