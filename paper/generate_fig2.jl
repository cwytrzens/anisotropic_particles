
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


using AnisotropicParticles: compute_ellipse_size, ellipse
using StaticArrays
lds = [compute_ellipse_size(p; chi = x, density = y) for x in X.range, y in Y.range] 

lds = Point2f.(Tuple.(lds))
lds = vec(lds)

chis = repeat(X.range, inner = length(Y.range))[:]

ls = LinRange(0.5, 2.5, 100)
ds = LinRange(0, 0.5, 100)
chis = @. (ls^2 - ds'^2) / (ls^2 + ds'^2)

pts = lds
faces = stack(vcat(
    [[i + 16*(j-1), i+1 + 16*(j-1), i+16 + 16*(j-1)] for i in 1:15 for j in 1:15],
    [[i + 1 + 16*(j-1), i+ 17 + 16*(j-1), i+16 + 16*(j-1)] for i in 1:15 for j in 1:15],
))'

begin 
    with_theme(manuscript_theme) do
        fig = Figure()
        ax = Axis(fig[1,1]; xlabel = L"l \text{ (long axis)}", ylabel = L"d \text{ (short axis)}", aspect = 1, title = L"\gamma_f \text{  at } t = %$(p.t_end)")
        rs = 18  # repeat(LinRange(5,20,16), 16)
        # sc = scatter!(lds, color = mean(data, dims=1)[1,:,:][:], colormap = :thermal, markersize = rs, colorrange = (0,1))

        mesh!(pts, faces, color = mean(data, dims=1)[1,:,:][:], colormap = :thermal, colorrange = (0,1), transparency = true)

        contour!(ls, ds, chis, levels = X.range, linewidth = 0.75, linestyle = :dash, color = (:gray, 0.5))


        scale = 60
        lx = (0.5, 3.0)
        ly = (0.1, 0.5)
        ratio = (lx[end]-lx[1]) / (ly[end] - ly[1])
        limits!(ax, lx..., ly...)
        
        for i in [1,8,13,16]
            chi = X.range[i]
            density = Y.range[end] + 0.78
            (l_, d_) = compute_ellipse_size(p; chi = chi, density = density)
            poly!(ellipse(SA[l_, d_], 0, (l=l_/scale*ratio, d=d_/scale)), color = (:black, 0.5))
        end 

        text!(1.15,0.45, text = L"\chi", rotation = pi/4.2)
        for i in [1,3,5,7,9,11,13,15,16]
            (l_, d_) = compute_ellipse_size(p; chi = X.range[i], density = Y.range[end] + 0.12)
            text!(l_, d_, text = L"%$(round(X.range[i],digits=2))", rotation = pi/4.2)
        end  


        Colorbar(fig[1,2], sc, label = L"\gamma_f")
        save(joinpath(plotdir, "havvas_fish.png"), fig)
        save(joinpath(plotdir, "havvas_fish.eps"), fig)
        fig
    end
end