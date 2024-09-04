using AnisotropicParticles
using AnisotropicParticles: Z_g, g, create_interpolations 

include("plotutils.jl")

datadir = "paper/data"
plotdir = "paper/plots"

p = loadparameters("paper/inputs/base.toml")

X = ParameterRange(:chi, L"\chi", [0.85, 0.9, 0.95, 0.99])
data_K = [
        create_interpolations(updateparameters(p, (X.sym => x,)); eta_range = 40.0, n_interp = 1000)
        for x in X.range
            ]

JLD2.@save "paper/data/massdensity_approx_traj.jld2" data_K X p

with_theme(manuscript_theme) do 
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = L"\eta", ylabel = L"K(\eta)", title = "Examples of mass density function")

    for (x, interp) in zip(X.range, data_K)
        lines!(interp.data.etas, interp.data.Ks, label = L"%$(X.name) = %$(x)")
    end

    hlines!(1 - minimum(X.range)^2, color = :gray, linestyle = :dash, linewidth = 2, label = L"1 - %$(minimum(X.range))^2")
    hlines!(1 - minimum(X.range)^2/2, color = :gray, linestyle = :dot, linewidth = 2, label = L"1 - %$(minimum(X.range))^2/2")

    Legend(fig[1,2], ax)

    save(joinpath(plotdir, "massdensity_approx_traj.png"), fig)
    save(joinpath(plotdir, "massdensity_approx_traj.eps"), fig)
    
    fig
end

using ProgressMeter
include("parameterscan.jl")
data = []

X = ParameterRange(:chi, L"\chi", LinRange(0.01, 0.99, 50))
Y = ParameterRange(:sigma, L"\sigma", 2 .^ LinRange(-8,8,4))

data = Matrix{Any}(undef, length(X.range), length(Y.range))
prog = Progress(length(X.range) * length(Y.range), 1, "Computing...")

vals = [(i,j,chi,sigma) for (i, chi) in enumerate(X.range) for (j, sigma) in enumerate(Y.range)]

Threads.@threads for (i, j, chi, sigma) in vals
    lambda = 1.0
    D_u = lambda * sigma
    interp = create_interpolations(
        updateparameters(p, (:D_u => D_u, :lambda => lambda, :chi => chi)); 
        n_interp = 2000, eta_range = 40 .* (1 - chi^2), deg = 1.0)
    data[i,j] = (;chi, sigma, interp)
    next!(prog)
end


mins = [minimum(d.interp.data.Ks) for d in data]
maxs = [maximum(d.interp.data.Ks) for d in data]

mins = minimum(mins, dims = 2)[:]
maxs = maximum(maxs, dims = 2)[:]


with_theme(manuscript_theme) do
    fig = Figure(size = (420, 320))
    ax = Axis(fig[1,1], xlabel = L"\chi", ylabel = L"K", title = L"\textbf{Range of mass density function } K(\eta)")

    # for d in data
    #     y = d.interp.data.Ks
    #     x = fill(d.chi, length(y))
    #     scatter!(x, y, label = L"K \text{ evaluated values}", color = (:gray, 0.5))
    # end

    band!(X.range, mins, maxs, label = L"\text{Observed } K \text{ values}", color = (:gray, 0.4))

    chis = LinRange(0,1,100)
    lines!(chis, 1 .- chis.^2, label = L"1 - \chi^2", linewidth = 2, linestyle = :dash, color = Cycled(1))
    lines!(chis, 1 .- chis.^2 ./ 2, label = L"1 - \chi^2/n", linewidth = 2, linestyle = :dot, color = Cycled(2))
    Legend(fig[2,1], ax, merge = true, orientation = :horizontal)

    limits!(0,1,0,1.01)
    save("paper/plots/massdensity_approx.eps", fig)
    save("paper/plots/massdensity_approx.png", fig)
    fig
end


using JLD2

JLD2.@save "paper/outputs/massdensity_approx.jld2" data







with_theme(manuscript_theme) do 
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = L"\rho", ylabel = L"K(\eta(\rho))", title = L"\text{Examples of diffusion parameter } K")

    for (x, interp) in zip(X.range, data_K)
        lines!(interp.data.rhos[2:end], interp.data.Ks[2:end], label = L"%$(X.name) = %$(x)")
    end

    #hlines!(1 - minimum(X.range)^2, color = :gray, linestyle = :dash, linewidth = 2, label = L"1 - %$(minimum(X.range))^2")
    #hlines!(1 - minimum(X.range)^2/2, color = :gray, linestyle = :dot, linewidth = 2, label = L"1 - %$(minimum(X.range))^2/2")

    xlims!(0,0.0001)
    Legend(fig[1,2], ax)

    save(joinpath(plotdir, "diffparam_rho_K.png"), fig)
    save(joinpath(plotdir, "diffparam_rho_K.eps"), fig)
    
    fig
end


with_theme(manuscript_theme) do 
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = L"\rho", ylabel = L"\eta(\rho)", title = L"\eta")

    for (x, interp) in zip(X.range, data_K)
        lines!(interp.data.rhos[2:end], interp.data.etas[2:end], label = L"%$(X.name) = %$(x)")
    end

    #hlines!(1 - minimum(X.range)^2, color = :gray, linestyle = :dash, linewidth = 2, label = L"1 - %$(minimum(X.range))^2")
    #hlines!(1 - minimum(X.range)^2/2, color = :gray, linestyle = :dot, linewidth = 2, label = L"1 - %$(minimum(X.range))^2/2")

    # xlims!(0,0.0001)
    Legend(fig[1,2], ax)

    save(joinpath(plotdir, "eta_rho.png"), fig)
    save(joinpath(plotdir, "eta_rho.eps"), fig)
    
    fig
end


X = ParameterRange(:sigma, L"\sigma", 2 .^ LinRange(-8,8,5))

Dus = X.range .* p.lambda 
data_sigma = [
        create_interpolations(updateparameters(p, (:D_u => x,)); eta_range = 40.0, n_interp = 1000)
        for x in Dus
            ]

JLD2.@save "paper/data/massdensity_approx_traj_sigma.jld2" data_sigma X p


with_theme(manuscript_theme) do 
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = L"\eta", ylabel = L"S_2(\eta)", title = L"S_2")
    interp = create_interpolations(p; eta_range = 40.0, n_interp = 1000)
    lines!(interp.data.etas[2:end], interp.data.S2s[2:end])
    
    #hlines!(1 - minimum(X.range)^2, color = :gray, linestyle = :dash, linewidth = 2, label = L"1 - %$(minimum(X.range))^2")
    #hlines!(1 - minimum(X.range)^2/2, color = :gray, linestyle = :dot, linewidth = 2, label = L"1 - %$(minimum(X.range))^2/2")

    # xlims!(0,0.0001)

    save(joinpath(plotdir, "S_2_eta.png"), fig)
    save(joinpath(plotdir, "S_2_eta.eps"), fig)
    
    fig
end
