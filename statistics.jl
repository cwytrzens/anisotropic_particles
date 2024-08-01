using KernelDensity, Statistics

function compute_density(s, p)
    data = reduce(hcat, s.X)'
    return kde(data, boundary = ((0,p.Lx),(0,p.Ly)))
end

function derivation_from_uniform(s, p)
    kd = compute_density(s, p)
    return std(kd.density) / mean(kd.density)
end


# create plot showing uniformness: 

dfn = [derivation_from_uniform(s, p) for s in sol]

begin 
    using GLMakie
    
    fig = Figure()
    ax = Axis(fig[1,1], title = "derivation form uniformness", xlabel = "time", ylabel = "std / mean (in %)")

    lines!(ax, ts, 100 * dfn)
    ylims!(ax, 0, 100)

    save("plots/derivation_from_uniformness.png", fig)

    fig
end