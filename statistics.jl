using KernelDensity, Statistics

#= 
1. Density of positions 

=#

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


#= 
    2. estimation of a vectorfield 

=# 



using Test 
@test norm(nematic_mean(fill(0.0, 10))) ≈ 1
@test norm(nematic_mean(rand(100_000) .* 2π)) < 1e-2



#avg_dir = [nematic_mean(s.theta) for s in sol]


begin 
    using GLMakie
    
    fig = Figure()
    ax = Axis(fig[1,1], title = "average direction", xlabel = "vx", ylabel = "vy")
    limits!(ax, -1.1, 1.1, -1.1, 1.1)

    ls = lines!(ax, getindex.(avg_dir, 1), getindex.(avg_dir, 2), color = ts, colormap = :thermal)

    Colorbar(fig[1,2], ls, label = "time")

    save("plots/average_direction.png", fig)

    fig
end