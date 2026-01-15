include("definitions.jl")
include("create_interpolations.jl")


p = loadparameters("params.toml")
cache = create_interpolations(p)

etas, rhos, Ks, S2s = cache.data

begin
    fig = Figure()

    inds = etas .< 40.0 
    rho_min = minimum(rhos)
    
    Axis(fig[1,1], title = "ρ(η)", xlabel = "η", ylabel = "ρ")
    lines!(etas[inds], rhos[inds])

    Axis(fig[2,1], title = "S₂(η)", xlabel = "η", ylabel = "S₂")
    lines!(etas[inds], S2s[inds])

    Axis(fig[1,2], title = "ρ ↦ K(η(ρ))", xlabel = "ρ", ylabel = "K")
    lines!(rhos, Ks)

    Axis(fig[2,2], title = "ρ ↦ Dₓ + μ K(η(ρ)) ρₘᵢₙ", xlabel = "ρ", ylabel = "effective diffusion")
    lines!(rhos, @. p.D_x + p.mu * Ks * rho_min)

    save("plots/interpolations_alpha=$(round(1 - p.chi^2,digits=4)).png", fig)

    fig
end

