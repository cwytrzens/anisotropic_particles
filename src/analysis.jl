function nematic_mean(vs::AbstractVector{<:SVector})
    dim = length(eltype(vs))
    avg_tensor = mean( w -> w*w', vs) - 1/dim * I 
    F = eigen(avg_tensor)
    return dim/(dim-1) * F.vectors[:,end] * F.values[end]
end

function nematic_mean(thetas::AbstractVector{Float64})
    vs = map(theta -> SA[cos(theta), sin(theta)], thetas)
    return nematic_mean(vs)
end



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
