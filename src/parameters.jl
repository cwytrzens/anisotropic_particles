function loadparameters(fn)
    d = TOML.parsefile(fn)

    ks = Symbol.(keys(d))
    # convert to GPU suitable types
    vs = map(v -> v isa String ? Symbol(v) : v isa Int ? IntT(v) : v isa Float64 ? FloatT(v) : v, values(d))

    p = NamedTuple(ks .=> vs)

    # postprocess parameters
    p = updateparameters(p)

    return p
end

function updateparameters(p, changes::Tuple = ())
    p = merge(p, NamedTuple(changes))

    (;l, d) = compute_ellipse_size(p)
    sigma = p.D_u / p.lambda
    nu = p.D_x / p.mu 
    cutoff = p.cutoff_factor * max(l, d)
    return (;p..., l, d, sigma, nu, cutoff)
end

function compute_ellipse_size(p; chi = p.chi, density = p.density)
    (;N, Lx, Ly) = p 

    C = pi*N/Lx/Ly
    l = ( density^2 * (1+chi) / (1-chi) / C^2 )^(1/4)
    d = density / l / C 

    return (;l, d)
end


# function rescaleparameters(p)
#     if hasproperty(p, :density)
#         cur_density = pi * p.l * p.d * p.N / (p.Lx * p.Ly)
#         factor = sqrt(p.density / cur_density)

#         p = @set p.l = p.l * factor 
#         p = @set p.d = p.d * factor 
#     end

#     if !hasproperty(p, :cutoff)
#         p = (; p..., cutoff = p.cutoff_factor * max(p.l, p.d))
#     end

#     return p 
# end


