function loadparameters(fn)
    d = TOML.parsefile(fn)

    # postprocess parameters
    d["chi"] = (d["l"]^2 - d["d"]^2) / (d["l"]^2 + d["d"]^2)
    d["sigma"] = d["D_u"] / d["lambda"]

    ks = Symbol.(keys(d))

    # convert to GPU suitable types
    vs = map(v -> v isa String ? Symbol(v) : v isa Int ? IntT(v) : v isa Float64 ? FloatT(v) : v, values(d))

    p = NamedTuple(ks .=> vs)
    p = rescale(p)

    return p
end


function rescale(p)
    if hasproperty(p, :density)
        cur_density = pi * p.l * p.d * p.N / (p.Lx * p.Ly)
        factor = sqrt(p.density / cur_density)

        p = @set p.l = p.l * factor 
        p = @set p.d = p.d * factor 
    end

    
    if !hasproperty(p, :cutoff)
        p = (; p..., cutoff = p.cutoff_factor * max(p.l, p.d))
    end

    return p 
end