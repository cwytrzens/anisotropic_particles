
mutable struct ParameterRange
    sym::Symbol
    name
    range 
end

modify(p, nt::Tuple{Pair{Symbol,<:Number}}) = rescale(merge(p, NamedTuple(nt)))


function parameterscan(p, X::ParameterRange, n_reps, fnc_analyze)
    data = Array{Any}(undef, length(X.range))
    for (i, x) in enumerate(X.range) 
        sols = simulate_ensemble(p, n_reps)
        data[i] = fnc_analyze(sols)
    end
    return data
end

S2(s::State) = nematic_mean(s.theta)

function analyze_S2_traj(ts) 
    return function (sols) 
        return [S2(State(sol, t)) for t in ts, sol in sols]
    end
end

function analyze_S2_fixedtime(t_end)
    return function (sols)
        return [S2(State(sol, t_end)) for sol in sols]
    end
end

