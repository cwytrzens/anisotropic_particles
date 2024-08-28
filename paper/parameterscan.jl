
mutable struct ParameterRange
    sym::Symbol
    name
    range 
end

modify(p, nt::Tuple) = rescale(merge(p, NamedTuple(nt)))


function parameterscan(p, X::ParameterRange, n_reps, fnc_analyze; changes = ())
    p = modify(p, changes)
    data = Array{Any}(undef, length(X.range))
    @showprogress for (i, x) in enumerate(X.range) 
        p_ = modify(p, (X.sym => x,))
        sols = simulate_ensemble(p_, n_reps)
        data[i] = fnc_analyze(sols)
    end
    return data
end


function parameterscan(p, X::ParameterRange, Y::ParameterRange, n_reps, fnc_analyze; changes = ())
    p = modify(p, changes)
    data = Array{Any}(undef, length(X.range))

    cases = [(x, y) for x in X.range, y in Y.range]

    @showprogress for ((i,j), (x,y)) in cases 
        p_ = modify(p, (X.sym => x, Y.sym => y))
        sols = simulate_ensemble(p, n_reps)
        data[i] = fnc_analyze(sols)
    end
    return data
end

S2(s::State) = nematic_mean(s.theta)

function analyze_S2_traj(ts) 
    return function (sols) 
        return [norm(S2(State(sol, t))) for t in ts, sol in sols]
    end
end

function analyze_S2_fixedtime(t_end)
    return function (sols)
        return [norm(S2(State(sol, t_end))) for sol in sols]
    end
end

