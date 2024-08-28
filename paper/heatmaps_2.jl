include("../definitions.jl")

# main parameter file 
p_base = loadparameters("paper/inputs/base.toml")
p_base = rescale(p_base)

# statistical number of repetitions 
n_reps = 5

mutable struct ParameterRange
    sym::Symbol
    name
    range 
end

function modify(p; kwargs...)
    return rescale( (;p..., kwargs...) )
end

function modify(p, nt::Tuple)
    modify(p; NamedTuple(nt)...)
end




X = ParameterRange(:D_x, L"D_x", [0.1, 0.001, 0.0001])
Y = ParameterRange(:D_u, L"D_u", [1e-5, 1e-8, 1e-11, 0.0])

data = Matrix{Any}(undef, length(X.range), length(Y.range))

analyze(sol) = norm(nematic_mean(State(sol, sol.t[end]).theta))

fn = "$(X.sym)_$(Y.sym)_alignment_2"
println("Start simulation for $(fn)")

@progress for (i, x) in enumerate(X.range) 
    for (j, y) in enumerate(Y.range) 

        p_ = modify(p_base, (X.sym => x, Y.sym => y))

        sols = simulate_ensemble(p_, n_reps)

        samples = analyze.(sols)
        
        avg_ = mean(samples)
        std_ = std(samples)

        println(X.sym, ": ", x, ", ", Y.sym, ": ",y, "\toutput: ", avg_)

        data[i, j] = (; avg = avg_, std = std_)
    end
end

data_avg = [d.avg for d in data]
data_std = [d.std for d in data]
JLD2.@save "paper/data/$(fn).jld2" data_avg data_std X Y p_base


# JLD2.@load "paper/data/$(fn).jld2" data_avg data_std X Y p_base

CairoMakie.activate!()

with_theme(merge(theme_latexfonts(), theme_minimal())) do
    fig = Figure(size=(420, 320))
    ax = Axis(fig[1,1], xlabel = X.name, ylabel = Y.name, title = "Global order parameter (at terminal state)") 
    hm = heatmap!(X.range, Y.range, data_avg, colormap = :plasma)
    Colorbar(fig[1,2], hm, label = L"\Vert \Omega \Vert")

    save("paper/plots/$(fn).png", fig)
    save("paper/plots/$(fn).eps", fig)
    fig
end
