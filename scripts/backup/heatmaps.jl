include("../definitions.jl")

# main parameter file 
p_base = loadparameters("paper/figure_1_alignment.toml")
p_base = rescale(p_base)

# statistical number of repetitions 
n_reps = 1

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




X = ParameterRange(:mu, L"\mu", LinRange(1, 200, 16))
Y = ParameterRange(:lambda, L"\lambda", LinRange(1, 400, 16))

data = Matrix{Any}(undef, length(X.range), length(Y.range))

analyze(sol) = norm(nematic_mean(State(sol, sol.t[end]).theta))

fn = "$(X.sym)_$(Y.sym)_alignment_2"

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
JLD2.@save "$(fn).jld2" data_avg data_std

begin
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = X.name, ylabel = Y.name, title = "Global Order Parameter") 
    hm = heatmap!(X.range, Y.range, data_avg, interpolate = true, colormap = :plasma)
    Colorbar(fig[1,2], hm, label = L"\Vert \Omega \Vert")

    save("$(fn).png", fig)
    fig
end


fn = "mu_lambda_alignment_2"
JLD2.@load "mu_lambda_alignment_2.jld2" data_avg data_std X Y p_base
CairoMakie.activate!()
with_theme(merge(theme_latexfonts(), theme_minimal())) do
    fig = Figure(size=(420, 320))
    ax = Axis(fig[1,1], xlabel = X.name, ylabel = Y.name, title = "Global order parameter (at terminal state)") 
    hm = heatmap!(X.range, Y.range, data_avg, colormap = :plasma)
    Colorbar(fig[1,2], hm, label = L"\Vert \Omega \Vert")

    save("$(fn).png", fig)
    save("$(fn).eps", fig)
    fig
end
