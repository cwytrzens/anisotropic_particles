

# main parameter file 
p = loadparameters("paper/figure_1_alignment.toml")

# statistical number of repetitions 
n_rep = 1


mutable struct ParameterRange
    sym::Symbol
    name
    range 
end

y = ParameterRange(:D_u, L"D_u", LinRange(0, 1, 10))


function create_movie(p) 
    
end

function time_vs_alignment_plot(p, prange::ParameterRange)

end