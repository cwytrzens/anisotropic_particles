using TOML 

#include("definitions.jl")

 #p = TOML.parsefile("params.toml")



data= (l = 1.5,
d = 0.3,
cutoff = 15.0,    
N = 10000,
t_step = 5.0,
t_save = 0.0,
t_start = 0.0,
t_end = 5000.0,
Lx = 100.0,
Ly = 100.0,
mu = 10.0,
lambda = 50.0,
D_x  = 0.01,
D_u =  0.001,
periodic = true,
t_write_to_file=1000.0
);

counter=0;

input_folder = "input"

vary_folder = "varyDuDx"

input_base = joinpath(input_folder, vary_folder)

mkpath(input_base)
redo = true 

for i in 0.000:0.001:0.01
    for j in 0.000:0.001:0.01
        counter=counter+1
        data = @set data.D_u = i
        data = @set data.D_x = j
        #println(data.D_u)
        sim_name="params$counter"
        input_file = joinpath(input_folder, vary_folder,sim_name)
        println(input_base)

        if redo #|| !isfile(sim_name)
            printstyled("generate parameterfile $sim_name\n", color=:green)
            datadict=Dict(pairs(data))
            open(input_file, "w") do io
                TOML.print(io, datadict, sorted=true)
            end
        end
    end
end