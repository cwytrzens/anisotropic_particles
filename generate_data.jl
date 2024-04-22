include("definitions.jl")



input_folder="input"
output_folder="output"


mkpath(output_folder)
redo=true

for (root,folder,files) in walkdir(input_folder)
    
    for fn in files

        param_fn=joinpath(root,fn)
        sim_name=splitext(fn)[1]
        output_base=joinpath(output_folder,relpath(root,input_folder),sim_name)
        mkpath(output_base)

        if !isdir(output_base) || redo
            cp(param_fn,joinpath(output_base,"param.toml"),force=true)

            p = loadparameters(param_fn)
            printstyled("generate $output_base, $sim_name\n",color=:green)
            #generate_data(p,output_base,sim_name)
        end
       
    end
end



p=loadparameters("params.toml")
output_base="output\\Du\\params"
sim_name="test"



Random.seed!(0)  # set random seed
s = init(p)
ts, sol = simulate(s, p)

#result_1 = deepcopy( (;ts, sol) )
#result_2 = deepcopy( (;ts, sol) )

# sum( norm(x - y) for (x,y) in zip(sol[end-1].X, result_2.sol[end-1].X))

####
# Single plot
####
# s_obs = Observable(s)
# fig = init_plot(s_obs, p)

####
# Make video 
####
# record(fig, "movie.mkv", eachindex(sol)) do i 
#     s_obs[] = sol[i]  # update state
# end


#####
# Interactive Window 
#####
fig = Figure()
sl = Slider(fig[2,1], range = eachindex(sol))
s_obs = @lift sol[$(sl.value)]

init_plot(s_obs, p, fig[1,1])
fig


#####
# Interactive Window, compare
#####
# fig = Figure()
# sl = Slider(fig[2,1], range = eachindex(sol))
# s_obs = @lift sol[$(sl.value)]

# init_plot(s_obs, p, fig[1,1])

# s_obs_2 = @lift result_1.sol[$(sl.value)]

# init_plot(s_obs_2, p, fig[1,2])
# fig