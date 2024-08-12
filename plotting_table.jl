#reading data and plotting multiple plots
include("definitions.jl")
using DelimitedFiles





input_folder = "output"
vary_folder = "varyDuDx"
output_folder = "plots"

mkpath(output_folder)
redo = false 




# use sortperm


for (root, folder) in walkdir(joinpath(input_folder,vary_folder))

#     Makie.inline!(true)
#     fig = Figure(size=(1600, 600))
#     println(lenght(folder))

#     fig_n=lenght(folder)/2
#     ax = Axis(fig[, ], aspect=DataAspect(), title="init   t="*string(p.t_start))


    for fn in folder

        sim_fn = joinpath(root, fn)

        p=loadparameters(joinpath(sim_fn,"param.toml"))
        push!(Du,p.Du)
        push!(Dx,p.Dx)  

        

        theta=readdlm(joinpath(sim_fn,"theta_t_end.csv"),',')
        push!(stheta, deepcopy(theta))
        
        X=readdlm(joinpath(sim_fn,"X_t_end.csv"),',')
        push!(sX, deepcopy(X))

#         E = [ellipse(X[i], theta[i], p) for i in 1:p.N]
#         poly!(ax, E, color = angles, colorrange = (0.0, π), colormap = :cyclic_mygbm_30_95_c78_n256_s25)

        
        
        
#       #  sim_name = splitext(fn)[1]
#       #  output_base = joinpath(output_folder, relpath(root, input_folder), sim_name)
#       #  println(output_base)

#      #   if !isdir(output_base) || redo
            
#      #       mkpath(output_base)
#      #       cp(param_fn, joinpath(output_base, "param.toml"), force=true)

#    #         p = loadparameters(param_fn)
#     #        printstyled("generate $output_base, $sim_name\n", color=:green)
#    #         generate_data(p,output_base,sim_name)
#    #     end

    end
    dataframe=(Du,Dx,stheta,sX)
end



# readdlm(joinpath(output_base,"theta_t_start.csv"), s.theta, ',')