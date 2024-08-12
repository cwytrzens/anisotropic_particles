include("definitions.jl")

using DelimitedFiles

function generate_data(p, output_base, sim_name)
    Random.seed!(0)  # set random seed
    s = init(p)
    ts, sol = simulate(s, p)

    at(ts, sol, t) = sol[findfirst(ts .>= t)]

    #result_1 = deepcopy( (;ts, sol) )
    #result_2 = deepcopy( (;ts, sol) )

    # sum( norm(x - y) for (x,y) in zip(sol[end-1].X, result_2.sol[end-1].X))

    ####
    # Single plot
    ####
    begin
        Makie.inline!(true)
        fig = Figure(size=(1600, 600))

        ax = Axis(fig[1, 1], aspect=DataAspect(), title="init   t="*string(p.t_start))

        s = at(ts, sol, p.t_start)
        X = Point2f.(s.X)
        angles = mod.(s.theta, π)
        U = Point2f.(sincos.(s.theta))
        E = [ellipse(s.X[i], s.theta[i], p) for i in 1:p.N]
        poly!(ax, E, color = angles, colorrange = (0.0, π), colormap = :cyclic_mygbm_30_95_c78_n256_s25)


        writedlm(joinpath(output_base,"theta_t_start.csv"), s.theta, ',')
        writedlm(joinpath(output_base,"sincostheta_t_start.csv"), U, ',')
        writedlm(joinpath(output_base,"X_t_start.csv"), s.X, ',')

        ax = Axis(fig[1, 2], aspect=DataAspect(), title="middle   t="*string(p.t_end/2))

        

        s = at(ts, sol, p.t_end / 2)
        X = Point2f.(s.X)
        theta_tmp=s.theta
        angles = mod.(s.theta, π)
        U = Point2f.(sincos.(s.theta))
        E = [ellipse(s.X[i], s.theta[i], p) for i in 1:p.N]
        poly!(ax, E, color = angles, colorrange = (0.0, π), colormap = :cyclic_mygbm_30_95_c78_n256_s25)

        writedlm(joinpath(output_base,"theta_t_middle.csv"), s.theta, ',')
        writedlm(joinpath(output_base,"sincostheta_t_middle.csv"), U, ',')
        writedlm(joinpath(output_base,"X_t_middle.csv"), s.X, ',')


        ax = Axis(fig[1, 3], aspect=DataAspect(), title="terminal   t="*string(p.t_end))

        s = at(ts, sol, p.t_end)
        angles = mod.(s.theta, π)
        X = Point2f.(s.X)
        U = Point2f.(sincos.(s.theta))
        E = [ellipse(s.X[i], s.theta[i], p) for i in 1:p.N]
        poly!(ax, E, color = angles, colorrange = (0.0, π), colormap = :cyclic_mygbm_30_95_c78_n256_s25)

        writedlm(joinpath(output_base,"theta_t_end.csv"), s.theta, ',')
        writedlm(joinpath(output_base,"sincostheta_t_end.csv"), U, ',')
        writedlm(joinpath(output_base,"X_t_end.csv"), s.X, ',')


        plot_label=sim_name * ", Du= "*string(p.D_u) *",Dx="*string(p.D_x)*", mu="*string(p.mu)*", lambda="*string(p.lambda)*", N="*string(p.N)*", l="*string(p.l)*", d="*string(p.d)
       # plot_label="IBM simulation with Du= "*string(p.D_u) *",Dx="*string(p.D_x)*", mu="*string(p.mu)*", lambda="*string(p.lambda)*", N="*string(p.N)*", l="*string(p.l)*", d="*string(p.d)

        Label(fig[0,1:3],plot_label, fontsize = 32)

        save(joinpath(output_base, "snapshots.png"), fig)
        fig
    end



    ####
    # Make video 
    ####

    Makie.inline!(false)
    s_obs = Observable(s)
    fig = init_plot(s_obs, p)
    ax = content(fig[1,1])
    record(fig, joinpath(output_base, "movie.mp4"), eachindex(sol)) do i 
        s_obs[] = sol[i]  # update state
        ax.title = sim_name * " t = $(round(ts[i], digits = 2))" 
        #ax.title = "IBM simulation, t = $(round(ts[i], digits = 2))" 

    end

end




input_folder = "input"
output_folder = "output"


mkpath(output_folder)
redo = false 

for (root, folder, files) in walkdir(input_folder)

    for fn in files

        param_fn = joinpath(root, fn)
        sim_name = splitext(fn)[1]
        output_base = joinpath(output_folder, relpath(root, input_folder), sim_name)
        println(output_base)

        if !isdir(output_base) || redo
            
            mkpath(output_base)
            cp(param_fn, joinpath(output_base, "param.toml"), force=true)

            p = loadparameters(param_fn)
            printstyled("generate $output_base, $sim_name\n", color=:green)
            generate_data(p,output_base,sim_name)
        end

    end
end



# p = loadparameters("params.toml")
# output_base = "output\\Du\\params"
# sim_name = "test"

# function sidebysidefig(folder_path)

#     input_folder = "input"
#     output_folder = "output"

#     vary_folder = "varyDuDx"

#     input_base = joinpath(input_folder, vary_folder)

#     mkpath(input_base)
#     mkpath(output_folder)
#     redo = false 
#     counter=0;

#     for (root, folder, files) in walkdir(input_folder)

#         for fn in files

#             param_fn = joinpath(root, fn)
#             sim_name = splitext(fn)[1]
#             output_base = joinpath(output_folder, relpath(root, input_folder), sim_name)
#             println(output_base)


#         Makie.inline!(true)
#         fig = Figure(size=(1600, 600))

#         ax = Axis(fig[1, 1], aspect=DataAspect(), title="init")

#         s = at(ts, sol, p.t_start)
#         X = Point2f.(s.X)
#         U = Point2f.(sincos.(s.theta))
#         E = [ellipse(s.X[i], s.theta[i], p) for i in 1:p.N]
#         poly!(ax, E)
#     end
# end