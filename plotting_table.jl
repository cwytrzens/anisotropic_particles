using DelimitedFiles
using DataFrames

#reading data and plotting multiple plots
include("definitions.jl")






input_folder = "output"
vary_folder = "varyDuDx"
output_folder = "plots"

output_folder_path=joinpath(input_folder,vary_folder,output_folder)

redo = false 


lambda= 10.0
mu=10.0

# use sortperm

function read_data_to_plot(input_folder,vary_folder)
  local Du
  local Dx
  local stheta
  local sX
  local dataframe

  c_folder=1
  for (root, folder) in walkdir(joinpath(input_folder,vary_folder))
    if c_folder==1

      
      counter=1

      theta=readdlm(joinpath(root,folder[1,1],"theta_t_end.csv"),',')
     
      Du=zeros(length(folder))
      l=zeros(length(folder))
      d=zeros(length(folder))
      Dx=zeros(length(folder))
      stheta=zeros(length(folder),length(theta))
      sX=zeros(length(folder),length(theta),2)
    
      for fn in folder

        sim_fn = joinpath(root, fn)
        p=loadparameters(joinpath(sim_fn,"param.toml"))
     
        Du[counter]=p.D_u
        Dx[counter]=p.D_x
        l[counter]=p.l
        d[counter]=p.d

        

        theta=readdlm(joinpath(sim_fn,"theta_t_end.csv"),',')
        X=readdlm(joinpath(sim_fn,"X_t_end.csv"),',')
      
 
        stheta[counter,:]=theta  
        sX[counter,:,:]=X
        
        counter=counter+1
 
      end

      c_folder +=1
    end
    #counter_folder +=1

    dataframe=(;Du, Dx, sX,stheta,l,d)
    return dataframe
  end
  
end


# plotting function
function ellipse2(X, theta, l, d )
  s, c = sincos(theta) 
  R = @SMatrix[ c -s; s c ]
  return Polygon(
  [Point2f(X + R * SVec2( l * cos(t + theta), d * sin(t + theta) )) for t in LinRange(0,2π,20)]
     )
 end

 function orderparameter(theta)
    N=length(theta)
    angles = mod.(theta, π)

    S= (3.0*(cos.(angles)).^2 .- 1.0)/2.0
    Sout = sum(S)/N 
  return Sout
end


function plots_table(output_folder_path,Du,Dx,X,theta,l,d)
  Makie.inline!(true)
  fig = Figure(size = (3000, 2000))

  #fig=Figure()


  counter_fig_Dx=1
  counter_fig_Du=1
  N=length(theta[1,:])

    

  for index in 1:length(Du)
    #println("fig="*string([counter_fig_Du,counter_fig_Dx])*"  Du Dx="*string([Du[index], Dx[index]]))

    ax = Axis(fig[counter_fig_Du, counter_fig_Dx], aspect=DataAspect(), title=" ")
    angles = mod.(theta[index,:], π)
    
    E=[ellipse2(X[index,i,:], theta[index,i], l[index], d[index]) for i in 1:N]
    poly!(ax, E, color = angles, colorrange = (0.0, π), colormap = :cyclic_mygbm_30_95_c78_n256_s25)
   
    if index<length(Du)
      if Du[index]==Du[index+1]
        counter_fig_Dx=counter_fig_Dx+1
      else
        counter_fig_Du +=1
        counter_fig_Dx=1;
      end   
    else 
        counter_fig_Dx=counter_fig_Dx+1
    end
  
  end

  plot_label="D_u Dx vary"
  # plot_label="IBM simulation with Du= "*string(p.D_u) *",Dx="*string(p.D_x)*", mu="*string(p.mu)*", lambda="*string(p.lambda)*", N="*string(p.N)*", l="*string(p.l)*", d="*string(p.d)

   #Label(fig[0,1:3],plot_label, fontsize = 32)

   save(joinpath(output_folder_path, "parameterscan.png"), fig)
   fig

  

 end


begin
    df=read_data_to_plot(input_folder,vary_folder)

    df_DuDx=DataFrame(Du=df.Du,Dx=df.Dx)
    index=sortperm(df_DuDx, [:Du, :Dx])
    Dusort=df.Du[index];
    Dxsort=df.Dx[index];
    Xsort=df.sX[index,:,:];
    thetasort=df.stheta[index,:];
    lsort=df.l[index]
    dsort=df.d[index]
    N=length(thetasort[1,:])

    mkpath(output_folder_path)
    plots_table(output_folder_path,Dusort,Dxsort,Xsort,thetasort,lsort,dsort)

end


#begin

  fig = Figure()
  ax = Axis(fig[1,1], title = "order parameter", xlabel = "sigma/nu", ylabel = "S")
  limits!(ax, 0.0, 10.0, 0.0, 1.0)


  N=length(Dxsort)
  avg_dir = [orderparameter(thetasort[i,:]) for i in 1:N]
  sigmamu=(Dxsort*lambda)./(Dusort*mu)
  
  index_order=findall(x-> x>1 , sigmamu)
  println(index_order)
  Dusort
  Dxsort
  scatter!(ax, sigmamu,avg_dir)

  #Colorbar(fig[1,2], ls, label = "time")

  save("plots/average_direction_varyDuDx_claudia.png", fig)

  fig
#end

typeof(thetasort[3,:])
orderparameter(thetasort[1,:])

a=[2,4,8]
b=[1,2,3]
a*4
b*3
(a*4)./(b*3)