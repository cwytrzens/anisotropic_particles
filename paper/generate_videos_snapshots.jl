include("generate_base.jl")

n_reps = 3
n_reps_heatmap = 2
p = loadparameters(fn_p)
ts = LinRange(0, p.t_end, 400)
redo = true

config = (;datadir, plotdir, videodir, redo)



using AnisotropicParticles: init_plot, plotparticles!

GLMakie.activate!()
Makie.inline!(false)

vid_theme = Theme(
    fontsize = 24,
    Axis = (
        leftspinevisible = false,
        rightspinevisible = false,
        topspinevisible = false,
        bottomspinevisible = false,
    )
)

set_theme!(merge(vid_theme, theme_black(), theme_latexfonts())) 
set_theme!(merge(vid_theme,manuscript_theme))

using GLMakie.Makie.LaTeXStrings: latexstring


vid_args = (; n_frames = 900 )

changes = (:D_x => 2^-5, :D_u => 2^-11, :chi => 0.92, :mu => 2^8, :lambda => 2^8)
fn = joinpath(config.videodir, "video_$(changes).mp4")
p_, sol = savevideo(fn, p; changes, vid_args...)



changes = (:chi => 0.98,)
changes_str = L"\chi = 0.98,\ \bar{\rho} = 1,\ D_x = 2^{-3},\ D_u = 2^{-11},\ \lambda = 150,\ \mu = 50,\ N = 10^5"
    

changes = ()
changes_str = L"\chi = 0.90,\ \bar{\rho} = 1,\ D_x = 2^{-3},\ D_u = 2^{-11},\ \lambda = 150,\ \mu = 50,\ N = 10^5"
    
changes = (:chi => 0.7,)
changes_str = L"\chi = 0.70,\ \bar{\rho} = 1,\ D_x = 2^{-3},\ D_u = 2^{-11},\ \lambda = 150,\ \mu = 50,\ N = 10^5"
    

all_stuff = [
    (
        changes = (:density => 0.5,), 
        changes_str = L"\chi = 0.90,\ \bar{\rho} = 0.5,\ D_x = 2^{-3},\ D_u = 2^{-11},\ \lambda = 150,\ \mu = 50,\ N = 10^5"
    ),
    (
        changes = (:density => 1.5,), 
        changes_str = L"\chi = 0.90,\ \bar{\rho} = 1.5,\ D_x = 2^{-3},\ D_u = 2^{-11},\ \lambda = 150,\ \mu = 50,\ N = 10^5"
    ),
    (
        changes = (:D_x => 2^-5,), 
        changes_str = L"\chi = 0.90,\ \bar{\rho} = 1.0,\ D_x = 2^{-5},\ D_u = 2^{-11},\ \lambda = 150,\ \mu = 50,\ N = 10^5"
    ),
    (
        changes = (:D_x => 2^-8,), 
        changes_str = L"\chi = 0.90,\ \bar{\rho} = 1.0,\ D_x = 2^{-8},\ D_u = 2^{-11},\ \lambda = 150,\ \mu = 50,\ N = 10^5"
    ),
    (
        changes = (:D_x => 2^-8, :mu => 2^12), 
        changes_str = L"\chi = 0.90,\ \bar{\rho} = 1.0,\ D_x = 2^{-8},\ D_u = 2^{-11},\ \lambda = 150,\ \mu = 2^12,\ N = 10^5"
    ),
    (
        changes = (:D_x => 2^-5, :mu => 2^12), 
        changes_str = L"\chi = 0.90,\ \bar{\rho} = 1.0,\ D_x = 2^{-5},\ D_u = 2^{-11},\ \lambda = 150,\ \mu = 2^12,\ N = 10^5"
    ),
    (
        changes = (:D_x => 2^-5, :mu => 2^14), 
        changes_str = L"\chi = 0.90,\ \bar{\rho} = 1.0,\ D_x = 2^{-5},\ D_u = 2^{-11},\ \lambda = 150,\ \mu = 2^14,\ N = 10^5"
    ),
    (
        changes = (:D_x => 2^-3, :mu => 2^14), 
        changes_str = L"\chi = 0.90,\ \bar{\rho} = 1.0,\ D_x = 2^{-3},\ D_u = 2^{-11},\ \lambda = 150,\ \mu = 2^14,\ N = 10^5"
    ),
    (
        changes = (:D_x => 2^-3, :mu => 2^14, :chi => 0.95), 
        changes_str = L"\chi = 0.95,\ \bar{\rho} = 1.0,\ D_x = 2^{-3},\ D_u = 2^{-11},\ \lambda = 150,\ \mu = 2^14,\ N = 10^5"
    ),
    (
        changes = (:chi => 0.98,), 
        changes_str = L"\chi = 0.98,\ \bar{\rho} = 1.0,\ D_x = 2^{-3},\ D_u = 2^{-11},\ \lambda = 150,\ \mu = 50,\ N = 10^5"
    ),
    (
        changes = (:D_u => 0.0,), 
        changes_str = L"\chi = 0.90,\ \bar{\rho} = 1.0,\ D_x = 2^{-3},\ D_u = 0,\ \lambda = 150,\ \mu = 50,\ N = 10^5"
    )
]

for (changes, changes_str) in all_stuff 
    @show changes_str 
    ts = [0, 10, 100, 1000, 2500, 5000, 7500, 10000]
    p_ = updateparameters(p, changes)
    sol = simulate(p_)
    fn_base = replace("sim_$(changes)", ":" => "", " => " => "=")
    fn_vid = joinpath(config.videodir, fn_base * ".mp4")
    #savevideo(fn_vid, p_; sol) 
    #saveplots(fn, p, sol, ts)  
        
        
    s_obs = Observable(State(sol, 0.0))


    begin
        fig = Figure(size = (768, 768))
        t_str = @lift lpad(round(Int, $s_obs.t), 5, " ")
        gamma_f_str = @lift(string(round(norm(AnisotropicParticles.nematic_mean($s_obs.theta)), digits=2))) 
        title_str = @lift latexstring("t = \\text{", $t_str, "},\\ \\gamma_f = ", $gamma_f_str)


        bg = fig.scene.backgroundcolor[]
        if bg == RGBAf(1,1,1,1)
            fg = :black
            alpha = 0.8
        else
            fg = :white
            alpha = 0.5
        end

        Label(fig[1,1], title_str, padding = (0,0,0,0), height = 32)
        ax = Axis(fig[2,1], xlabel = "", ylabel = "", aspect = DataAspect())

        for ix in -1:1, iy in -1:1
            plotparticles!(ax, p_, s_obs; offset = (ix, iy), alpha = alpha)
        end

        poly!(Rect(0, 0, p.Lx, p.Ly), alpha = 0.0, strokecolor = (fg, 0.8), strokewidth = 3, linestyle = :dash)

        margin = 0.1
        poly!(Rect(-margin*p.Lx, -margin*p.Ly, (margin*p.Lx), (1+2*margin)*p.Ly), 
                color = [(bg,1.0),(bg,0.0),(bg,1.0),(bg,0.0)])

        poly!(Rect(p.Lx, -margin*p.Ly, (margin*p.Lx), (1+2*margin)*p.Ly), 
            color = [(bg,0.0),(bg,1.0),(bg,0.0),(bg,1.0)])


        poly!(Rect(-margin*p.Lx, -margin*p.Ly, (1+2*margin)*p.Lx, margin*p.Ly), 
        color = [(bg,1.0),(bg,1.0),(bg,0.0),(bg,0.0)])

        poly!(Rect(-margin*p.Lx, p.Ly, (1+2*margin)*p.Lx, margin*p.Ly), 
                color = [(bg,0.0),(bg,0.0),(bg,1.0),(bg,1.0)])

        limits!(ax, (-margin * p.Lx, (1+margin) * p.Lx), (-margin * p.Ly, (1 + margin) * p.Ly))
        Label(fig[3,1], changes_str, padding = (0,0,0,0), height = 32)

        display(fig)
    end


    # screenshots 
    mkpath(joinpath(config.videodir, "snapshots"))
    for t in ts 
        s_obs[] = State(sol, t)
        save(joinpath(config.videodir, "snapshots", fn_base * "_$(lpad(t,5,"0")).png"), fig)
    end

end


AnisotropicParticles.init_plot(s_obs, p_, fig[1,1])
display(fig)

    prog = Progress(length(frames), 1)
    record(fig, fn, frames) do t 
        s_obs[] = State(sol, t)
        next!(prog)
    end            
end











using DiffEqCallbacks
GLMakie.activate!()
Makie.inline!(false)

using AnisotropicParticles: initstate

let 
    changes = (:D_x => 2^-5, :D_u => 2^-11, :chi => 0.95, :mu => 50)
    p = loadparameters(fn_p)
    p_ = updateparameters(p, changes)
    with_theme(video_theme) do
        fig = Figure() 
        s_obs = Observable(initstate(p))
        traj = Observable(Point2f[])

        AnisotropicParticles.init_plot(s_obs, p, fig[1,1])
        lines!(traj, color = :orange, linewidth = 2)
        display(fig)

        function update_plot!(inte)
            t = inte.t 
            x, theta = AnisotropicParticles.unpack(inte.u)
            s_obs[] = State(x, theta, t)
            push!(traj[], Point2f(x[1][1], x[1][2]))
            traj[] = traj[]
        end

        cb = PeriodicCallback(update_plot!, 50.0)
        sol = simulate(p_; callbacks = (cb,))
    end
end



let 
    GLMakie.activate!()
    Makie.inline!(false)
    
    p = loadparameters(fn_p)
    xi = 1.0
    changes = (:xi => xi, :lambda => p.lambda / (1 - p.chi^2)^xi )
    p_ = updateparameters(p, changes)
    with_theme(video_theme) do
        fig = Figure() 
        s_obs = Observable(initstate(p))
        traj = Observable(Point2f[])

        AnisotropicParticles.init_plot(s_obs, p, fig[1,1])
        lines!(traj, color = :orange, linewidth = 2)
        display(fig)

        function update_plot!(inte)
            t = inte.t 
            x, theta = AnisotropicParticles.unpack(inte.u)
            s_obs[] = State(x, theta, t)
            push!(traj[], Point2f(x[1][1], x[1][2]))
            traj[] = traj[]
        end

        cb = PeriodicCallback(update_plot!, 50.0)
        sol = simulate(p_; callbacks = (cb,))
    end
end


# p_, sol = savevideo("paper/videos/mu_2^12.mp4", p; changes = (:mu => 2^12, :t_end => 3000,), n_frames = 900)
# p_, sol = savevideo("paper/videos/mu_2^14.mp4", p; changes = (:mu => 2^14, :t_end => 3000,), n_frames = 900)
# p_, sol = savevideo("paper/videos/mu_2^14.mp4", p; changes = (:mu => 2^14, :t_end => 3000,), n_frames = 900)
# p_, sol = savevideo("paper/videos/mu=2^14_Dx=2^-5.mp4", p; changes = (:mu => 2^14, :D_x => 2^-5, :t_end => 3000,), n_frames = 900)
