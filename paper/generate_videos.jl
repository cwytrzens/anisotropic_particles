include("generate_base.jl")

n_reps = 3
n_reps_heatmap = 2
p = loadparameters(fn_p)
ts = LinRange(0, p.t_end, 400)
redo = true

config = (;datadir, plotdir, videodir, redo)

vid_args = (; n_frames = 900 )

changes = (:D_u => 0.0, :D_x => 2^-5)
fn = joinpath(config.videodir, "video_$(changes).mp4")
p_, sol = savevideo(fn, p; changes, vid_args...)


using DiffEqCallbacks
GLMakie.activate!()
Makie.inline!(false)

using AnisotropicParticles: initstate

let 
    changes = (:D_x => 2^-4, :D_u => 2^-11, :chi => 0.95, :mu => 50)
    p = loadparameters(fn_p)
    p_ = updateparameters(p, changes)
    with_theme(video_theme) do
        fig = Figure() 
        s_obs = Observable(initstate(p))
        traj = Observable(Point2f[])

        AnisotropicParticles.init_plot(s_obs, p, fig[1,1])
        lines!(traj, color = :white, linewidth = 3)
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
