# to install, go into REPL, press ']' and then pkg> add GLMakie StaticArrays
const use_gpu = false

if use_gpu 
    using CUDA
    const FloatT = Float32 
    const IntT = Int32  
    const IndexVecT = CuVector{IntT}
    const DataArray = CuArray
    const DataVector = CuVector
else
    const FloatT =Float64 
    const IntT = Int64 
    const IndexVecT = Vector{IntT}

    const DataArray = Array
    const DataVector = Vector
end

using LinearAlgebra, Random, StaticArrays
using GLMakie, ForwardDiff, ProgressMeter, Accessors
using GLMakie.GeometryBasics
using SpatialHashTables, KernelAbstractions
using TOML
using CellListMap
using SpatialHashTables: dist, dist²



#= 

A few notes 
- using `global` is not good for performance. try to pass all parameters (can account for 10 to 100 times faster!)
  if globals are really needed, use `const` instead.
  - to avoid globals, one can use NamedTuples to pass many values around 
    define: nt = (a = 1, b = 2)
    use: nt.a, nt.b, 
    unpacl: (;b, a) = nt

- an operation like `x = [R[1] R[2]` creates a new vector and is hence a bit slow, consider using `x = (R[1], R[2])` or `(@view R[1:2])^T`
  - it is better to use StaticArrays instead. For example, 
    inv(::SMatrix) is very fast, but inv(::Matrix) needs to allocated each step 

- I would not use DifferentialEquations.jl here:
  - we don't need to mash X and theta together in one array, can use two separate arrays instead 
  - we might needs some cache variables
  - the time stepping is simple, so, ODEProblems do not really accelerate anything

- for the plot, one can use Makie, which actually also supports ellipses! 
  the animation works like this: 
  - first create a plot, but use Observable and @lift to define how the state 
    translates into the data which needs to be plotted.
  - then, we just need to update the state to record the movie 
  - one could also have a slider in the window for interactive use 
=#  


const SVec2 = SVector{2, FloatT}

const backend = use_gpu ? CUDA.CUDABackend() : KernelAbstractions.CPU()


#############################
# code to define parameters
# #############################
# p = let 
#     l = 1.5
#     d = 0.2 
#     cutoff = l * 2
    
#     (
#         N = 2000,  # total number of particles 
#         l = l,  # length of major axis of ellipse 
#         d = d,  # length of minor axis of ellipse
#         t_step = 1.0,
#         t_save = 0.0,
#         t_start = 0.0,
#         t_end = 5000.0,
#         Lx = 50.0,
#         Ly = 50.0,
#         mu = 100.0,
#         lambda = 1.0,
#         D_x  = 0 * 0.01,
#         D_u =  0 * 0.001,
#         periodic = true,
#         cutoff = cutoff
#     )
# end 


function loadparameters(fn)
    d = TOML.parsefile(fn)

    # postprocess parameters
    d["chi"] = (d["l"]^2 - d["d"]^2) / (d["l"]^2 + d["d"]^2)
    d["sigma"] = d["D_u"] / d["lambda"]

    return NamedTuple{Tuple(Symbol.(keys(d)))}(Tuple(values(d)))
end



# plotting function
function ellipse(X, theta, p)
    (;l, d) = p
    s, c = sincos(theta) 
    R = @SMatrix[ c -s; s c ]
    return Polygon(
        [Point2f(X + R * SVec2( l * cos(t + theta), d * sin(t + theta) )) for t in LinRange(0,2π,20)]
        )
end

@inline function wrap(p, dx; dom = SA[p.Lx, p.Ly], inv_dom = inv.(dom))
    if p.periodic
        return @. dx - dom * round(dx * inv_dom)
    else 
        return dx
    end
end


#############################
#initial random values for X and theta
#############################
function init(p, rng = Random.default_rng())
    (;N, Lx, Ly) = p

    # create initial positions and angles

    X = [SVec2(rand(rng) * Lx, rand(rng) * Ly) for _ in 1:N] 
    theta = [rand(rng, FloatT) * pi for _ in 1:N]

    if hasproperty(p, :init) 
        if p.init == "Havva_1"
            for i in eachindex(X)
                if X[i][1] < Lx/2 
                    theta[i] = 0.0
                else 
                    theta[i] = pi/2.0
                end
            end
        elseif p.init == "Antoine_1"
            for i in eachindex(X)
                if sin(X[i][1]/Lx * p.freq * 2π) < 0 
                    theta[i] = 0.0
                else 
                    theta[i] = pi/2.0
                end
            end
        end
        

    end

    
    s = (;X, theta)

    if use_gpu
        s = (; X = cu(X), theta = cu(theta))
    end

    return s
end 



function potential(x, p)
    (;l, d) = p 

    R = @SVector[x[1], x[2]]
    alpha = x[3]
    beta = x[4]

    sa, ca = sincos(alpha)
    sb, cb = sincos(beta)
    
    gamma_i = (l^2 - d^2) * @SMatrix[ca^2 ca*sa; ca*sa sa^2] + d^2 * I
    gamma_j = (l^2 - d^2) * @SMatrix[cb^2 cb*sb; cb*sb sb^2] + d^2 * I
    Σ = gamma_i + gamma_j
    
    # expo = hasproperty(p, :expo) ? p.expo : FloatT(0.5)
    expo = FloatT(0.5)

    return FloatT(1/(4π)) * det(Σ)^(expo) * exp(-dot(R, inv(Σ) * R))
end

@kernel function particle_interaction!(p, s, grid, dX, dTheta, dom, inv_dom)

    i = @index(Global) 
    Xi = s.X[i]
    dXi = zero(Xi)
    dThetai = FloatT(0.0)

    N = length(s.X)

    for j in neighbours(grid, Xi, p.cutoff)
        Xij = Xi - s.X[j]

        d² = dist²(Xij)
        if 0 < d² < p.cutoff 
            R = wrap(p, Xij; dom, inv_dom)
            alpha = s.theta[i]
            beta  = s.theta[j]

            z = SA[R[1], R[2], alpha, beta]
            dzij = ForwardDiff.gradient(z -> potential(z, p), z)
            
            dXi     += -1/N * SA[dzij[1], dzij[2]]
            dThetai += -1/N * dzij[3]
        end
    end

    dX[i] = dXi
    dTheta[i] = dThetai
end 


#############################
# Solve ODE 
##########################
function simulate(s_init, p, rng = Random.default_rng()) 

    (;N, mu, lambda) = p
    
    p_gpu = map(x -> (x isa IntT || x isa FloatT || x isa Bool) ? x : nothing, p)

    ts = [p.t_start]
    sol = [s_init]

    s = deepcopy(s_init)


    dX = similar(s.X)
    dTheta = similar(s.theta)

    grid = BoundedGrid(p.cutoff, SA[p.Lx, p.Ly], s.X, IndexVecT)
    particle_interaction_kernel = particle_interaction!(get_backend(s.X), 64)


    n_steps = ceil(Int64, (p.t_end - p.t_start) / p.t_step)
    t_last_save = 0.0

    # rename for convenience
    dt = p.t_step 
    t = p.t_start 


    # define function of potential, capture parameters
    pot(z) = potential(z, p)

    dom = SVec2(p.Lx, p.Ly)
    inv_dom = FloatT.(1 ./ dom)
  

    @showprogress for k in 1:n_steps

        updatecells!(grid, s.X)

        particle_interaction_kernel(p_gpu, s, grid, dX, dTheta, dom, inv_dom, ndrange = length(s.X))

        # integrate forces
        for i in 1:N
            s.X[i]     += dt * mu * dX[i]         + sqrt(dt) * sqrt(2 * p.D_x) * randn(rng, SVec2)
            s.theta[i] += dt * lambda * dTheta[i] + sqrt(dt) * sqrt(2 * p.D_u) * randn(rng, FloatT)
        end

        if p.periodic
            for i in 1:N
                s.X[i] = mod.(s.X[i], SVec2(p.Lx, p.Ly))
            end
        end
        t += dt

        # save, if needed 
        if t_last_save >= p.t_save 
            push!(sol, deepcopy(s))
            push!(ts, t)
            t_last_save = 0.0
        else
            t_last_save += dt
        end
    end

    return ts, sol
end


############################################
# Create plot 
############################################
init_plot(s, p, fig_pos = Figure()[1,1]) = init_plot(Observable(s), p, fig)
fli(x) = SA[x[2],-x[1]]
function init_plot(s::Observable, p, fig_pos = Figure()[1,1])
    ax = Axis(fig_pos, aspect = DataAspect())

    X = @lift Point2f.($s.X) 
    angles = @lift mod.($s.theta, π)
    U = @lift Point2f.(sincos.($s.theta))

    E = @lift [ ellipse(fli($s.X[i]), $s.theta[i], p) for i in 1:p.N ]

    poly!(E, color = angles, colorrange = (0.0, π), colormap = :cyclic_mygbm_30_95_c78_n256_s25) # :phase, :twilight,  :cyclic_mygbm_30_95_c78_n256_s25   
    # scatter!(ax, X, markersize = 1.6)
    #arrows!(ax, X, U, lengthscale = 0.2)

    current_figure() 
end
