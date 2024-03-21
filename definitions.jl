# to install, go into REPL, press ']' and then pkg> add GLMakie StaticArrays

using LinearAlgebra, Random, StaticArrays
using GLMakie, ForwardDiff, ProgressMeter, Accessors
using GLMakie.GeometryBasics
using SpatialHashTables
using TOML

#= A few notes 
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


const SVec2 = SVector{2, Float64}




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
    return NamedTuple{Tuple(Symbol.(keys(d)))}(Tuple(values(d)))
end



# plotting function
function ellipse(X, theta, p)
    (;l, d) = p
    s, c = sincos(theta) 
    R = @SMatrix[ c -s; s c ]
    Polygon(
        [Point2f(X + R * SVec2( l * cos(t + theta), d * sin(t + theta) )) for t in LinRange(0,2π,20)]
        )
end

@inline function neighbours_bc(p, st, pos, r)
    if p.periodic
        return periodic_neighbours(st, pos, r)
    else
        return ( (i, zero(SVec2)) for i in neighbours(st, pos, r) )
    end
end

@inline function wrap(p, dx)
    if p.periodic
        dom = SVec2(p.Lx, p.Ly)
        return @. dx - dom * round(dx ./ dom)
    else 
        return dx
    end
end


#############################
#initial random values for X and theta
#############################
function init(p)
    (;N, Lx, Ly) = p

    # create initial positions and angles
    X = [SVec2(rand() * Lx, rand() * Ly) for _ in 1:N] 
    theta = [rand() * pi for _ in 1:N]
    
    s = (;X, theta)

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
    
    return 1/(4π)^1 * sqrt(det(Σ)) * exp(-dot(R, inv(Σ) * R))
end

# begin 
#     i = 1
#     j = 2

#     R = s.X[j] - s.X[i]
#     alpha = s.theta[i]
#     beta  = s.theta[j]

#     x = [R..., alpha, beta]
#     potential(x, p)

#     (dX1, dX2, dalpha, dbeta) = ForwardDiff.gradient( (x) -> potential(x,p), x)
# end 



#############################
# Solve ODE 
##########################
function simulate(s_init, p) 
    (;N, mu, lambda) = p


    ts = [p.t_start]
    sol = [s_init]

    s = deepcopy(s_init)

    dX = similar(s.X)
    dTheta = similar(s.theta)

    bht = BoundedHashTable(s.X, p.cutoff, (0.0, 0.0), (p.Lx, p.Ly))

    n_steps = ceil(Int64, (p.t_end - p.t_start) / p.t_step)
    t_last_save = 0.0

    # rename for convenience
    dt = p.t_step 
    t = p.t_start 

    @showprogress for k in 1:n_steps

        # set forces to zero 
        for i in 1:N
            dX[i] = zero(SVec2)
            dTheta[i] = 0.0
        end
            
        # compute forces
        for i in 1:N
            for (j, offset) in neighbours_bc(p, bht, s.X[i], p.cutoff)
            #for j in 1:N 
                if i != j
                    
                    R = wrap(p, s.X[i] - s.X[j])
                    #R = wrap(p, s.X[i] - s.X[j])

                    if norm(R) < p.cutoff
                        alpha = s.theta[i]
                        beta  = s.theta[j]

                        z = @SVector[R[1], R[2], alpha, beta]
                        dz = ForwardDiff.gradient( (z) -> potential(z,p), z)
                        
                        dX[i]     += -1/N * SVec2(dz[1], dz[2])
                        dTheta[i] += -1/N * dz[3]
                    end
                end
            end
        end

        # integrate forces
        for i in 1:N
            s.X[i]     += dt * mu * dX[i]         + sqrt(dt) * sqrt(2 * p.D_x) * randn(SVec2)
            s.theta[i] += dt * lambda * dTheta[i] + sqrt(dt) * sqrt(2 * p.D_u) * randn()
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
init_plot(s, p, fig = Figure()) = init_plot(Observable(s), p, fig)

function init_plot(s::Observable, p, fig = Figure())
    ax = Axis(fig[1,1], aspect = DataAspect())

    X = @lift Point2f.($s.X) 
    U = @lift Point2f.(sincos.($s.theta))

    E = @lift [ ellipse($s.X[i], $s.theta[i], p) for i in 1:p.N ]

    poly!(E)
    #scatter!(ax, X)
    #arrows!(ax, X, U, lengthscale = 0.2)

    fig 
end
