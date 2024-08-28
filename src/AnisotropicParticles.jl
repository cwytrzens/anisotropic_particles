# for now, we just include the files manually 

# to install, go into REPL, press ']' and then pkg> add GLMakie StaticArrays
const use_gpu = false

using LinearAlgebra, Random, StaticArrays, TOML, Statistics

using GLMakie, ForwardDiff, ProgressMeter, Accessors, JLD2, ProgressLogging
using GLMakie.GeometryBasics

using DiffEqCallbacks, StochasticDiffEq, OrdinaryDiffEq
using SpatialHashTables, KernelAbstractions
using SpatialHashTables: dist, dist²

using DataInterpolations, Integrals
using KernelDensity, Statistics


@static if use_gpu 
    using CUDA
    const FloatT = Float32 
    const IntT = Int32  
    const IndexVecT = CuVector{IntT}
    const DataArray = CuArray
    const DataVector = CuVector
    const backend = CUDA.CUDABackend()
else
    const FloatT =Float64 
    const IntT = Int64 
    const IndexVecT = Vector{IntT}

    const DataArray = Array
    const DataVector = Vector
    const backend = KernelAbstractions.CPU()
end
const SVec2 = SVector{2, FloatT}


include("parameters.jl")
include("analysis.jl")
include("particle_model.jl")
include("plotting.jl")

include("massdensity_approx.jl")
include("pde_model.jl")


# module AnisotropicParticles 

    # We make a module later...

# end