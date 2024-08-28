# Anisotropic Particles


A particle model with aligment as a result of repulsive forces between anisotropic particles.

## Example

TODO: Insert video.

## Usage

To run the simulations, start a Julia session in the repository folder and activate this project, either in the REPL
via `]` to enter Pkg mode and `pkg> activate` to active the current project and `pkg> instantiate` to install all required dependencies. 

Alternatively, use the Julia code `import Pkg; Pkg.activate(); Pkg.instantiate()`.

### Generating a video

```julia
using AnisotropicParticles

p = loadparameters("parameters.toml")
sol = simulate(p)

createvideo("test.mp4", sol, p)
```

You can also create an interactive window via
```julia
interact(sol, p)
```

### Interacting with the solution

The solution is the return of the `SDEProblem` which solves the problem for `u = [vec(X); theta]` where `X` are the 2D positions of the particles and `theta` are the nematic orientations. You can use the function 
```julia
t = 1.0
s = State(sol, t)
# s.X, s.theta, s.t
```
to extract the data at a given time `t`. 

### Plotting states

Use 
```julia
plotstate(s, p)
```
to plot the current state.

For more access, one can use the function 
```julia
using GLMakie 

s_obs = Observable(State(sol, 0.0))

fig = Figure()
ax = Axis(fig[1,1])

plotstate!(ax, s_obs, p)

display(fig)
```



