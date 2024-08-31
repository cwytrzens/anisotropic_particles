using Test 
using AnisotropicParticles

@test norm(nematic_mean(fill(0.0, 10))) ≈ 1
@test norm(nematic_mean(rand(100_000) .* 2π)) < 1e-2
