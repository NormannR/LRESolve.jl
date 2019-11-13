using LRESolve, Test

@testset "sims" begin
    # Following Structural Macroeconometrics (Dejong and Cheetan) 2nd ed. p. 27
    Γ₀ = [ 0. 1. 0. 0.33 -1. -0.33 -1. ; 0. 1.5 0. -0.06 0.5 0.06 -0.08 ; 1. 0. 0. -0.67 0. -0.33 -1. ; 1. -0.77 -0.23 0. 0. 0. 0. ; 0. 0. 0. 0. 0. 1. 0. ; 0. 0. 0. -0.47 -0.53 0. 0. ; 0. 0. 0. 0. 0. 0. 1.]
    Γ₁ = zeros(7,7)
    Γ₁[1,6] = 0.33
    Γ₁[2,2] = 1.5
    Γ₁[2,5] = 0.5
    Γ₁[5,3] = 0.06
    Γ₁[5,6] = 0.94
    Γ₁[7,7] = 0.9

    C = zeros(7)

    Ψ = zeros(7,1)
    Ψ[7,1] = 1.

    Π = zeros(7,1)
    Π[2,1] = 1.

    M = ModelSims(Γ₀,Γ₁,C,Ψ,Π)
    @inferred solve_sims(M::ModelSims)
end
