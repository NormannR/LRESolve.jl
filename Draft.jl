
using BenchmarkTools
# %%
include("./src/LRESolve.jl")
using .LRESolve
# %%
#Modèle p 27
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
Π[2,1] = 1.]
# %%

# function sims_to_uhlig(M::ModelSims)
#
#     n = size(M.Γ₀,1)
#     ns = size(ind_z,1)
#     m = n-ns
#
#     ind_z = find(M.)
#     ind = setdiff(1:n, ind_z)
#
#     F = zeros(m,m)
#     G = zeros(m,m)
#     H = zeros(m,m)
#     L = zeros(m,ns)
#     M = zeros(m,ns)
#     N = zeros(ns,ns)
#
#     F .= M.Γ₀[ind,ind]
#     G .= M.Γ₁[ind,ind]
#     H .= M.Γ₁[ind,ind]
#
#     return ModelUhlig(F,G,H,L,M,N)
# end

# %%
M = ModelSims(Γ₀,Γ₁,C,Ψ,Π)
@btime solve_sims(M::ModelSims)
# %%

# %%

struct ModelUhlig
    F::Array{Float64,2}
    G::Array{Float64,2}
    H::Array{Float64,2}
    L::Array{Float64,2}
    M::Array{Float64,2}
    N::Array{Float64,2}
end

function quadeq(Ψ,Γ,Θ)
    m = size(Ψ,2)

    Ξ = zeros(2*m,2*m)
    Δ = zeros(2*m,2*m)

    Ξ[1:m,1:m] .= Γ
    Ξ[1:m,m+1:end] .= Θ
    Ξ[m+1:end,1:m] .= I(m)

    Δ[1:m,1:m] .= Ψ
    Δ[m+1:end,m+1:end] .= I(m)

    F = schur(Ξ,Δ)
    return F

end
