module LRESolve

    using LinearAlgebra

    export ModelSims, solve_sims

    struct ModelSims
        Γ₀::Array{Float64,2}
        Γ₁::Array{Float64,2}
        C::Vector{Float64}
        Ψ::Array{Float64,2}
        Π::Array{Float64,2}
    end

    function buildΦ(A,B)
        F = svd(A)
        return B*transpose(F.Vt)*Diagonal(1 ./ F.S)*transpose(F.U)
    end

    function solve_sims(M::ModelSims)::Tuple{Array{Float64,1},Array{Float64,2},Array{Float64,2}}

        F = schur(M.Γ₀,M.Γ₁)
        select = abs.(F.α./F.β) .> 1
        ordschur!(F, select)

        ns = count(select)
        n = size(M.Γ₀,1)
        k = n-ns

        Q1 = transpose(F.Q)[1:ns,:]
        Q2 = transpose(F.Q)[ns+1:end,:]
        Φ = buildΦ(Q2*M.Π,Q1*M.Π)

        invΛ11 = inv(F.S[1:ns,1:ns])
        Λ12 = F.S[1:ns,ns+1:end]
        Λ22 = F.S[ns+1:end,ns+1:end]

        Ω11 = F.T[1:ns,1:ns]
        Ω12 = F.T[1:ns,ns+1:end]
        Ω22 = F.T[ns+1:end,ns+1:end]

        H = zeros(n,n)
        H[1:ns,1:ns] .= invΛ11
        H[1:ns,ns+1:end] .= -invΛ11*(Λ12 - Φ*Λ22)
        H[ns+1:end,ns+1:end] .= I(n-ns)
        H = F.Z*H

        Θ₁ = F.Z[:,1:ns]*invΛ11*hcat(Ω11, Ω12 - Φ*Ω22)*transpose(F.Z)
        topΘ = Q1-Φ*Q2
        Θ = H*vcat(topΘ, (Ω22-Λ22)\Q2)*M.C
        Θ₀ = H*vcat(topΘ, zeros(n-ns,n))*M.Ψ

        return (Θ, Θ₀, Θ₁)

    end

end # module
