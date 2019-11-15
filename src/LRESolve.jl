module LRESolve

    using LinearAlgebra

    export ModelSims, solve_sims, ModelUhlig, solve_uhlig

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

    function quadeq(Ψ::Array{Float64,2},Γ::Array{Float64,2},Θ::Array{Float64,2})

        m = size(Ψ,2)
        Ξ = zeros(2*m,2*m)
        Δ = zeros(2*m,2*m)

        Ξ[1:m,1:m] .= Γ
        Ξ[1:m,m+1:end] .= Θ
        Ξ[m+1:end,1:m] .= I(m)

        Δ[1:m,1:m] .= Ψ
        Δ[m+1:end,m+1:end] .= I(m)

        V = eigen(Ξ,Δ)

        perm = sortperm(abs.(V.values))
        λ = V.values[perm[1:m]]

        if any(abs.(λ) .> 1)
            warning("The solution is unstable")
        end

        Λ = diagm(λ)
        Ω = V.vectors[m+1:end,perm[1:m]]
        return Ω*Λ*inv(Ω)

    end

    struct ModelUhlig
        F::Array{Float64,2}
        G::Array{Float64,2}
        H::Array{Float64,2}
        L::Array{Float64,2}
        M::Array{Float64,2}
        N::Array{Float64,2}
    end

    function solve_uhlig(X::ModelUhlig)

        P = quadeq(X.F,-X.G,-X.H)
        V = kron(transpose(X.N),X.F) + kron(I(size(X.N,1)),X.F*P+X.G)
        Q = V \ (-vec(X.L*X.N+X.M))

        return (P,Q)

    end


end # module
