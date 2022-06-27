#!/usr/local/bin/julia
#title           :SimilarityMeasures.jl
#description     :Collection of similarity and distance measures
#author          :Carlos Vigil Vásquez
#date            :20220511
#version         :20220511a
#notes           :
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
#license         :Permission to copy and modify is granted under the MIT license

using Base
Base.setdiff(A::AbstractVector{Bool}, B::AbstractVector{Bool}) = @. A & !B

module SimilarityMeasures
export Tanimoto, Euclidean, Tversky

function _compare_matrices(M::AbstractMatrix, N::AbstractMatrix, ƒ::Function)
    # Check if both matrices have the same number of columns
    Mₘ, Mₙ = size(M)
    Nₘ, Nₙ = size(N)
    @assert Mₙ == Nₙ

    # Create measurement matrix
    MN = zeros(Mₘ, Nₘ)
    for i in 1:Mₘ
        for j in 1:Nₘ
            MN[i, j] = ƒ(M[i, :], N[j, :])
        end
    end
    return MN
end
_compare_matrix(M::AbstractMatrix, ƒ::Function) = _compare_matrices(M, M, ƒ)

Tanimoto(X::AbstractVector{Bool}, Y::AbstractVector{Bool}) = sum(X .& Y) / sum(X .| Y)
Tanimoto(X::AbstractVector{<:Number}, Y::AbstractVector{<:Number}) = sum(X .* Y) / (sum(X .^ 2) + sum(Y .^ 2) - sum(X .* Y))
Tanimoto(M::AbstractMatrix) = _compare_matrix(M, Tanimoto)
Tanimoto(M::AbstractMatrix, N::AbstractMatrix) = _compare_matrices(M, N, Tanimoto)

Euclidean(X::AbstractVector{Bool}, Y::AbstractVector{Bool}) = ((sum(X) .- sum(Y)) .^ 2)^0.5
Euclidean(X::AbstractVector{<:Number}, Y::AbstractVector{<:Number}) = (sum(X .- Y)^2)^0.5
Euclidean(M::AbstractMatrix) = _compare_matrix(M, Euclidean)
Euclidean(M::AbstractMatrix, N::AbstractMatrix) = _compare_matrices(M, N, Euclidean)

Tversky(X::AbstractVector{Bool}, Y::AbstractVector{Bool}, α::AbstractFloat, β::AbstractFloat) = sum(X .& Y) / (sum(X .& Y) + α * sum(setdiff(X, Y)) + β * sum(setdiff(Y, X)))
function Tversky(M::AbstractMatrix{Bool}, N::AbstractMatrix{Bool}, α::AbstractFloat, β::AbstractFloat)
    # Check if both matrices have the same number of columns
    Mₘ, Mₙ = size(M)
    Nₘ, Nₙ = size(N)
    @assert Mₙ == Nₙ

    # Create measurement matrix
    MN = zeros(Mₘ, Nₘ)
    for i in 1:Mₘ
        for j in 1:Nₘ
            MN[i, j] = Tversky(M[i, :], N[j, :], α, β)
        end
    end
    return MN
end
Tversky(M::AbstractMatrix{Bool}, α::AbstractFloat, β::AbstractFloat) = Tversky(M, M, α, β)
end
