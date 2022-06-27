#!/usr/local/bin/julia
#title           :TanimotoMatrix.jl
#description     :Calculate Tanimoto coefficient similarity matrix from feature matrix
#author          :Carlos Vigil Vásquez
#date            :20220511
#version         :20220511a
#notes           :Requires ArgParse.jl and ProgressMeter.jl
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
#license         :Permission to copy and modify is granted under the MIT license

using Base
using ArgParse
using DelimitedFiles
using ProgressMeter

Base.setdiff(A::AbstractVector{Bool}, B::AbstractVector{Bool}) = @. A & !B

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

function main(args)
    configs = ArgParseSettings()

    add_arg_group!(configs, "I/O option:")
    @add_arg_table! configs begin
        "--MD", "-i"
        arg_type = String
        action = :store_arg
        help = "M x F feature matrix file path"
        required = true
        "--ND"
        arg_type = String
        action = :store_arg
        help = "N x F feature matrix file path"
        required = false
        "--MN", "-o"
        arg_type = String
        action = :store_arg
        help = "M x N similarity matrix file path"
        required = true
        "--delimiter", "-d"
        arg_type = Char
        action = :store_arg
        help = "Feature matrix delimiter"
        required = false
        default = ' '
        "--binary"
        action = :store_true
        help = "Calculate similarity assuming binary feature"
    end

    # Argument parsing
    args = parse_args(args, configs)
    args["ND"] = args["ND"] === nothing ? args["MD"] : args["ND"]
    descriptor_type = args["binary"] ? Bool : Float64

    # Load matrices and get number of rows and columns per feature matrix
    println("Calculating Tanimoto coefficient matrix from file $(args["ND"])")
    MD = readdlm(args["MD"], args["delimiter"], descriptor_type)
    ND = readdlm(args["ND"], args["delimiter"], descriptor_type)
    @assert size(MD, 2) == size(ND, 2)

    # Calculate similarity matrix and save
    MN = Tanimoto(MD, ND)
    writedlm(args["MN"], MN, args["delimiter"])
end

main(ARGS)
