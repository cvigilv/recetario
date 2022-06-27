#!/usr/local/bin/julia
#title           :TverskyMatrix.jl
#description     :Calculate Tversky index similarity matrix from feature matrix
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

function main(args)
    configs = ArgParseSettings()

    add_arg_group!(configs, "I/O option:")
    @add_arg_table! configs begin
        "--MD", "-i"
        arg_type = String
        action = :store_arg
        help = "M x D matrix"
        required = true
        "--ND"
        arg_type = String
        action = :store_arg
        help = "N x D matrix"
        required = false
        "--MN", "-o"
        arg_type = String
        action = :store_arg
        help = "M x N similarity matrix"
        required = true
        "--alpha"
        arg_type = Float64
        action = :store_arg
        help = "α coefficient"
        default = 0.5
        "--beta"
        arg_type = Float64
        action = :store_arg
        help = "β coefficient"
        default = 0.5
        "--delimiter", "-d"
        arg_type = Char
        action = :store_arg
        help = "Delimiter character between descriptor bits"
        required = false
        default = ' '
    end

    # Argument parsing
    args = parse_args(args, configs)
    args["ND"] = args["ND"] === nothing ? args["MD"] : args["ND"]

    # Load matrices and get number of drugs, compounds, substructures and targets
    println("Calculating Tversky index matrix from file $(args["ND"])")
    MD = readdlm(args["MD"], args["delimiter"], Bool)
    ND = readdlm(args["ND"], args["delimiter"], Bool)
    @assert size(MD, 2) == size(ND, 2)

    # Calculate similarity matrix and save
	MN = Tversky(MD, ND, args["alpha"], args["beta"])
    writedlm(args["MN"], MN, args["delimiter"])
end

main(ARGS)
