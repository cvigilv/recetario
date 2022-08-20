#!/usr/local/bin/julia
#title           :TanimotoMatrix.jl
#description     :Calculate Tanimoto coefficient similarity matrix from feature matrix
#author          :Carlos Vigil Vásquez
#date            :20220810
#version         :20220810a
#notes           :Requires ArgParse.jl and ProgressMeter.jl
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
#license         :Permission to copy and modify is granted under the MIT license

using Base
using ArgParse
using DelimitedFiles
using ProgressMeter
using LinearAlgebra
using NamedArrays
using .Threads
include("../julia/named_array_helper.jl")

Base.setdiff(A::AbstractVector{Bool}, B::AbstractVector{Bool}) = @. A & !B

 """
    Tanimoto(X::AbstractVector{Bool}, Y::AbstractVector{Bool})
"""
function Tanimoto(X::T, Y::T) where {T<:AbstractVector{Bool}}
    # Check if both matrices have the same number of columns
    @assert length(X) == length(Y) "Vector have different sizes!"
    return sum(X .& Y) / sum(X .| Y)
end

 """
    Tanimoto(X::AbstractVector{Number}, Y::AbstractVector{Number})
"""
function Tanimoto(X::T, Y::T) where {T<:AbstractVector{Number}}
    # Check if both matrices have the same number of columns
    @assert length(X) == length(Y) "Vector have different sizes!"
    return sum(X .* Y) / (sum(X .^ 2) + sum(Y .^ 2) - sum(X .* Y))
end

"""
    Tanimoto(M::AbstractMatrix{Bool}, N::AbstractMatrix{Bool})
"""
function Tanimoto(M::T, N::T) where {T<:AbstractMatrix}
    # Check if both matrices have the same number of columns
    Mₘ, Mₙ = size(M)
    Nₘ, Nₙ = size(N)
    @assert Mₙ == Nₙ "Matrices should have the same amount of columns!"

    # Create measurement matrix
    MN = zeros(Mₘ, Nₘ)
    pbar = Progress(Mₘ; showspeed=true)
    @threads for i in 1:Mₘ
        for j in 1:Nₘ
            @inbounds MN[i, j] = Tanimoto(M[i, :], N[j, :])
        end
        next!(pbar)
    end

    return MN
end

"""
    Tanimoto(M::T) where {T<:AbstractMatrix}
"""
function Tanimoto(M::T) where {T<:AbstractMatrix}
    # Check if both matrices have the same number of columns
    Mₘ, _ = size(M)

    # Create measurement matrix
    MM = zeros(Mₘ, Mₘ)
    pbar = Progress(Mₘ; showspeed=true)
    @threads for i in 1:Mₘ
        for j in i:Mₘ
            @inbounds MM[i, j] = Tanimoto(M[i, :], M[j, :])
        end
        next!(pbar)
    end
    MM .+= MM'
    MM[diagind(MM)] ./= 2

    return MM
end

function main(args)
    # Argument parsing {{{
    configs = ArgParseSettings()

    add_arg_group!(configs, "I/O option:")
    @add_arg_table! configs begin
        "--MD", "-i"
        arg_type = String
        action = :store_arg
        help = "M × D matrix"
        required = true
        "--ND"
        arg_type = String
        action = :store_arg
        help = "N × D matrix"
        required = false
        "--MN", "-o"
        arg_type = String
        action = :store_arg
        help = "M × N similarity matrix"
        required = true
        "--named"
        action = :store_true
        help = "Matrix has index/names in first column"
        "--delimiter", "-d"
        arg_type = Char
        action = :store_arg
        help = "Delimiter character between bits"
        required = false
        default = ' '
    end

    args = parse_args(args, configs)
    args["ND"] = args["ND"] === nothing ? args["MD"] : args["ND"]
    # }}}

    # Load matrices and get number of drugs, compounds, substructures and targets
    println("Calculating Tversky index matrix from file $(args["ND"])")
    MD = parse_matrix(readdlm(args["MD"], args["delimiter"], String), true, false; type=Bool)
    ND = parse_matrix(readdlm(args["ND"], args["delimiter"], String), true, false; type=Bool)

    # Calculate similarity matrix and save
    if MD == ND
        MN = Tanimoto(MD)
    else
        MN = Tanimoto(MD, ND)
    end
    if args["named"]
        MN = NamedArray(MN)
        setnames!(MN, names(MD, 1), 1)
        setnames!(MN, names(ND, 1), 2)
    end
    writedlm(args["MN"], MN, args["delimiter"])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
