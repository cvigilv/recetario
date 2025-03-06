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
using ProgressMeter
using LinearAlgebra
using NamedArrays
using DelimitedFiles
using .Threads

Base.setdiff(A::AbstractVector{Bool}, B::AbstractVector{Bool}) = @. A & !B

DelimitedFiles.writedlm(io::IO, x::NamedMatrix{T} where {T}) = writedlm(io, vcat(["" names(x, 2)...], hcat(names(x, 1), x)))
DelimitedFiles.writedlm(io::IO, x::NamedMatrix{T} where {T}, delimiter::Char) = writedlm(io, vcat(["" names(x, 2)...], hcat(names(x, 1), x)), delimiter)
DelimitedFiles.writedlm(io::AbstractString, x::NamedMatrix{T} where {T}, delimiter::Char) = writedlm(io, vcat(["" names(x, 2)...], hcat(names(x, 1), x)), delimiter)

"""
```
parse_matrix(M::AbstractMatrix, rows::Bool, cols::Bool; type::Type)
```

Convert matrix with row/column names to NamedMatrix (assumes names are in firt row/column).

# Arguments

- `M::AbstractMatrix` : Matrix to parse
- `rows::Bool` : Matrix has row names (default = false)
- `cols::Bool` : Matrix has column names (default = false)
- `type::Type` : Type of matrix values (default = Any)
"""
function parse_matrix(M::AbstractMatrix, rows::Bool=false, cols::Bool=false; type::Type=Any)
    # Extract values from matrix
    c_idx = rows ? 2 : 1
    r_idx = cols ? 2 : 1
    values = parse.(type, M[r_idx:end, c_idx:end])

    # Extract dimensions names
    row_names = rows ? [i for i in String.(M[r_idx:end, 1])] : ["R#$i" for i in 1:size(values, 1)]
    col_names = cols ? [i for i in String.(M[1, c_idx:end])] : ["C#$i" for i in 1:size(values, 2)]

    namedM = NamedArray(values, (row_names, col_names))
    namedM = namedM[sort(row_names), sort(col_names)]

    return namedM
end

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
    Tanimoto(M::AbstractMatrix)
"""
function Tanimoto(M::AbstractMatrix)
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
	MD = parse_matrix(readdlm(args["MD"], args["delimiter"], String), args["named"], false; type=Bool)
	ND = parse_matrix(readdlm(args["ND"], args["delimiter"], String), args["named"], false; type=Bool)

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
	print(typeof(MN))
    writedlm(args["MN"], MN, args["delimiter"])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
