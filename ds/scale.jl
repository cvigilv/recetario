#!/usr/local/bin/julia
#title           :scale.jl
#description     :Scale matrix
#author          :Carlos Vigil Vásquez
#date            :20220511
#version         :20220511a
#notes           :Requires ArgParse.jl and StatsBase.jl
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
#license         :Permission to copy and modify is granted under the MIT license

using ArgParse
using LinearAlgebra
using DelimitedFiles
using StatsBase

function minmax(A::AbstractMatrix)
    A_copy = deepcopy(A)

    min = minimum(A_copy)
    max = maximum(A_copy)

    scaler(x) = (x - min) / (max - min)
    scaler.(A_copy)
end

function minmax(A::AbstractMatrix; dims::Int64 = 1)
    A_copy = deepcopy(A)
    StatsBase.transform!(fit(UnitRangeTransform, A; dims = dims, unit = true), A_copy)
    return A_copy
end

function main(args)
    # Parse command-line arguments
    parser = ArgParseSettings()
    add_arg_group!(parser, "I/O options:")
    @add_arg_table! parser begin
        "-i"
        arg_type = String
        action = :store_arg
        help = "Matrix file path"
        required = true
        "-o"
        arg_type = String
        action = :store_arg
        help = "Scaled matrix file path"
        required = true
    end

    add_arg_group!(parser, "Scaling options:")
    @add_arg_table! parser begin
        "--scaler"
        arg_type = Function
        action = :store_arg
        help = "Scaling function"
        required = false
        default = minmax
        "--vmin"
        action = :store_arg
        help = "Minimum value after scaling"
        required = false
        default = 0
        "--vmax"
        action = :store_arg
        help = "Maximum value after scaling"
        required = false
        default = 1
        "--column"
        action = :store_true
        help = "Column-wise scaling"
        "--row"
        action = :store_true
        help = "Column-wise scaling"
        "--element"
        action = :store_true
        help = "Element-wise scaling"
    end
    arguments = parse_args(args, parser)
    arguments["column"] = arguments["column"] == arguments["row"] == arguments["element"] == false ? true : arguments["column"]

    # Load and scale matrix
    M = readdlm(arguments["i"])
    if arguments["element"]
        M_scaled = arguments["scaler"](M)
    elseif arguments["column"]
        M_scaled = arguments["scaler"](M; dims = 1)
    elseif arguments["row"]
        M_scaled = arguments["scaler"](M; dims = 2)
    end

    replace!(M_scaled, NaN => 0)

    writedlm(arguments["o"], M_scaled, ' ')

end

main(ARGS)
