#!/usr/local/bin/julia
#title           :lazyUCLUST.jl
#description     :UCLUST-like clustering without reading to RAM
#author          :Carlos Vigil Vásquez
#date            :20220823
#version         :20220823a
#notes           :Requires ArgParse.jl, StatsBase.jl, NamedArrays.jl & ProgressMeter.jl
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
#license         :Permission to copy and modify is granted under the MIT license

using ArgParse
using DelimitedFiles
using LinearAlgebra
using NamedArrays
using ProgressMeter
using StatsBase

const _TITLE = "lazyUCLUST.jl"
const _DESCRIPTION = "UCLUST-like clustering without reading to RAM"
const _VERSION = "20220823a"

function main(args)
    configs = ArgParseSettings()

    add_arg_group!(configs, "I/O option:")
    @add_arg_table! configs begin
        "-i"
        arg_type = String
        action = :store_arg
        help = "Matrix to cluster"
        required = true
        "-o"
        arg_type = String
        action = :store_arg
        help = "Clustered matrix"
        required = true
        "--delimiter", "-d"
        arg_type = Char
        action = :store_arg
        help = "Delimiter character"
        required = false
        default = ' '
    end

    add_arg_group!(configs, "UCLUST-like parameters:")
    @add_arg_table! configs begin
        "--cutoff"
        arg_type = Float64
        action = :store_arg
        help = "Similarity cutoff"
        required = true
    end

    # Argument parsing
    args = parse_args(args, configs)

    # Run clustering algorithm
    names, seeds = open(args["i"]) do io
        # Extract names of each element in matrix
        names = string.(split(readline(io), args["delimiter"])[begin+1:end])

        # Initialize UCLUST-like algorithm
        initial_seed, _... = string.(split(readline(io), args["delimiter"]))
        seeds = Dict(initial_seed => [])

        # Iterate over every element in matrix
        pbar = Progress(
            length(names);
            desc="Clustering matrix (N seeds = $(length(keys(seeds))))",
            showspeed=true
        )
        i = 1
        while true
            clustered = false
            i += 1

            # Read line and check if it has contents
            contents = readline(io)
            if contents == ""
                break
            end

            # Cluster element
            element, weights... = string.(split(contents, args["delimiter"])[begin:i])
            weights = parse.(Float64, weights)
            for seed in keys(seeds)
                seed_idx = findfirst(x -> x == seed, names)
                if weights[seed_idx] ≥ args["cutoff"]
                    push!(seeds[seed], element)
                    clustered = true
                    break
                end
            end

            # If clustering fails, add new seed
            if !clustered
                seeds[element] = []
            end

            next!(pbar; desc="Clustering matrix (N seeds = $(length(keys(seeds))))")
        end
        ProgressMeter.finish!(
            pbar;
            desc="Clustering complete! (Total seeds = $(length(keys(seeds))))"
        )

        # Filter all empty seeds
        filter!(kv -> !isempty(kv.second), seeds)

        return names, seeds
    end

    seeds = keys(seeds)
    seeds_idx = sort([findfirst(x -> x == seed, names) for seed in seeds])

    # Save clustered matrix
    M = open(args["i"]) do io
        # Create clustered matrix
        M = NamedArray(zeros(length(seeds), length(seeds)))
        setnames!(M, Vector([seeds...]), 1)
        setnames!(M, Vector([seeds...]), 2)
        M = M[sort(string.(seeds)), sort(string.(seeds))]

        # Drop header
        readline(io)

        # Initialize matrix constructor
        i = 0
        pbar = Progress(
            length(names);
            desc="Constructing clustered matrix ($(i)/$(length(names)))",
            showspeed=true
        )

        # Populate clustered matrix
        while true
            i += 1
            i < length(names) ? nothing : break
            contents = readline(io)
            element, weights... = string.(split(contents, args["delimiter"])[begin:i+1])
            if element in seeds
                if contents == ""
                    break
                end

                weights = parse.(Float64, weights[filter(x -> x <= i, seeds_idx)])
                W = zeros(length(seeds))
                for (i, val) in enumerate(weights)
                    W[i] = val
                end
                M[element, :] = W
            end
            next!(pbar; desc="Constructing clustered matrix ($(i)/$(length(names)))")
        end
        ProgressMeter.finish!(pbar; desc="Constructed clustered matrix!")

        # Symmetrize matrix
        M .+= M'
        M[diagind(M)] .= 1

        return M
    end

    # Save clustered matrix
    writedlm(args["o"], M, args["delimiter"])
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("$(_TITLE) v$(_VERSION) - $(_DESCRIPTION)")

    main(ARGS)
end
