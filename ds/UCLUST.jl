#!/usr/local/bin/julia
#title           :UCLUST.jl
#description     :CLuster matrix in a UCLUST-like manner
#author          :Carlos Vigil Vásquez
#date            :20220713
#version         :20220713a
#notes           :Requires ArgParse.jl & StatsBase.jl
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
#license         :Permission to copy and modify is granted under the MIT license

using ArgParse
using DelimitedFiles
using LinearAlgebra
using Random
using StatsBase
using ProgressMeter

"""
    UCLUST(M::AbstractMatrix{T}, cutoff::T; initial_seed::Any = nothing, rng::T = 42) where {T<:Float64}

Greedy clustering algorithm based in "Search and clustering orders of magnitude faster than
BLAST"¹.

# Arguments
- `M::AbstractMatrix{Float64}` = Matrix to cluster.
- `cutoff::Float64` = Similarity cutoff to use for clustering.
- `initial_seed::Any` = Index of element to use as initial seed (`nothing` for random initial seed).
- `rng::Float64` = Seed used for matrix shuffling.

# References
1. Edgar, R. C. (2010). Search and clustering orders of magnitude faster than BLAST. 
Bioinformatics, 26(19), 2460–2461. https://doi.org/10.1093/bioinformatics/btq461
"""
function UCLUST(M::AbstractMatrix{T}, cutoff::T; initial_seed::Any = nothing, rng::Int64 = 42) where {T<:Float64}
    # Get indices of matrix
    idx = 1:size(M, 1)
    shuf_idx = shuffle(MersenneTwister(rng), idx)

    # Initialize UCLUST-like algorithm
    if initial_seed ≠ nothing
        filter!(i -> i ≠ initial_seed, shuf_idx)
    else
        initial_seed = pop!(shuf_idx)
    end
    seeds = Dict(initial_seed => [])

    # Iterate over each element
    @showprogress for element in shuf_idx
        clustered = false

        # Check if element belongs to a cluster
        for seed in keys(seeds)
            if M[seed, element] ≥ cutoff
                push!(seeds[seed], element)
                clustered = true
                break
            end
        end

        # Create new cluster is element doesn't belong to any cluster
        if !clustered
            seeds[element] = []
        end
    end

    # Filter out orphan centroids from matrix
    filter!(kv -> !isempty(kv.second), seeds)
    idx2keep = [keys(seeds)...]

    return M[idx2keep, idx2keep]
end

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
        "--seed"
        arg_type = Any
        action = :store_arg
        help = "Initial seed node index"
        required = false
        default = nothing
        "--rng"
        arg_type = Int64
        action = :store_arg
        help = "Matrix shuffling seed"
        required = false
        default = 42
    end

    # Argument parsing
    args = parse_args(args, configs)

    # Run clustering algorithm
    M = readdlm(args["i"], args["delimiter"], Float64)
    if M ≠ M'
        @warn "Matrix not symmetric, may produce unwanted results"
    end

    M′ = UCLUST(M, args["cutoff"]; initial_seed = args["seed"], rng = args["rng"])

    # Save clustered matrix
    writedlm(args["o"], M′, args["delimiter"])
end
main(ARGS)
