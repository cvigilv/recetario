#!/usr/local/bin/julia
#title           :double_edge_swap.jl
#description     :Double edge swap function for graph randomization
#author          :Carlos Vigil Vásquez
#date            :20211121
#version         :20211121a
#notes           :Requires StatsBase
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
#license         :Permission to copy and modify is granted under the MIT license

using StatsBase
using Random

k(vᵢ::Integer, G::AbstractMatrix) = count(!iszero, G[vᵢ, :])
k(eᵢ::AbstractVector) = count(!iszero, eᵢ)
k(G::AbstractMatrix) = mapslices(k, G; dims = 2)

function double_edge_swap(A::AbstractMatrix; nswap = 1, max_tries = 100)
    A′ = deepcopy(A)

    swapcount = 0
    ntries = 0

    degrees = k(A′)
    vertices = collect(1:size(A′, 1))

    while swapcount < nswap
        u, x = sample(vertices, 2)
        if u == x
            continue
        end

        E_v = findall(!iszero, A′[u, :])
        E_y = findall(!iszero, A′[x, :])
        if (length(E_v) == 0) || (length(E_y) == 0)
            continue
        end

        v = sample(E_v)
        y = sample(E_y)
        if v == y
            continue
        end

        if (x ∉ findall(!iszero, A′[u, :])) && (y ∉ findall(!iszero, A′[v, :]))
            uv = A′[u, v]
            xy = A′[x, y]
            ux = A′[u, x]
            vy = A′[v, y]

            A′[u, x] = uv
            A′[v, y] = xy
            A′[x, u] = uv
            A′[y, v] = xy

            A′[u, v] = ux
            A′[x, y] = vy
            A′[v, u] = ux
            A′[y, x] = vy

            swapcount += 1
        end
        if (max_tries > 0) && (ntries >= max_tries)
            throw(ErrorException("Reached maximum number of tries for edge swapping"))

            break
        end
        ntries += 1
    end
    return A′
end
