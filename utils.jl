using PyCall
using SparseArrays

function scipyCSC_to_julia(A)
    m, n = A.shape
    colPtr = Int[i+1 for i in PyArray(A."indptr")]
    rowVal = Int[i+1 for i in PyArray(A."indices")]
    nzVal = Vector{Float64}(PyArray(A."data"))
    B = SparseMatrixCSC{Float64,Int}(m, n, colPtr, rowVal, nzVal)
    return PyCall.pyjlwrap_new(B)
end

function edgelist_to_adjacency(edgelist::Array{Int,2}; ispycall::Bool=false)
    Is = [edgelist[:,1]; edgelist[:,2]]
    Js = [edgelist[:,2]; edgelist[:,1]]
    A = sparse(Is, Js, ones(Int, length(Is)))
    return ispycall ? PyCall.pyjlwrap_new(A) : A
end
