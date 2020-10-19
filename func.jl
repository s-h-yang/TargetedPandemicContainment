function get_edges(ll::Vector{Vector{Int}}, n::Int, m::Int)
    edges = Vector{Tuple{Int,Int}}(undef, m)
    ct = 0
    for i = 1:n
        for j in ll[i]
            if j > i
                ct += 1
                edges[ct] = (i,j)
            end
        end
    end
    @assert(ct == m)
    return edges
end

function get_list_of_weights(ll::Vector{Vector{Int}},
            edges::Vector{Tuple{Int,Int}}, edge_betweenness::Vector{Float64},
            edge_to_loc::Dict{Tuple{Int,Int},Int}, target_perc::Float64;
            weight::Float64=0.1)

    lw::Vector{Vector{Float64}} = [ones(Float64,length(l)) for l in ll]
    m = length(edges)
    sorted_indices = sortperm(edge_betweenness, rev=true)
    ct = 0
    for i in sorted_indices
        u,v = edges[i]
        lw[u][edge_to_loc[(u,v)]] = weight
        lw[v][edge_to_loc[(v,u)]] = weight
        ct += 1
        if ct/m > target_perc
            break
        end
    end
    return lw
end

function get_list_of_weights_sorted_edges(ll::Vector{Vector{Int}},
            sorted_edges::Array{Int,2}, edge_to_loc::Dict{Tuple{Int,Int},Int},
            target_perc::Float64; weight::Float64=0.1)

    lw::Vector{Vector{Float64}} = [ones(Float64,length(l)) for l in ll]
    m = size(sorted_edges,1)
    targeted = 0
    for i in 1:m
        u = sorted_edges[i,1]
        v = sorted_edges[i,2]
        lw[u][edge_to_loc[(u,v)]] = weight
        lw[v][edge_to_loc[(v,u)]] = weight
        targeted += 1
        if targeted/m > target_perc
            break
        end
    end
    return lw
end

function get_list_of_weights_degree_dist(ll::Vector{Vector{Int}},
            edge_to_loc::Dict{Tuple{Int,Int},Int}, target_perc::Float64;
            descending::Bool=true, weight::Float64=0.1)

    lw::Vector{Vector{Float64}} = [ones(Float64,length(l)) for l in ll]
    degree = Int[length(ll[i]) for i in 1:length(ll)]
    m = sum(degree)/2
    order = sortperm(degree, rev=descending)
    ct = 0
    for i in order
        for j in ll[i]
            if lw[i][edge_to_loc[(i,j)]] == weight
                continue
            else
                lw[i][edge_to_loc[(i,j)]] = weight
                lw[j][edge_to_loc[(j,i)]] = weight
                ct += 1
            end
            if ct/m > target_perc
                break
            end
        end
        if ct/m > target_perc
            break
        end
    end
    return lw
end

function get_list_of_weights_node_centrality(ll::Vector{Vector{Int}},
            node_centrality::Vector{Float64}, num_edges::Int,
            edge_to_loc::Dict{Tuple{Int,Int},Int}, target_perc::Float64;
            weight::Float64=0.1)

    lw::Vector{Vector{Float64}} = [ones(Float64,length(l)) for l in ll]
    sorted_indices = sortperm(node_centrality, rev=true)
    ct = 0
    for i in sorted_indices
        for j in ll[i]
            if lw[i][edge_to_loc[(i,j)]] == weight
                continue
            else
                lw[i][edge_to_loc[(i,j)]] = weight
                lw[j][edge_to_loc[(j,i)]] = weight
                ct += 1
            end
            if ct/m > target_perc
                break
            end
        end
        if ct/num_edges > target_perc
            break
        end
    end
    return lw
end

function get_list_of_weights_target_nodes(ll::Vector{Vector{Int}},
            edge_to_loc::Dict{Tuple{Int,Int},Int},
            sorted_node_indices::Vector{Int}, target_perc::Float64;
            weight::Float64=0.1)

    lw::Vector{Vector{Float64}} = [ones(Float64,length(l)) for l in ll]
    n = length(ll)
    ct = 0
    for i in sorted_node_indices
        ct += 1
        for j in ll[i]
            lw[i][edge_to_loc[(i,j)]] = weight
            lw[j][edge_to_loc[(j,i)]] = weight
        end
        if ct/n > target_perc
            break
        end
    end
    return lw
end

function get_list_of_weights_degree_dist_target_nodes(
            ll::Vector{Vector{Int}}, edge_to_loc::Dict{Tuple{Int,Int},Int},
            target_perc::Float64; descending::Bool=true, weight::Float64=0.1)

    lw::Vector{Vector{Float64}} = [ones(Float64,length(l)) for l in ll]
    degree = Int[length(ll[i]) for i in 1:length(ll)]
    n = length(degree)
    order = sortperm(degree, rev=descending)
    ct = 0
    for i in order
        ct += 1
        for j in ll[i]
            lw[i][edge_to_loc[(i,j)]] = weight
            lw[j][edge_to_loc[(j,i)]] = weight
        end
        if ct/n > target_perc
            break
        end
    end
    return lw
end
