using Random

function agentSEIR(ll::Vector{Vector{Int}};
            weighted::Bool=false, rand_init::Bool=true, days::Int=50,
            lw::Vector{Vector{Float64}}=[[1.0]],
            rand_infectious_ratio::Float64=0.001,
            init_infectious_nodes::Vector{Int}=[1],
            init_state::Vector{String}=["N"],
            sigma::Float64=0.4, gamma::Float64=0.2, beta::Float64=0.1)

    n = length(ll)
    curr_infected_neighbors = zeros(Int, n)
    curr_state = String["S" for i in 1:n]
    if rand_init
        init_infectious_nodes = rand(1:n, Int(ceil(n*rand_infectious_ratio)))
        for i in init_infectious_nodes
            curr_state[i] = "I"
        end
    elseif length(init_state) == n && init_state[1] != "N"
        for i in 1:n
            curr_state[i] = init_state[i]
        end
        init_infectious_nodes = findall(x -> x == "I", curr_state)
    else
        for i in init_infectious_nodes
            curr_state[i] = "I"
        end
    end
    for i in init_infectious_nodes
        for j in ll[i]
            curr_infected_neighbors[j] += 1
        end
    end
    next_state = copy(curr_state)
    next_infected_neighbors = copy(curr_infected_neighbors)

    degree = Int[length(ll[i]) for i in 1:n]
    max_d = maximum(degree)
    probs = Float64[1.0-(1.0-beta)^i for i in 1:max_d]
    rand01 = Vector{Float64}(undef, n)

    sum_S = Vector{Int}(undef, days+1)
    sum_E = Vector{Int}(undef, days+1)
    sum_I = Vector{Int}(undef, days+1)
    sum_R = Vector{Int}(undef, days+1)
    sum_S[1] = count(s->(s=="S"), curr_state)
    sum_E[1] = count(s->(s=="E"), curr_state)
    sum_I[1] = count(s->(s=="I"), curr_state)
    sum_R[1] = count(s->(s=="R"), curr_state)
    for t in 1:days
        rand!(rand01)
        for i in 1:n
            if curr_state[i] == "S"
                if curr_infected_neighbors[i] == 0
                    continue
                else
                    p = weighted ? get_prob(ll,lw,degree,curr_state,beta,i) :
                                    probs[curr_infected_neighbors[i]]
                    if p > rand01[i]
                        next_state[i] = "E"
                    end
                end
            elseif curr_state[i] == "E"
                if sigma > rand01[i]
                    next_state[i] = "I"
                    for j in ll[i]
                        next_infected_neighbors[j] += 1
                    end
                end
            elseif curr_state[i] == "I"
                if gamma > rand01[i]
                    next_state[i] = "R"
                    for j in ll[i]
                        next_infected_neighbors[j] -= 1
                    end
                end
            end
        end
        for i in 1:n
            curr_state[i] = next_state[i]
            curr_infected_neighbors[i] = next_infected_neighbors[i]
        end
        sum_S[t+1] = count(s->(s=="S"), curr_state)
        sum_E[t+1] = count(s->(s=="E"), curr_state)
        sum_I[t+1] = count(s->(s=="I"), curr_state)
        sum_R[t+1] = count(s->(s=="R"), curr_state)
    end

    return curr_state, sum_S, sum_E, sum_I, sum_R
end

function get_prob(ll, lw, degree, state, beta, i)
    p = 1.0
    for (idx,j) in enumerate(ll[i])
        if state[j] == "I"
            p *= (1.0 - lw[i][idx]*beta)
        end
    end
    return 1.0 - p
end
