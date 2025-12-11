using CellAdhesion
using Plots
using Statistics
using NonlinearSolve


function f(u, p)
    return (1.0 .- u) .- p.k_vector .* exp.(p.f_vector ./u ) .* u
end

function df(u, p)
    return 1.0 .+ p.k_vector .* exp.(p.f_vector ./u ) .* (1.0 .- p.f_vector ./ u) 
end


# Looks up te stable n solutions for a given k_sim and f_sim
function critical_n(k_sim, f_sim, K, F, n_stable_store, n_unstable_store, return_unstable = false)

    # find the indices of the closest k_vector and f_vector to the given k_sim and f_sim
    k_idx = findmin(abs.(K[:,1] .- k_sim))[2]
    f_idx = findmin(abs.(F[1,:] .- f_sim))[2]

    # return the critical n for the given k_sim and f_sim
    if return_unstable
        return n_stable_store[k_idx, f_idx], n_unstable_store[k_idx, f_idx]
    else
        return n_stable_store[k_idx, f_idx]
    end
end 



# Given a k_vector, what is the maximum load I can apply before the system becomes unstable
function critical_f_vector(k_sim, K, F, n_stable_store, n_unstable_store)

    # find the indices of the closest k_vector to the given k_sim
    k_idx = findmin(abs.(K[:,1] .- k_sim))[2]

    # find the critical n_value
    n_critical = (minimum(n for n in n_stable_store[k_idx, :] if !isnan(n)) +maximum(n for n in n_unstable_store[k_idx, :] if !isnan(n)))/2

    # find the index of the critical n_value
    n_idx = findmin(abs.(n_stable_store[k_idx, :] .- n_critical))[2]

    # return the critical f_vector for the given k_sim
    return F[1, n_idx], n_critical
end


function setup_stores(k_ratio)
    k_min = k_ratio
    k_max = k_ratio
    k_steps = 1
    k_vectors = collect(range(k_min, k_max, k_steps))
    f_min = 0.0
    f_max = 4.0
    f_steps = 10000
    f_vectors = collect(range(f_min, f_max, f_steps))
        
    u0 = [0.001,0.1,0.8] # initial guesses for the solver
    
    n_stable_store = zeros(length(k_vectors), length(f_vectors));
    n_unstable_store = zeros(length(k_vectors), length(f_vectors));
    for (k_idx, k_vector) in enumerate(k_vectors)
    
        for (f_idx, f_vector) in enumerate(f_vectors)
            sols = [] #Store all solutions
            p = (k_vector=k_vector, f_vector = f_vector) #k_vector = k_off_0/k_on_0
    
            # Loop to store all solutions
            for u in u0
                prob = NonlinearProblem(f, u, p);
                sol = solve(prob, NewtonRaphson()); 
                if abs(f(sol.u, p)) > 0.005 # if the solution is not close to the root, skip it, the dfference between the two curves needs to be max 0.05 to be considered a root
                    continue
                end
                push!(sols, sol.u)
            end
    
            # remove duplicates in sols after rounding them to 5 decimal places
            sols = unique(round.(sols; digits = 5))
    
            #Loop to separate stable and unstable solutions
            for sol in sols
                if sol > 1.0 || sol < 0.0
                    continue
                end
                if df(sol, p) > 0.0
                    n_stable_store[k_idx, f_idx] = sol
                else
                    n_unstable_store[k_idx, f_idx] = sol
                end
            end
        end
    end
    
    
    K = [k_vectors[i] for i in 1:length(k_vectors), j in 1:length(f_vectors)];
    F = [f_vectors[j] for i in 1:length(k_vectors), j in 1:length(f_vectors)];
    
    n_stable_store[n_stable_store .== 0] .= NaN
    n_unstable_store[n_unstable_store .== 0] .= NaN
        
    return K, F, n_stable_store, n_unstable_store
end

