using CellAdhesion
using Plots
using Statistics
using NonlinearSolve


function f(u, p)
    return p.k_on .*(1.0 .- u) .- p.k_off_s .* exp.(p.f_vector ./(u.*p.f_1s) ) .* u .- p.k_off_c .* exp.(-p.f_vector ./(u.*p.f_1c) ) .* u
end

function df(u, p)
    return p.k_on .+ p.k_off_s .* exp.(p.f_vector ./(u.*p.f_1s) ) .* (1.0 .- p.f_vector ./ (u.*p.f_1s)) .+ p.k_off_c .* exp.(-p.f_vector ./(u.*p.f_1c) ) .* (1.0 .+ p.f_vector ./ (u.*p.f_1c))
end


# Looks up te stable n solutions for a given k_sim and f_sim
function critical_n(tension, F, n_stable_store, n_unstable_store, return_unstable = false)

    # find the indices of the closest k_vector and f_vector to the given k_sim and f_sim
    f_idx = findmin(abs.(F .- tension))[2]

    # return the critical n for the given k_sim and f_sim
    if return_unstable
        return n_stable_store[f_idx], n_unstable_store[f_idx]
    else
        return n_stable_store[f_idx]
    end
end 



# Given a k_vector, what is the maximum load I can apply before the system becomes unstable
function critical_f_vector(F, n_stable_store, n_unstable_store)
    # find the critical n_value
    last_idx_stable = findlast(n_stable_store .==n_stable_store)
    last_idx_unstable = findlast(n_unstable_store .==n_unstable_store)

    n_critical_idx = max(last_idx_stable, last_idx_unstable)
    
    n_unstable_critical = n_unstable_store[last_idx_unstable]
    n_stable_critical = n_stable_store[last_idx_stable]
    n_critical = (n_unstable_critical + n_stable_critical)/2

    return F[n_critical_idx], n_critical
end


function setup_stores(k_on_0, k_off_s, k_off_c, f_1s, f_1c)
    f_min = 0.0
    f_max = 6.0 * minimum([f_1s, f_1c])
    f_steps = 15000
    f_vectors = collect(range(f_min, f_max, f_steps))
        
    u0 = [0.001,0.1,0.8] # initial guesses for the solver
    
    n_stable_store = zeros(length(f_vectors));
    n_unstable_store = zeros(length(f_vectors));
    
    for (f_idx, f_vector) in enumerate(f_vectors)
        sols = [] #Store all solutions
        p = (k_on=k_on_0, k_off_s=k_off_s, k_off_c=k_off_c, f_vector = f_vector, f_1s=f_1s, f_1c=f_1c)

        # Loop to store all solutions
        for u in u0
            prob = NonlinearProblem(f, u, p);
            sol = solve(prob, NewtonRaphson()); 
            if abs(f(sol.u, p)) < 1e-6
                push!(sols, sol.u)
            end
        end

        # remove duplicates in sols after rounding them to 5 decimal places
        sols = unique(round.(sols; digits = 5))

        #Loop to separate stable and unstable solutions
        for sol in sols
            if sol > 1.0 || sol < 0.0
                continue
            end
            if df(sol, p) > 0.0
                n_stable_store[f_idx] = sol
            else
                n_unstable_store[f_idx] = sol
            end
        end
    end

    F = [f_vectors[j] for j in 1:length(f_vectors)];
    
    n_stable_store[n_stable_store .== 0] .= NaN
    n_unstable_store[n_unstable_store .== 0] .= NaN
        
    return F, n_stable_store, n_unstable_store
end



# ## tests
# k_on_0 = 3.47e-4
# k_off_s = 2.68e-4 
# f_1s = 2.49e-1
# k_off_c = k_off_s * 10
# f_1c = f_1s *0.1

# # test setup stores
# F, n_stable_store, n_unstable_store = setup_stores(k_on_0, k_off_s, k_off_c, f_1s, f_1c)

# # test critical n
# tension = 0.05
# n_stable, n_unstable = critical_n(tension, F, n_stable_store, n_unstable_store, true)
# p = (k_on=k_on_0, k_off_s=k_off_s, k_off_c=k_off_c, f_vector = tension, f_1s=f_1s, f_1c=f_1c)
# plot(0:0.01:1.0, f.(0:0.01:1.0, Ref(p)), label="n_stable_store", color=:green)
# vline!([n_stable], label="n_stable = $(round(n_stable, digits=3))", color=:green)
# vline!([n_unstable], label="n_unstable = $(round(n_unstable, digits=3))", color=:red)
# ylims!(-0.0001, 0.0001)
# u = n_stable
# f(u, p)

# # test critical f_vector
# f_critical, n_critical = critical_f_vector(F, n_stable_store, n_unstable_store)

# plot(n_stable_store, F, label="n_stable_store", color=:green)
# plot!(n_unstable_store, F, label="n_unstable_store", color=:red)
# hline!([f_critical], label="f_critical = $(round(f_critical, digits=3))", color=:blue)
# vline!([n_critical], label="n_critical = $(round(n_critical, digits=3))", color=:black)
# ylabel!("Tension N/m")
# xlabel!("n_ratio")


# prob = NonlinearProblem(f, 0.1, p);
# sol = solve(prob, NewtonRaphson()); 
# sol.u




