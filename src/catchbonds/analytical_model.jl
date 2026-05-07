using OrdinaryDiffEq
using Suppressor
using LsqFit


# function dn/dt 
function dn_dt_catch(n_ratio, p,t)
    f1 = (p.k_on_0) .* (1 .-n_ratio)
    f2s = exp.(p.Gamma(t).s ./(n_ratio .* p.f_1s)) .* p.k_off_s
    f2c = exp.(-p.Gamma(t).s ./(n_ratio .* p.f_1c)) .* p.k_off_c
    f2 = (f2s .+ f2c) .* n_ratio
    return (f1 .- f2)
end

function dn_dt_catch_separate(n_ratio, p,t)
    f1 = (p.k_on_0) .* (1 .-n_ratio)
    f2s = exp.(p.Gamma(t).s ./(n_ratio .* p.f_1s)) .* p.k_off_s
    f2c = exp.(-p.Gamma(t).s ./(n_ratio .* p.f_1c)) .* p.k_off_c
    f2 = (f2s .+ f2c) .* n_ratio
    return f1, f2
end

function solve_dn_dt_catch(n_ratio, k_on_0, k_off_s, f_1s, k_off_c, f_1c, Gamma, duration, dt=0.01, verbose=false)
    p = (k_on_0 = k_on_0, k_off_s= k_off_s, f_1s=f_1s, k_off_c= k_off_c, f_1c=f_1c, Gamma=Gamma)
    tspan = (0.0, duration)
    prob = ODEProblem(dn_dt_catch, n_ratio, tspan, p, reltol=1e-8, abstol=1e-8, saveat=dt)    # initialisation [s, e, sb, edb; edot] = [0, 0, [0_1,...0_n], [0_1,...0_n], 0.01] for n branches
    sol = solve(prob, Tsit5(), dtmax=dt, maxiters=4e8, verbose=verbose)   # solve ODE
    # unpack sol
    t = sol.t
    n_ratio = [i for i in sol.u]
    return( (t=t, n_ratio=n_ratio) )
end 


function find_fracture(sol)
    # when solution became unstable which only happens at failure
    if (sum(sol.n_ratio.<= 0) == 0)  
        return sol.t[end]
    
    # when simulation reaches <= 0 
    else
        find = findfirst(sol.n_ratio .<= 0)
        return sol.t[find]
    end    
end



# Find the t* and t_r for a range of Γ_high and Γ_low
function get_fracture_time_catch(gamma_high, k_on_0, k_off_s, f_1s, k_off_c, f_1c, duration, gamma_low=nothing, dt=0.01, verbose=true)
    loading_creep = constant_loading(gamma_high)
    n_init = maximum(n_stable_store[n_stable_store.==n_stable_store])    
    if ! isnothing(gamma_low)
        n_stable, n_unstable = critical_n(gamma_low, F, n_stable_store, n_unstable_store, true)
        n_init = n_stable
    end
    # suppress warnings/console output from solve_dn_dt
    sol_creep = solve_dn_dt_catch(n_init, k_on_0, k_off_s, f_1s, k_off_c, f_1c, loading_creep, duration, dt, verbose)
    fracture_time = find_fracture(sol_creep)
    if fracture_time >= duration
        @warn "analytical simulation did not fracture. Consider increasing duration or changing gamma_high/gamma_low"
    end
    return fracture_time
end




function get_recovery_time_catch(gamma_low, k_on_0, k_off_s, f_1s, k_off_c, f_1c, duration, is_forth_regime, dt = 0.01)
    """
        n_init: initial value of n_ratio if not provided, it will be set to n_unstable + 0.001
    """
    n_stable, n_unstable = critical_n(gamma_low, F, n_stable_store, n_unstable_store, true)
    if n_unstable != n_unstable
        n_unstable = 0
    end    
    # n_init = (n_stable + n_unstable) / 2
    # if is_forth_regime
    #     _, n_critical = critical_f_vector(F, n_stable_store, n_unstable_store)
    #     n_init = n_critical
    # else
        n_init = n_unstable + 0.001
    # end
    if (n_init < n_unstable) #|| (n_init > n_stable) 
        println("Warning: n_init should be larger than $n_unstable")
    end

    loading_recovery = constant_loading(gamma_low)
    sol_recovery = solve_dn_dt_catch(n_init, k_on_0, k_off_s, f_1s, k_off_c, f_1c, loading_recovery, duration, dt)

    t = sol_recovery.t
    recovery_curve = - sol_recovery.n_ratio .+ n_stable
    recovery_curve = recovery_curve ./ recovery_curve[1]  
    # fit exponental decay to the recovery curve
    if length(recovery_curve) > 1
        model(x, p) = exp.(-p[1] .* x)
        fit = LsqFit.curve_fit(model, t, recovery_curve, [0.1])
        mean_time_to_recovery = 1 / fit.param[1]
    else
        mean_time_to_recovery = 0
    end
    if mean_time_to_recovery == 0
        println("Simulations did not recover. change n_init/gamma_low or the simulation will likely fail")        
        return NaN
    end
    return mean_time_to_recovery
end




# ###################### debugging plots #########################

# using Plots
# function constant_loading(amp)
#     return t -> (s=amp, sdot=0)
# end


# n_init = 0.5;

# k_on_0 = 3.47e-4
# k_off_s = 2.68e-4 
# f_1s = 2.49e-1
# k_off_c = k_off_s * 2
# f_1c = f_1s *0.1
# # need k_off_c/f_1c > k_off_s/f_1s for catch bond behaviour
# # with k_off c > k_off s, f_1c >~ f_1s
# duration = 10000.0

# sol = solve_dn_dt_catch(n_init, k_on_0, k_off_s, f_1s, k_off_c, f_1c, constant_loading(0.06), duration, 0.01, true)

# plot(sol.t, sol.n_ratio, xlabel="Time", ylabel="Bond Ratio", title="Catch Bond Dynamics under Constant Loading", ylim=(0,1))

# # find Equilibrium n for constant loading Range
# tension = 0.0:0.005:0.15
# n_eqs = Float64[]
# for t in tension
#     sol_eq = solve_dn_dt_catch(n_init, k_on_0, k_off_s, f_1s, k_off_c, f_1c, constant_loading(t), 10000.0, 0.01, true)
#     t, n_ratio = sol_eq.t[end], sol_eq.n_ratio[end]
#     n_eq = t < duration ? 0.0 : n_ratio
#     # if sol_eq.t[end] < 10000.0
#     #     @warn "Solution did not reach full duration"
#     # end
#     push!(n_eqs, n_eq)
# end
# plot(tension, n_eqs, xlabel="Tension", ylabel="Equilibrium Bond Ratio", title="Equilibrium Bond Ratio vs Tension")


# # Plot kon(n) and koff(n)
# n_ratio = 0:0.01:1
# tension = 0.01
# Gamma = constant_loading(tension)
# p = (k_on_0 = k_on_0, k_off_s= k_off_s, f_1s=f_1s, k_off_c= k_off_c, f_1c=f_1c, Gamma=Gamma)
# f1_results=[]
# f2_results=[]

# for n in n_ratio
#     f1,f2 = dn_dt_catch_separate(n, p,0.)
#     push!(f1_results,f1)
#     push!(f2_results,f2)
# end

# plot(n_ratio,f1_results, label="k_on")
# plot!(n_ratio,f2_results, xlabel="Bond ratio", ylabel="dn/dt", ylims=(0,0.001), title= "Tension @ $tension", label="k_off")

    




# ft = get_fracture_time_catch.(0.0:0.001:0.1, k_on_0, k_off_s, f_1s, k_off_c, f_1c, 10000.0, nothing, true)
# plot(0.0:0.001:0.1, ft, xlabel="Gamma_high", ylabel="Fracture Time", title="Fracture Time vs Loading Level for Catch Bonds")


# ht = get_recovery_time_catch.(0.00:0.01:0.5, k_on_0, k_off_s, f_1s, k_off_c, f_1c, 20000.0, false, 0.01)
# plot(0.00:0.01:0.5, ht, xlabel="Gamma_low", ylabel="Recovery Time", title="Recovery Time vs Unloading Level for Catch Bonds")