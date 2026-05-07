# Import the necessary libraries and dependencies
using Plots
using OrdinaryDiffEq
using Suppressor
using Measures
using LsqFit

# function dn/dt 
function dn_dt(n_ratio, p,t)
    f1 = (p.k_on_0) .* (1 .-n_ratio)
    if n_ratio <=0
        return f1
    end
    f2 = exp.(p.Gamma(t).s ./(n_ratio .* p.f_1e)) .* n_ratio .* p.k_off_0
    return (f1 .- f2)
end

function solve_dn_dt(n_ratio, k_on_0, k_off_0, f_1e, Gamma, duration, dt=0.01)
    p = (k_on_0 = k_on_0, k_off_0= k_off_0, f_1e=f_1e, Gamma=Gamma)
    tspan = (0.0, duration)
    prob = ODEProblem(dn_dt, n_ratio, tspan, p, reltol=1e-8, abstol=1e-8, saveat=dt)    # initialisation [s, e, sb, edb; edot] = [0, 0, [0_1,...0_n], [0_1,...0_n], 0.01] for n branches
    sol = solve(prob, Tsit5(), dtmax=dt, maxiters=4e8, verbose=false)   # solve ODE
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
function get_fracture_time(gamma_high, k_on_0, k_off_0, f_1e, duration, gamma_low=nothing)
    loading_creep = constant_loading(gamma_high)
    n_init = k_on_0/(k_on_0+k_off_0)
    if ! isnothing(gamma_low)
        n_stable, n_unstable = critical_n(k_off_0/k_on_0, gamma_low/f_1e, K, F, n_stable_store, n_unstable_store, true)
        n_init = n_stable
    end
    # suppress warnings/console output from solve_dn_dt
    sol_creep = @suppress solve_dn_dt(n_init, k_on_0, k_off_0, f_1e, loading_creep, duration)
    fracture_time = find_fracture(sol_creep)
    return fracture_time
end

function get_recovery_time(gamma_low, k_on_0, k_off_0, f_1e, duration, is_forth_regime, dt = 0.01)
    """
        n_init: initial value of n_ratio if not provided, it will be set to n_unstable + 0.001
    """
    n_stable, n_unstable = critical_n(k_off_0/k_on_0, gamma_low/f_1e, K, F, n_stable_store, n_unstable_store, true)
    if n_unstable != n_unstable
        n_unstable = 0
    end    
    if is_forth_regime
        _, n_critical = critical_f_vector(k_off_0/k_on_0, K, F, n_stable_store, n_unstable_store)
        n_init = n_critical
    else
        n_init = n_unstable + 0.001
    end
    if (n_init < n_unstable) || (n_init > n_stable) 
        println("Warning: n_init should be between $n_unstable and $n_stable")
    end

    loading_recovery = constant_loading(gamma_low)
    sol_recovery = solve_dn_dt(n_init, k_on_0, k_off_0, f_1e, loading_recovery, duration, dt)

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
