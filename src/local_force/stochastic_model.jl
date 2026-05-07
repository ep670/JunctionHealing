using Statistics
using CellAdhesion
using Plots
using StatsBase
using Base.Threads
using LsqFit



include("../loading_patterns.jl")

# Find the t* and t_r for a range of Γ_high and Γ_low
function get_fracture_time_stochastic_local(gamma_high, k_on_0, k_off_0, f_1e, duration, gamma_low=nothing, num_sim=100, n=100, is_return_mean=true)
    l = 1
    tension_bond = gamma_high*n*l
    n_timesteps = duration*5

    # set n_init to the stable equilibrium of gamma_low if provided, otherwise set to equilibrium at zero tension
    n_init = k_on_0/(k_on_0+k_off_0)
    if gamma_low != nothing
        n_stable, _ = critical_n(k_off_0/k_on_0, gamma_low/f_1e, K, F, n_stable_store, n_unstable_store, true)
        n_init = n_stable
    end
    rupture_time = zeros(num_sim) # storage array for rupture time values
    model = SlipBondModel((k_on_0=k_on_0,), (k_off_0=k_off_0, f_1e=f_1e)) 
    @threads for sim = 1:1:num_sim   # number of simulations
        x = Cluster(n, l, model, :force_local)
        x = initialise_bonds_state(x, n_init)
        ## Run the Montecarlo simulation until it breaks or it reaches the maximum number of iterations
        _, _, rupture_time[sim], _ = runcluster(x, tension_bond, duration/n_timesteps, max_steps = n_timesteps, verbose=false)
    end
    mean_time_to_fracture = mean(rupture_time)
    if mean_time_to_fracture >= duration
        @warn "Stochastic simulation did not fracture. Consider increasing duration or changing gamma_high/gamma_low"
    end
    if is_return_mean
        return mean_time_to_fracture
    else
        return rupture_time
    end
end

# Function to calculate mean along the first dimension, ignoring NaNs
function mean_along_dim1_ignore_nan(A)
    means = [mean(filter(!isnan, A[:, i])) for i in 1:size(A, 2)]
    return means
end

function std_along_dim1_ignore_nan(A)
    stds = [std(filter(!isnan, A[:, i])) for i in 1:size(A, 2)]
    return stds
end


function get_recovery_time_stochastic_local(gamma_low, k_on_0, k_off_0, f_1e, duration, is_forth_regime, return_bond_info=false, dt = nothing, num_sim = 100, n=100)    
    l = 1
    if isnothing(dt)
        n_timesteps = duration*5
        dt = duration/n_timesteps
    else
        n_timesteps = Int(round(duration/dt))
    end

    tension_bond = gamma_low*n*l
    model = SlipBondModel((k_on_0=k_on_0,), (k_off_0=k_off_0, f_1e=f_1e)) 

    n_stable, n_unstable = critical_n(k_off_0/k_on_0, gamma_low/f_1e, K, F, n_stable_store, n_unstable_store, true)
    if is_forth_regime  # default initialisation
        _, n_critical = critical_f_vector(k_off_0/k_on_0, K, F, n_stable_store, n_unstable_store)
        n_init = n_critical
    else
        if n_unstable != n_unstable # check if n_unstable is NaN then set to zero
            n_unstable = 0
        end
        n_init = n_unstable + 0.01  # start slightly above unstable equilibrium
    end
    if (n_init < n_unstable) || (n_init > n_stable)
        println("Warning: n_init should be between $n_unstable and $n_stable")
    end


    idx_missing = 1:num_sim
    recovery_time = zeros(num_sim)
    bond_memories = zeros(num_sim, n_timesteps+1)
    @threads for sim in idx_missing 
        x = Cluster(n, l, model, :force_local)
        x = initialise_bonds_state(x, n_init)
        _, _, recovery_time[sim], _, bond_memories[sim,:] = runcluster_with_bond_history(x, tension_bond, dt, max_steps = n_timesteps, verbose=false, return_bond_history=true)
    end
    
    # find simulations that did not break
    idx_recovered = bond_memories[:,end] .== bond_memories[:,end]
    pct_recovery = 100*sum(idx_recovered)/num_sim
    # check if there are any simulations that recovered
    if pct_recovery == 0.
        println("No simulations recovered")        
        return NaN
    end
    
    # find mean signal of recovered simulations
    mean_signal = mean_along_dim1_ignore_nan(bond_memories)     # bond_mories: [num_sim x n_timesteps]
    std_signal = std_along_dim1_ignore_nan(bond_memories)

    # check if initialisation is too close to the stable number of bonds
    if (n_stable - n_init) < std_signal[1,end]/n
        println("Warning: the initial number of bonds is too close to the stable number of bonds (less than 1 std ($(n_stable-std_signal[1,end]/n)) away)")
    end

    # Fit exponential to the recovery curve to find the recovery time
    t = range(0, stop=duration, length=n_timesteps+1)
    # transform signal to be able to fit an exponential decay
    recovery_curve = - mean_signal .+ n_stable*n
    recovery_curve = recovery_curve ./ recovery_curve[1]
    if length(recovery_curve) > 1
        model_for_fit(x, p) = exp.(-p[1] .* x)
        fit = LsqFit.curve_fit(model_for_fit, t, recovery_curve, [0.1])
        recovery_time = 1 / fit.param[1]
    else
        recovery_time = 0
    end
    
    if return_bond_info == true
        return recovery_time, pct_recovery, bond_memories
    else
        return recovery_time
    end
end
