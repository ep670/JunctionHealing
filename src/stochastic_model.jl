using Statistics
using CellAdhesion
using Plots
using StatsBase
using Base.Threads
using LsqFit

print("Number of threads: ", nthreads(), "\n")
if nthreads() == 1
    @warn "You are running with a single thread. For faster performance, run Julia with multiple threads by setting the JULIA_NUM_THREADS environment variable."
end


include("loading_patterns.jl")

const tol = (eps(CellAdhesion.CellAdhesionFloat))^(0.125)

# Function to draw rectangles with specified center, width, and height
function draw_square(x, y, width, height)
  for (xi, yi) in zip(x, y)
      x_coords = [xi - width/2, xi + width/2, xi + width/2, xi - width/2, xi - width/2]
      y_coords = [yi - height/2, yi - height/2, yi + height/2, yi + height/2, yi - height/2]
      plot!(x_coords, y_coords, seriestype=:shape, lw=0, c=:black, label="")
  end
end

import CellAdhesion: runcluster, BondModel, distributeforce!
function runcluster_plot(v::Cluster, force::Vector{Float64}, dt::Float64; max_steps::Integer = 1000, verbose::Bool = false)
    # Arbitrary force history applied to the junction
    n = length(force)
    if max_steps > n
           @warn max_steps<=n "Maximum number of steps exceed force vector length"
         max_steps = n
           print("\n Maximum number of steps = ", max_steps, "\n")
    end
    step = 1
    force = convert(Vector{CellAdhesionFloat},force)
    dt = convert(CellAdhesionFloat, dt)
    bond_states = []
    while (step <= max_steps) && (v.state == true)
        F = force[step]
        setforce!(v, F)
        update!(v, dt)
        step = step + 1
        push!(bond_states, deepcopy(v))
    end
    if verbose == true
        if v.state == false
            print("\nJunction broken")
        elseif step > max_steps
            print("Maximum number of iterations reached")
        end
    end
    
    if verbose
      print("\nPloting bond attachement and detachement")
    end
    p = plot()
    for bond_index in 1:bond_states[1].n
      history = [state.u[bond_index].state for state in bond_states]
      # find timepoints where bond state changes
      bond_break = findall(diff(history) .== -1)
      bond_make = findall(diff(history) .== 1)
      
      # check that last timestep has a broken bond and add a bond break at the end
      if history[end] == false
        push!(bond_make, length(history))
      end

      # check if first timestep has a broken bond, if so add a bond break at time 0
      if history[1] == false
        bond_break = vcat(0, bond_break)
      end

      # check that you have the same number of bond breaks and makes
      if length(bond_break) != length(bond_make)
        println("Warning - number of bond breaks and makes do not match")
      end

      # draw the boxes
      for i in 1:length(bond_break)
        draw_square((bond_break[i]+bond_make[i])/2, bond_index, bond_make[i]-bond_break[i], 1)
      end
    end
    
    return v.state, force[step-1], dt*(step-1), (step-1), p
  end


  function initialise_bonds_state(cluster, n_init)
    n_closed = maximum([1,Int(round(n_init*cluster.n))])
    closed_indices = sample(1:cluster.n, n_closed, replace=false)
    for i = 1:1:cluster.n
        if i in closed_indices
            cluster.u[i].state = true
        else
            cluster.u[i].state = false
        end
    end
    return cluster
  end



function count_bonds(x)
    count = 0
    for i = 1:1:x.n
        if x.u[i].state == true
            count += 1
        end
    end
    return count
end

# version of runcluster that returns how the state of bonds changes over time
function runcluster_with_bond_history(v::Cluster, force::Float64, dt::Float64; max_steps::Integer = 1000, verbose::Bool = false, return_bond_history::Bool = false)
    step = 0
 
    force = convert(CellAdhesionFloat,force)
    dt = convert(CellAdhesionFloat, dt)
 
    bond_history = zeros(max_steps+1)
    # changes zeros to nan
    bond_history .= NaN
    bond_history[1] = count_bonds(v)
    while (step < max_steps) && (bond_history[step+1] > 0)
        step = step + 1
        setforce_v2!(v, force)
        update_v2!(v, dt)
        bond_history[step+1] = count_bonds(v)
    end
 
    if verbose == true
        if v.state == false
            println("Junction broken")
        elseif step > max_steps
            println("Maximum number of iterations reached")
        end
    end
    if return_bond_history == true
        return v.state, force, dt*step, step, bond_history
    else 
        return v.state, force, dt*step, step
    end
end

# version of runcluster that returns how the state of bonds changes over time
function runcluster_with_bond_history(v::Cluster, force::Vector{Float64}, dt::Float64; max_steps::Integer = 1000, verbose::Bool = false, return_bond_history::Bool = false)
    step = 0
 
    force = convert(Vector{CellAdhesionFloat},force)
    dt = convert(CellAdhesionFloat, dt)
 
    bond_history = zeros(max_steps+1)
    # changes zeros to nan
    bond_history .= NaN
    bond_history[1] = count_bonds(v)
    while (step < max_steps) && (bond_history[step+1] > 0)
        step = step + 1
        F = force[step]
        setforce_v2!(v, F)
        update_v2!(v, dt)
        bond_history[step+1] = count_bonds(v)
    end
 
    if verbose == true
        if v.state == false
            println("Junction broken")
        elseif step > max_steps
            println("Maximum number of iterations reached")
        end
    end
    if return_bond_history == true
        return v.state, force, dt*step, step, bond_history
    else 
        return v.state, force, dt*step, step
    end
end


function update_v2!(v::Cluster, dt::CellAdhesionFloat)
    for i = 1:1:v.n 
        k = v.u[i]
        update!(k, dt)  
    end

    #Get the state value for each bond
    interface_v = getfield.(v.u, :state);

    # If the sum of the state values is 0, the junction is broken 
    sum_v = sum(interface_v);
    state = isequal(sum_v,0);

    # Update the state value of the junction
    setfield!(v, :state, !state)
end

function setforce_v2!(v::Cluster{Bond{T}}, F::CellAdhesionFloat) where T <:BondModel
    setfield!(v, :f, F)
    distributeforce!(v)
end

function setforce_v2!(v::Cluster{Bond{T}}) where T <:BondModel
    distributeforce!(v)
end

function setforce_v2!(v::Cluster, F::CellAdhesionFloat)
    setfield!(v, :f, F)
    setforce_v2!(v)
end

function setforce_v2!(v::Cluster)
    distributeforce!(v)
    for i = 1:1:v.n 
        k = v.u[i]
        setforce_v2!(k)  
    end
end

# Find the t* and t_r for a range of Γ_high and Γ_low
function get_fracture_time_stochastic(gamma_high, k_on_0, k_off_0, f_1e, duration, gamma_low=nothing, num_sim=100, n=100, is_return_mean=true)
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
        x = Cluster(n, l, model, :force_global)
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


function get_recovery_time_stochastic(gamma_low, k_on_0, k_off_0, f_1e, duration, is_forth_regime, return_bond_info=false, dt = nothing, num_sim = 100, n=100)    
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
        x = Cluster(n, l, model, :force_global)
        x = initialise_bonds_state(x, n_init)
        _, _, recovery_time[sim], _, bond_memories[sim,:] = runcluster_with_bond_history(x, tension_bond, dt, max_steps = n_timesteps, verbose=false, return_bond_history=true)
    end
    
    # find simulations that did not break
    idx_recovered = bond_memories[:,end] .== bond_memories[:,end]
    pct_recovery = 100*sum(idx_recovered)/num_sim
    # check if there are any simulations that recovered
    if pct_recovery == 0.
        println("No simulations recovered")        
        if return_bond_info == true
            return NaN, 0., bond_memories
        else
            return NaN
        end
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
