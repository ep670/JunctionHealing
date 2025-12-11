function constant_loading(amp)
    return t -> (s=amp, sdot=0)
end

function square_cycle_uneven_new(gamma_1, gamma_2, t_1, t2)
    t_innit = 0 # initialise tension to gamma_1
    return t -> begin
        tm = (t+t_innit) % (t_1 + t2)
        if tm <= t_1
            (s = gamma_1, sdot = 0)
        elseif tm > t_1
            (s = gamma_2, sdot = 0)
        end
    end
end


function loading_stochastic_model(loading, n, duration, dt=1.0)
    l = 1.0 # Distance between bonds (unit distance)
    interface_length = n*l # (m)

    n_timesteps =  round(Int, duration/dt) - 1
    time = collect(LinRange(0, duration, n_timesteps))
    loading_result = loading.(time)
    tension = [loading_res.s for loading_res in loading_result]
    force = tension .* (interface_length)    # force = tension * meters
    return force, n_timesteps
end