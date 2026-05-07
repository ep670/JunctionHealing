function square_cycle_with_ramp(gamma_1, gamma_2, t_1, t2, t_ramp)
    t_innit = 0 # initialise tension to gamma_1
    gamma_slope = (gamma_1 - gamma_2) / t_ramp
    return t -> begin
        tm = (t+t_innit) % (t_1 + t2 + 2*t_ramp)
        if tm <= t_1
            (s = gamma_1 , sdot = 0)
        elseif tm <= t_1 + t_ramp
            (s = gamma_1 - gamma_slope * (tm - t_1), sdot = -gamma_slope)
        elseif tm <= t_1 + t2 + t_ramp
            (s = gamma_2, sdot = 0)
        elseif tm <= t_1 + t2 + 2*t_ramp
            (s = gamma_2 + gamma_slope * (tm - t_1 - t2 - t_ramp), sdot = gamma_slope)
        end
    end
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

