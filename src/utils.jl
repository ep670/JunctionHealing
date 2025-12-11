# util functions
using Plots
function setup_fig()
    plot(size=(400, 600),linewidth=3, legend=:topright, guidefontsize=16, tickfontsize=14, legendfontsize=14, margin=3Plots.mm)
end

function fracture_to_t_high(time_to_fracture, t_high, t_low)
    period = t_high + t_low
    n_cycles = floor(time_to_fracture / period)
    remainder = time_to_fracture - n_cycles*period
    t_high_total = n_cycles * t_high + min(remainder, t_high)
    return t_high_total
end