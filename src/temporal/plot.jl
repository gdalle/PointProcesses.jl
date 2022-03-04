"""
    plot_events(h)

Plot an event history with its marks if those are of a real subtype.
"""
function plot_events(h::History{M}) where {M<:Real}
    scatter(
        h.times,
        h.marks,
        xlim = (h.tmin, h.tmax),
        title = "Event history",
        xlabel = "Time",
        ylabel = "Marks",
        label = nothing,
    )
end

"""
    plot_intensity(pp, h; npoints)

Plot the conditional intensity function of a temporal point process along a given history.
"""
function plot_intensity(
    pp::TemporalPointProcess{M},
    h::History{M};
    npoints = 100,
) where {M}
    curve_times = h.tmin:(duration(h)/npoints):h.tmax
    curve_vals = [ground_intensity(pp, h, t) for t in curve_times]
    history_times = h.times
    history_vals = [ground_intensity(pp, h, t) for t in history_times]
    plot(curve_times, curve_vals, label = "intensity")
    scatter!(history_times, history_vals, label = "events")
    plot!(
        xlim = (h.tmin, h.tmax),
        xlabel = "Time",
        ylabel = "Ground intensity",
        title = "Intensity evolution",
    )
end

"""
    qqplot_interevent_times(h)

Compare the distribution of inter-event times in a given history to the exponential distribution.
"""
function qqplot_interevent_times(h::History)
    plot(
        qqplot(Exponential, diff(vcat(0., event_times(h)))),
        xlabel = "Quantiles of exponential distribution",
        ylabel = "Quantiles of inter-event times",
        title = "Residual analysis",
    )
end

@doc raw"""
    check_residuals(pp, h)

Check whether the point process `pp` is a good fit for history `h` by applying Ogata's time rescaling method: if $(t_i)_i$ is a temporal point process with intensity $\lambda(t)$, then $(\Lambda(t_i))_i$ is a standard temporal Poisson process.
"""
function check_residuals(pp::TemporalPointProcess, h::History)
    Λ(t) = integrated_ground_intensity(pp, h, t)
    h_rescaled = time_change(h, Λ)
    return qqplot_interevent_times(h_rescaled)
end
