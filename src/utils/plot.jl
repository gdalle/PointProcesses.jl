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

function plot_intensity(pp::TemporalPointProcess, h::History; npoints = 100)
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

function qqplot_interevent_times(h::History)
    plot(
        qqplot(Exponential, diff(vcat(0.0, h.times))),
        xlabel = "Quantiles of exponential distribution",
        ylabel = "Quantiles of inter-event times",
        title = "Residual analysis",
    )
end
