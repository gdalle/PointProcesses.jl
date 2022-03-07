"""
    plot_events(h)

Plot an event history with its marks if those are of a real subtype.
"""
function plot_events(h::History{M}) where {M<:Real}
    return scatterplot(
        h.times,
        h.marks;
        xlim=(h.tmin, h.tmax),
        title="Event history",
        xlabel="Time",
        ylabel="Marks",
        label=nothing,
    )
end

"""
    plot_intensity(pp, h; npoints)

Plot the conditional intensity function of a temporal point process along a given history.
"""
function plot_intensity(pp::TemporalPointProcess{M}, h::History{M}; npoints=100) where {M}
    curve_times = (h.tmin):(duration(h) / npoints):(h.tmax)
    curve_vals = [ground_intensity(pp, h, t) for t in curve_times]
    history_times = h.times
    history_vals = [ground_intensity(pp, h, t) for t in history_times]
    plt = lineplot(
        curve_times,
        curve_vals;
        label="intensity",
        xlim=(h.tmin, h.tmax),
        xlabel="Time",
        ylabel="Ground intensity",
        title="Intensity evolution",
    )
    scatterplot!(plt, history_times, history_vals; label="events")
    return plt
end
