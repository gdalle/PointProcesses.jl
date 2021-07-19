function Plots.scatter(history::History{M}) where {M}
    scatter(
        history.times,
        history.marks,
        xlim = (history.tmin, history.tmax),
        title = "Event history",
        xlabel = "Time",
        ylabel = "Marks",
        label = nothing,
    )
end
