using Plots
import Plots: plot

function plot(history::History{M}) where {M}
    scatter(
        get_times(history),
        get_marks(history),
        title = "Event history",
        xlabel = "Time",
        ylabel = "Marks",
        label = nothing,
    )
end
