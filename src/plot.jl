using Plots
import Plots: plot

function plot(history::History)
    scatter(history.times, history.marks, title="Event history", xlabel="Time", ylabel="Marks", label=nothing)
end