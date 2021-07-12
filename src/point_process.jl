import Distributions: logpdf

abstract type PointProcess{M} end

function ground_intensity(pp::PointProcess{M}, history::History{M}, t) where {M}
    error("not implemented")
end

function intensity(pp::PointProcess{M}, history::History{M}, t, m::M) where {M}
    error("not implemented")
end

function mark_distribution(pp::PointProcess{M}, history::History{M}, t) where {M}
    error("not implemented")
end

function ground_intensity_bound(pp::PointProcess{M}, history::History{M}, t) where {M}
    error("not implemented")
end

function ground_intensity_bound_validity(
    pp::PointProcess{M},
    history::History{M},
    t,
) where {M}
    error("not implemented")
end