abstract type PointProcess{M} end

function intensity(pp::PointProcess{M}, history::History{M}, t, m::M) where {M}
    error("not implemented")
end

function all_marks(pp::PointProcess{M}) where {M}
    error("not implemented")
end

function mark_distribution(pp::PointProcess{M}, history::History{M}, t) where {M}
    error("not implemented")
end

function ground_intensity_bound(pp::PointProcess{M}, history::History{M}, t) where {M}
    error("not implemented")
end

function default_param(::Type{<:PointProcess{M}}, history::History{M}) where {M}
    error("not implemented")
end

## Default implementations

function ground_intensity(pp::PointProcess{M}, history::History{M}, t) where {M}
    return sum(intensity(pp, history, t, m) for m in all_marks(pp))
end