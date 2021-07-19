abstract type PointProcess{M} end

function get_θ(pp::PointProcess)
    return ComponentVector{Float64}(ntfromstruct(pp))
end

function intensity(pp::PointProcess{M}, history::History{M}, t, m::M) where {M}
    return intensity(typeof(pp), get_θ(pp), history, t, m)
end

function mark_distribution(pp::PointProcess{M}, history::History{M}, t) where {M}
    return mark_distribution(typeof(pp), get_θ(pp), history, t)
end

function ground_intensity(pp::PointProcess{M}, history::History{M}, t) where {M}
    return ground_intensity(typeof(pp), get_θ(pp), history, t)
end

function ground_intensity_bound(pp::PointProcess{M}, history::History{M}, t) where {M}
    return ground_intensity_bound(typeof(pp), get_θ(pp), history, t)
end
