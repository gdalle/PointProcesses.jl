abstract type PointProcess{MarkType} end

function ground_intensity(pointprocess::PointProcess, history, t)
    error("not implemented")
end

function ground_intensity_bound(pointprocess::PointProcess, history, t)
    error("not implemented")
end

function ground_intensity_bound_validity_duration(pointprocess::PointProcess, history, t)
    error("not implemented")
end

function mark_distribution(pointprocess::PointProcess, history, t)
    error("not implemented")
end

function mark_density(pointprocess::PointProcess, history, t, m)
    return pdf(mark_distribution(pointprocess, history, t), m)
end

function intensity(pointprocess::PointProcess, history, t, m)
    return ground_intensity(pointprocess, history, t) * mark_density(pointprocess, history, t, m)
end
