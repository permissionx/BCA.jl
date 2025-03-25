using QuadGK

function g(r::Float64, p_squared::Float64, E_r::Float64, type_p::Int64, type_t::Int64, constantsByType::ConstantsByType)
    value = 1 - p_squared / (r * r) - V(r, type_p, type_t, constantsByType) / E_r
    return value > 0.0 ? sqrt(value) : 0.0
end