module BCA
using LinearAlgebra
using QuadGK
using Base.MathConstants
using Main: ConstantsByType

export CollisionParams

qe_squared = Float64(14.399764) # square of element charge, unit: eV*angstrom

function Φ(x::Float64)
    A = Vector{Float64}([0.1818, 0.5099, 0.2802, 0.02817])
    B = Vector{Float64}([3.2, 0.9423, 0.4028, 0.2016])
    result = 0.0
    for i in 1:4
        result += A[i] * exp(-B[i] * x)
    end
    return result
end

function V(r::Float64, type_p::Int64, type_t::Int64, constantsByType::ConstantsByType)
    return constantsByType.V_upterm[[type_p, type_t]] / r * Φ(r / constantsByType.a_U[[type_p, type_t]])  
end

function E_r(energy_p::Float64, mass_p::Float64, mass_t::Float64)
    return mass_t * energy_p / (mass_p + mass_t)
end


function g(r::Float64, p_squared::Float64, E_r::Float64, type_p::Int64, type_t::Int64, constantsByType::ConstantsByType)
    value = 1 - p_squared / (r * r) - V(r, type_p, type_t, constantsByType) / E_r
    return value > 0.0 ? sqrt(value) : 0.0
end


function FindTurningPoint(p_squared::Float64, E_r::Float64,type_p::Int64, type_t::Int64,  guessStart::Float64, constantsByType::ConstantsByType)
    rLeft = guessStart
    rRight = 10.0
    while (rRight - rLeft) > 1e-15
        rMid = (rLeft + rRight) / 2
        g_squared = 1 - p_squared / (rMid * rMid) - V(rMid, type_p, type_t, constantsByType) / E_r
        if g_squared < 0
            rLeft = rMid
        elseif g_squared > 0
            rRight = rMid
        else
            return rMid+1e-15
        end
    end
    return rRight
end


function Integrate_g_θ(p_squared::Float64, E_r::Float64, type_p::Int64, type_t::Int64, rStart::Float64, constantsByType::ConstantsByType)
    result = quadgk(r -> 1 / (r * r * g(r, p_squared, E_r, type_p, type_t, constantsByType)),
        rStart, rStart+10000,
        rtol=1e-8)[1]
    return result
end


function θ(p::Float64, p_squared::Float64, E_r::Float64, type_p::Int64, type_t::Int64, rStart::Float64, constantsByType::ConstantsByType)
    return π - 2 * p * Integrate_g_θ(p_squared, E_r, type_p, type_t, rStart, constantsByType)
end


function Integrate_g_τ(p_squared::Float64, E_r::Float64, type_p::Int64, type_t::Int64, rStart::Float64, constantsByType::ConstantsByType)
    result = quadgk(r -> 1 / g(r, p_squared, E_r, type_p, type_t, constantsByType) - 1 / sqrt(1 - p_squared / (r * r)),
        rStart, rStart+1000,
        rtol=1e-8)[1]
    return result
end


function τ(p_squared::Float64, type_p::Int64, type_t::Int64, E_r::Float64, rStart::Float64, constantsByType::ConstantsByType)
    return sqrt(rStart * rStart - p_squared) - Integrate_g_τ(p_squared, E_r, type_p, type_t, rStart, constantsByType)
end


function tanφ(mass_p::Float64, mass_t::Float64, θ::Float64)
    return mass_t * sin(θ) / (mass_p + mass_t * cos(θ))
end

function tanψ(θ::Float64)
    return sin(θ) / (1 - cos(θ))
end

function E_t(mass_p::Float64, mass_t::Float64, energy_p::Float64, θ::Float64)
    return 4 * mass_p * mass_t / (mass_p + mass_t)^2 * energy_p * sin(θ / 2)^2
end

function E_p(energy_p::Float64, E_t::Float64, Q::Float64)
    return energy_p - E_t - Q
end

function x_p(mass_p::Float64,mass_t::Float64, p::Float64, θ::Float64, τ::Float64)
    return (2 * mass_p * τ + (mass_t - mass_p) * p * tan(θ / 2)) / (mass_p + mass_t)
end

function x_t(p::Float64, θ::Float64, x_p::Float64)
    return p * tan(θ / 2) - x_p
end

function θτ(energy_p::Float64, mass_p::Float64, mass_t::Float64, type_p::Int64, type_t::Int64,
    p::Float64, constantsByType::ConstantsByType)
    p_squared = p * p
    E_r_v = E_r(energy_p, mass_p, mass_t)
    rStart = FindTurningPoint(p_squared, E_r_v, type_p, type_t, p, constantsByType)
    θ_v = θ(p, p_squared, E_r_v, type_p, type_t, rStart, constantsByType)
    τ_v = τ(p_squared, type_p, type_t, E_r_v, rStart, constantsByType)
    return θ_v, τ_v
end


function CollisionParams(energy_p::Float64, mass_p::Float64, mass_t::Float64, type_p::Int64, type_t::Int64,
                         p::Float64, pL::Float64, N::Float64, constantsByType::ConstantsByType,
                         θFunction::Function, τFunction::Function)
    #=
    energy_p = 100.0
    mass_p = 4.0
    mass_t = 12.0
    type_p = 1
    type_t = 2
    p = 2.0
    pL = 1.0
    N = 1.0
    =#
    #p_squared = p * p
    E_r_v = E_r(energy_p, mass_p, mass_t)
    #rStart = FindTurningPoint(p_squared, E_r_v, type_p, type_t, p, constantsByType)
    #θ_v = θ(p, p_squared, E_r_v, type_p, type_t, rStart, constantsByType)
    #τ_v = τ(p_squared, type_p, type_t, E_r_v, rStart, constantsByType)
    θ_v = θFunction(energy_p, p)
    τ_v = τFunction(energy_p, p)
    tanφ_v = tanφ(mass_p, mass_t, θ_v)
    tanψ_v = tanψ(θ_v)
    E_t_v = E_t(mass_p, mass_t, energy_p, θ_v)
    x_p_v = x_p(mass_p, mass_t, p, θ_v, τ_v)
    x_t_v = x_t(p, θ_v, x_p_v)
    Q_v = QLoss.Q(energy_p, type_p, type_t, E_r_v, p, pL, N, constantsByType)
    #println("E_r: ", E_r_v,"\n", "rStart: ", rStart,"\n", "θ: ", θ_v,"\n","τ: ", τ_v,"\n", "tanφ: ", tanφ_v,"\n", "tanψ: ", tanψ_v,"\n", "E_t: ", E_t_v,"\n", "x_p: ", x_p_v,"\n", "x_t: ", x_t_v,"\n", "Q: ", Q_v)
    #println("\n")
    return tanφ_v, tanψ_v,  E_t_v, x_p_v, x_t_v, Q_v
end




module QLoss
using Main: ConstantsByType
# 
# constants for different types:
function S_e(energy_p::Float64, type_p::Int64, type_t::Int64, constantsByType::ConstantsByType)
    δ_1_2 = Float64(0.7125) # δ^(1/2)
    δ_1_2_ = Float64(-0.7125) # -δ^(1/2)
    δ_1 = Float64(0.7017543859649122) # 1/δ
    A = energy_p / constantsByType.E_m[type_p]
    termLeftDown = (A / log(A + 1 / A + ℯ - 2))^δ_1_2 
    termRightDown = A ^ δ_1_2_
    termUp = constantsByType.S_e_upTerm[[type_p, type_t]]
    return termUp / (termLeftDown + termRightDown)^δ_1
end


function x_nl(type_p::Int64, type_t::Int64, E_r::Float64, constantsByType::ConstantsByType)
    # termOthers 
    termE_r = E_r^0.075
    termOthers = constantsByType.x_nl[[type_p, type_t]]
    return termE_r * termOthers
end 

function x_loc(x_nl::Float64)
    return 1 - x_nl
end


function Q_nl(type_p::Int64, type_t::Int64, S_e::Float64, x_nl::Float64, x_loc::Float64, pL::Float64, N::Float64, constantsByType::ConstantsByType)
    # constant:  pMax (half of lattice constant)
    termRight = x_nl + x_loc * constantsByType.Q_nl[[type_p, type_t]]
    return S_e * N * termRight * pL 
end 

function Q_loc(type_p::Int64, type_t::Int64, S_e::Float64, x_loc::Float64, p::Float64, constantsByType::ConstantsByType)
    termUp = x_loc * S_e * exp(-p / constantsByType.a[[type_p, type_t]])
    termDown = constantsByType.Q_loc[[type_p, type_t]]
    return termUp * termDown
end


function Q(energy_p::Float64, type_p::Int64, type_t::Int64, E_r::Float64, p::Float64, pL::Float64, N::Float64, constantsByType::ConstantsByType)
    S_e_v = S_e(energy_p, type_p, type_t, constantsByType)
    x_nl_v = x_nl(type_p, type_t, E_r, constantsByType)
    x_loc_v = x_loc(x_nl_v)
    Q_nl_v = Q_nl(type_p, type_t, S_e_v, x_nl_v, x_loc_v, pL, N, constantsByType)
    Q_loc_v = Q_loc(type_p, type_t, S_e_v, x_loc_v, p, constantsByType)
    # println("S_e: ", S_e_v, "\n", "x_nl: ", x_nl_v, "\n", "x_loc: ", x_loc_v, "\n", "Q_nl: ", Q_nl_v, "\n", "Q_loc: ", Q_loc_v)
    return Q_nl_v + Q_loc_v 
end 

end

module ConstantFunctions
using Base.MathConstants
using ..BCA: qe_squared


function V_upterm(Z_p::Float64, Z_t::Float64)
    return Z_p * Z_t * qe_squared
end

function a_U(Z_p::Float64, Z_t::Float64)
    a_0 = 0.529177210903::Float64  # Bohr radius, unit: angstrom
    return 0.8854 * a_0 / (Z_p^0.23 + Z_t^0.23)
end

function E_m(Z_p::Float64, mass_p::Float64)
    v_b_squared = Float64(49708.0)  # square of Bohr velocity, unit: eV/u need to check carefully
    return 0.5 * mass_p * v_b_squared * Z_p^(4/3)
end

function S_e_upTerm(type_p::Int64,Z_p::Float64, Z_t::Float64, m_p::Float64, α_p::Float64)
    # In initializations
    k_LS = 1.212 * Z_p^(7/6) * Z_t / ((Z_p^(2/3) + Z_t^(2/3))^(3/2) * m_p^(1/2))
    return α_p * k_LS * sqrt(E_m(Z_p, m_p))
end

function x_nl(type_p::Int64, Z_p::Float64, Z_t::Float64, β_p::Float64)
    return β_p * (a_U(Z_p, Z_t) / (Z_p * Z_t * qe_squared))^0.075
end

function a(Z_p::Float64, Z_t::Float64)
    return 1.45 * a_U(Z_p, Z_t) / 0.3 / Z_p^0.4
end

function Q_nl(Z_p::Float64, Z_t::Float64, pMax::Float64)
    return (1 + pMax / a(Z_p, Z_t)) * exp(-pMax / a(Z_p, Z_t))
end

function Q_loc(Z_p::Float64, Z_t::Float64)
    return 1 / (2 * π * a(Z_p, Z_t)^2)
end

end


end
