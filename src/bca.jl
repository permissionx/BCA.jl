module BCA
using LinearAlgebra
using QuadGK
using Base.MathConstants
using Main: Atom, Simulator
using .QLoss

export Collision

qe_squared = Float64(1.4399764) # square of element charge, unit: eV*angstrom

function Φ(x::Float64)
    A = Vector{Float64}([0.1818, 0.5099, 0.2802, 0.02817])
    B = Vector{Float64}([3.2, 0.9423, 0.4028, 0.2016])
    result = 0.0
    for i in 1:4
        result += A[i] * exp(-B[i] * x)
    end
    return result
end

function V(r::Float64, Z_p::Float64, Z_t::Float64, type_P::Int64, type_T::Int64, simulator::Simulator)
    return Z_p * Z_t * qe_squared / r * Φ(r / simulator.constantsByType.a_U[[type_P, type_T]])
end

function E_r(atom_p::Atom, atom_t::Atom)
    return atom_t * atom_p.energy / (atom_p.mass + atom_t.mass)
end

function g(r::Float64, p_squared::Float64, E_r::Float64, Z_p::Float64, Z_t::Float64, type_P::Int64, type_T::Int64, simulator::Simulator)
    value = 1 - p_squared / (r * r) - V(r, Z_p, Z_t, type_P, type_T, simulator) / E_r
    return value > 0.0 ? sqrt(value) : 0.0
end

function FindTurningPoint(p_squared::Float64, E_r::Float64, Z_p::Float64, Z_t::Float64,type_P::Int64, type_T::Int64,  guessStart::Float64, simulator::Simulator)
    rLeft = guessStart
    rRight = 10.0
    while (rRight - rLeft) > 1e-10
        rMid = (rLeft + rRight) / 2
        if (1 - p_squared / (rMid * rMid) - V(rMid, Z_p, Z_t, type_P, type_T, simulator) / E_r) < 0
            rLeft = rMid
            elses
            rRight = rMid
        end
    end
    return rRight + 1e-6
end

function Integrate_g_θ(p_squared::Float64, E_r::Float64, Z_p::Float64, Z_t::Float64, type_P::Int64, type_T::Int64, rStart::Float64, simulator::Simulator)::Float64
    result = quadgk(r -> 1 / (r * r * g(r, p_squared, E_r, Z_p, Z_t, type_P, type_T, simulator)),
        rStart, 1000.0,
        rtol=1e-8)[1]
    return result
end


function θ(p::Float64, p_squared::Float64, atom_p::Atom, atom_t::Atom, rStart::Float64, simulator::Simulator)
    E_r = E_r(atom_p, atom_t)
    Z_p = atom_p.Z
    Z_t = atom_t.Z
    type_P = atom_p.type
    type_T = atom_t.type
    return π - 2 * p * Integrate_g_θ(p_squared, E_r, Z_p, Z_t, type_P, type_T, rStart, simulator)
end


function Integrate_g_τ(p_squared::FLoat64, E_r::Float64, Z_p::Float64, Z_t::Float64, type_P::Int64, type_T::Int64, rStart::Float64, simulator::Simulator)::Float64
    result = quadgk(r -> 1 / g(r, p_squared, E_r, Z_p, Z_t, type_P, type_T, simulator) - 1 / sqrt(1 - p_squared / (r * r)),
        rStart, 1000.0,
        rtol=1e-8)[1]
    return result
end


function τ(p_squared::Float64, atom_p::Atom, atom_t::Atom, rStart::Float64, simulator::Simulator)::Float64
    E_r = E_r(atom_p, atom_t)
    Z_p = atom_p.Z
    Z_t = atom_t.Z
    type_P = atom_p.type
    type_T = atom_t.type
    return sqrt(rStart * rStart - p_squared) - Integrate_g_τ(p_squared, E_r, Z_p, Z_t, type_P, type_T, rStart, simulator)
end


function tanφ(atom_p::Atom, atom_t::Atom, θ::Float64)
    return atom_t.mass * sin(θ) / (atom_p.mass + atom_t.mass * cos(θ))
end

function tanψ(θ::Float64)
    return sin(θ) / (1 - cos(θ))
end

function E_t(atom_p::Atom, atom_t::Atom, θ::Float64)
    return 4 * atom_p.mass * atom_t.mass / (atom_p.mass + atom_t.mass)^2 * atom_p.energy * sin(θ / 2)^2
end

function E_p(atom_p::Atom, E_t::Float64, Q::Float64)
    return atom_p.energy - E_t - Q
end

function x_p(atom_p::Atom, atom_t::Atom, p::Float64, θ::Float64, τ::Float64)
    return (2 * atom_p.mass * τ + (atom_t.mass - atom_p.mass) * p * tan(θ / 2)) / (atom_p.mass + atom_t.mass)
end

function x_t(p::Float64, θ::Float64, x_p::Float64)
    return p * tan(θ / 2) - x_p
end

function CollisionParams(atom_p::Atom, atom_t::Atom, p::Float64, p_squared::Float64, simulator::Simulator)
    E_r = E_r(atom_p, atom_t)
    rStart = FindTurningPoint(p_squared, E_r, atom_p.Z, atom_t.Z, atom_p.type, atom_t.type, p, simulator)
    θ = θ(p, p_squared, atom_p, atom_t, rStart, simulator)
    τ = τ(p_squared, atom_p, atom_t, rStart, simulator)
    tanφ = tanφ(atom_p, atom_t, θ)
    tanψ = tanψ(θ)
    E_t = E_t(atom_p, atom_t, θ)
    Q = QLoss.Q(atom_p, atom_t, E_r, p, simulator)
    E_p = E_p(atom_p, E_t, Q)
    x_p = x_p(atom_p, atom_t, p, θ, τ)
    x_t = x_t(p, θ, x_p)
    return tanφ, tanψ, E_t, E_p, x_p, x_t, Q
end





module QLoss
# atomic_density = 0.25::Float64 # atomic density, unit: 1/angstrom^3
# 
# constants for different types:
function S_e(atom_p::Atom, atom_t::Atom, simulator::Simulator)
    δ_1_2 = Float64(1.1937336386313322) # δ^(1/2)
    δ_1_2_ = Float64(-1.1937336386313322) # -δ^(1/2)
    δ_1 = Float64(0.7017543859649122) # 1/δ
    A = atom_p.energy / simulator.constantsByType.E_m[atom_p.type]
    termLeftDown = (A / log(A + 1 / A + ℯ - 2))^δ_1_2 
    termRightDown = A ^ δ_1_2_
    termUp = simulator.constantsByType.S_e_upTerm[[atom_p.type, atom_t.type]]
    return termUp / (termLeftDown + termRightDown)^δ_1
end


function x_nl(atom_p::Atom, atom_t::Atom,  E_r::Float64, simulator::Simulator)
    # termOthers 
    termE_r = E_r^0.075
    termOthers = simulator.constantsByType.x_nl[[atom_p.type, atom_t.type]]
    return termE_r * termOthers
end 

function x_loc(x_nl::Float64)
    return 1 - x_nl
end


function Q_nl(atom_p::Atom, S_e::Float64, x_nl::Float64, x_loc::Float64, simulator::Simulator)
    # constant: N_density, p_max (half of lattice constant)
    termRight = x_nl + x_loc * simulator.constantsByType.Q_nl[atom_p.type]
    return S_e * simulator.constantsByType.N_density * termRight
end 

function Q_loc(atom_p::Atom, S_e::Float64, x_loc::Float64, p::Float64, simulator::Simulator)
    termUp = x_loc * S_e * exp(-p / simulator.constantsByType.a[atom_p.type])
    termDown = simulator.constantsByType.Q_loc[satom_p.type]
    return termUp * termDown
end


function Q(atom_p::Atom, atom_t::Atom, E_r::Float64, p::Float64, simulator::Simulator)
    S_e = S_e(atom_p, atom_t, simulator)
    x_nl = x_nl(atom_p, atom_t, E_r, simulator)
    x_loc = x_loc(x_nl)
    Q_nl = Q_nl(atom_p, S_e, x_nl, x_loc, simulator)
    Q_loc = Q_loc(atom_p, S_e, x_loc, p, simulator)
    return Q_nl + Q_loc 
end 

end

module Constants
using Base.MathConstants

function α(type_p::Int64)
    return type_p / 137.035999074
end

function β(type_p::Int64)
    return type_p / 137.035999074
end

function a_U(Z_p::Float64, Z_t::Float64)
    a_0 = 0.529177210903::Float64 # Bohr radius, unit: angstrom
    return 0.8854 * a_0 / (Z_p^0.23 + Z_t^0.23)
end

function E_m(Z_p::Float64, mass_p::Float64)
    v_b_squared = Float64(0.113796)  # square of Bohr velocity, unit: eV need to check carefully
    return 0.5 * mass_p * v_b_squared * Z_p^(4/3)
end

function S_e_UpTerm(type_p::Int64,Z_p::Float64, Z_t::Float64, m_p::Float64)
    # In initializations
    k_LS = 1.212 * Z_p^(7/6) * Z_t / ((Z_p^(3/2) + Z_t^(3/2))^(3/2) * m_p^(1/2))
    return α(type_p) * k_LS * E_m(Z_p, m_p)
end

function x_nl(type_p::Int64, Z_p::Float64, Z_t::Float64)
    using ..BCA: qe_squared
    return β(type_p) * (a_U(Z_p, Z_t) / (Z_p * Z_t * qe_squared))^0.075
end

function a(Z_p::Float64)
    return 1.45 * a_U(Z_p, Z_t) / 0.3 / Z_p^0.4
end

function Q_nl(Z_p::Float64, p_max::Float64)
    a = a(Z_p)
    return (1 + p_max / a) * exp(-p_max / a)
end

function Q_loc(Z_p::Float64)
    a = a(Z_p)
    return 1 / (2 * π * a * a)
end

end
end
