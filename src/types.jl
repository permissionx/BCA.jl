mutable struct Box
    vectors::Matrix{Float64}
    reciprocalVectors::Matrix{Float64}
    isOrthogonal::Bool
end


mutable struct Atom
    index::Int64  # never change
    isAlive::Bool
    type::Int64
    coordinate::Vector{Float64}
    cellIndex::Vector{Int64}
    radius::Float64
    mass::Float64
    velocityDirection::Vector{Float64}
    energy::Float64
    Z::Float64

    dte::Float64
    bde::Float64

    # for atom_t
    pValue::Dict{Int64, Float64}
    pPoint::Dict{Int64, Vector{Float64}}
    pVector::Dict{Int64, Vector{Float64}}
    pL::Dict{Int64, Float64}

    # for atom_p
    lastTargets::Vector{Int64}


    latticePointIndex::Int64 # -1 for off lattice
end

mutable struct LatticePoint
    index::Int64
    type::Int64  # Initial Type, will not change
    coordinate::Vector{Float64}
    cellIndex::Vector{Int64}

    atomIndex::Int64 # -1 for vacancy
end


struct NeighborCellInfo
    index::Vector{Int64}
    cross::Vector{Int64} # 0 for no cross, 1 for hi, -1 for lo, eg. [0,0,1] for top 
end



mutable struct GridCell
    # only for orthogonal box
    index::Vector{Int64}
    atoms::Vector{Int64}  
    latticePoints::Vector{Int64}  
    ranges::Matrix{Float64}
    centerCoordinate::Vector{Float64}
    neighborCellsInfo::Dict{Vector{Int64}, NeighborCellInfo}
    isExplored::Bool
    atomicDensity::Float64
end


mutable struct CellGrid
    cells::Array{GridCell, 3}
    vectors::Matrix{Float64}
    sizes::Vector{Int64}      
    cellVolume::Float64
end 


struct ConstantsByType
    V_upterm::Dict{Vector{Int64}, Float64}
    a_U::Dict{Vector{Int64}, Float64}
    E_m::Dict{Int64, Float64}
    S_e_upTerm::Dict{Vector{Int64}, Float64}
    S_e_downTerm::Dict{Vector{Int64}, Float64}
    x_nl::Dict{Vector{Int64}, Float64}
    a::Dict{Vector{Int64}, Float64}
    Q_nl::Dict{Vector{Int64}, Float64}  
    Q_loc::Dict{Vector{Int64}, Float64}
    qMax_squared::Dict{Vector{Int64}, Float64}
end


struct Constants
    pMax::Float64
    stopEnergy::Float64
    vacancyRecoverDistance_squared::Float64
    pLMax::Float64

    isDumpInCascade::Bool
    isLog::Bool
end


mutable struct Simulator
    atoms::Vector{Atom}
    latticePoints::Vector{LatticePoint}
    box::Box
    cellGrid::CellGrid
    periodic::Vector{Bool}
    isOrthogonal::Bool

    maxAtomID::Int64
    numberOfAtoms::Int64


    constantsByType::ConstantsByType
    constants::Constants

    isStore::Bool
    displacedAtoms::Vector{Int64}
    atomNumberWhenStore::Int64

    nIrradiation::Int64

    isDumpInCascade::Bool
    isLog::Bool
end


struct Parameters
    pMax::Float64
    stopEnergy::Float64
    vacancyRecoverDistance_squared::Float64
    pLMax::Float64

    isDumpInCascade::Bool
    isLog::Bool

    typeDict::Dict{Int, NamedTuple{(:radius, :mass, :Z, :dte, :bde, :alpha, :beta), Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64}}}
end

