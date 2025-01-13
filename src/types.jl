# need to combine the cell and box structs ?
mutable struct Box
    vectors::Matrix{Float64}
    reciprocalVectors::Matrix{Float64}
    isOrthogonal::Bool
end


mutable struct Atom
    index::Int64
    type::Int64
    coordinate::Vector{Float64}
    cellIndex::Vector{Int64}
    radius::Float64
    mass::Float64
    velocityDirection::Vector{Float64}
    energy::Float64
    Z::Float64

    pValue::Float64
    pPoint::Vector{Float64}
    pVector::Vector{Float64}
end


struct NeighborCellInfo
    index::Vector{Int64}
    cross::Vector{Int64} # 0 for no cross, 1 for hi, -1 for lo, eg. [0,0,1] for top 
end


mutable struct GridCell
    # only for orthogonal box
    index::Vector{Int64}
    atoms::Vector{Atom}
    ranges::Matrix{Float64}
    centerCoordinate::Vector{Float64}
    neighborCellsInfo::Dict{Vector{Int64}, NeighborCellInfo}
    isExplored::Bool
end


mutable struct CellGrid
    cells::Array{GridCell, 3}
    vectors::Matrix{Float64}
    sizes::Vector{Int64}      
end 

struct ConstantsByType
    a_U::Dict{Vector{Int64}, Float64}
    E_m::Dict{Int64, Float64}
    S_e_UpTerm::Dict{Vector{Int64}, Float64}
    S_e_DownTerm::Dict{Vector{Int64}, Float64}
    x_nl::Dict{Vector{Int64}, Float64}
    a::Dict{Int64, Float64}
    Q_nl::Dict{Int64, Float64}  
    Q_loc::Dict{Int64, Float64}
end

struct Constants
    p_max::Float64
    q_max::Float64
    q_max_squared::Float64
end

mutable struct Simulator
    atoms::Vector{Atom}
    box::Box
    cellGrid::CellGrid
    periodic::Vector{Bool}
    isOrthogonal::Bool

    maxAtomID::Int64
    numberOfAtoms::Int64

    types::Vector{Int64}

    constantsByType::ConstantsByType
    constants::Constants
end


