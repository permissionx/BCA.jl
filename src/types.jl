# need to combine the cell and box structs ?
mutable struct Box
    vectors::Matrix{Float64}
    reciprocalVectors::Matrix{Float64}
    isOrthogonal::Bool
end


mutable struct Atom
    id::Int64
    type::Int64
    coordinate::Vector{Float64}
    cellIndex::Vector{Int64}
    radius::Float64
    mass::Float64
    velocityDirection::Vector{Float64}
    energy::Float64
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

mutable struct Simulator
    atoms::Vector{Atom}
    box::Box
    cellGrid::CellGrid
    periodic::Vector{Bool}
    isOrthogonal::Bool

    maxAtomID::Int64
    numberOfAtoms::Int64
end
