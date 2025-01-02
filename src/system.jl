function Box(Vectors::Matrix{Float64})
    # need to improve to detact if it is orithogonal. 
    return Box(Vectors, inv(Vectors'), true)
end 

function CreateBoxByPrimaryVectors(primaryVectors::Matrix{Float64}, sizes::Vector{Int64})
    vectors = primaryVectors .* sizes
    return Box(vectors)
end


function Atom(type::Int64, coordinate::Vector{Float64})
    id = 0
    cellIndex = Vector{Int64}()
    velocityDirection = Float64[0.0,0.0,0.0]  # length： 0 or one
    energy = 0.0
    radius, mass = TypeToRadiusMass(type)
    return Atom(id, type, coordinate, cellIndex, radius, mass, velocityDirection, energy)
end

function AssignVelocity(atom::Atom, energy::Float64, velocityDirection::Vector{Float64})
    atom.energy = energy
    atom.velocityDirection = velocityDirection
end

function IterPushCellNeighbors!(cellGrid::CellGrid, gridCell::GridCell, 
                                neighborKeys::Vector{Int64}, neighborIndex::Vector{Int64}, neighborCross::Vector{Int64}, 
                                nd::Int64)
    # Including self cell
    if nd <= 3
        for delta in [-1,0,1]
            index = gridCell.index[nd] + delta
            cross = 0
            if index < 1
                index += cellGrid.sizes[nd]
                cross = -1
            elseif index > cellGrid.sizes[nd]
                index -= cellGrid.sizes[nd]
                cross = 1
            end
            IterPushCellNeighbors!(cellGrid, gridCell, 
                                   push!(copy(neighborKeys), delta), 
                                   push!(copy(neighborIndex), index), 
                                   push!(copy(neighborCross), cross),
                                   nd+1)
        end
    else
        neighborCellInfo = NeighborCellInfo(neighborIndex, neighborCross)
        gridCell.neighborCellsInfo[neighborKeys] = neighborCellInfo
    end
end

function CreateCellGrid(box::Box, inputVectors::Matrix{Float64}, isOrthogonal::Bool)
    if !isOrthogonal        
        error("The box is not orthogonal, please use the orthogonal box.")
    end
    sizes = Vector{Int64}(undef, 3)
    vectors = Matrix{Float64}(undef, 3, 3)
    for d in 1:3
        sizes[d] = Int64(floor(box.vectors[d,d] / inputVectors[d,d]))
        vectors[d,d] = box.vectors[d,d] / sizes[d]
    end
    cells = Array{GridCell, 3}(undef, sizes[1], sizes[2], sizes[3])
    for x in 1:sizes[1]
        for y in 1:sizes[2]
            for z in 1:sizes[3]
                ranges = Matrix{Float64}(undef, 3, 2)
                ranges[1,1] = (x-1) * vectors[1,1]
                ranges[1,2] = x * vectors[1,1]
                ranges[2,1] = (y-1) * vectors[2,2]
                ranges[2,2] = y * vectors[2,2]
                ranges[3,1] = (z-1) * vectors[3,3]
                ranges[3,2] = z * vectors[3,3]  
                centerCoordinate = Vector{Float64}(undef, 3)  
                for d in 1:3
                    centerCoordinate[d] = (ranges[d,1] + ranges[d,2])/2
                end
                cells[x, y, z] = GridCell(Vector{Int64}([x,y,z]),Vector{Atom}(), ranges, centerCoordinate, Dict{Vector{Int64}, Vector{Int64}}(), false)
            end
        end    
    end
    cellGrid = CellGrid(cells, vectors, sizes) 
    for cell in cellGrid.cells
        IterPushCellNeighbors!(cellGrid, cell, Vector{Int64}(), Vector{Int64}(), Vector{Int64}(), 1)
    end
    return cellGrid
end

function Simulator(box::Box, inputGridVectors::Matrix{Float64}, periodic::Vector{Bool})
    cellGrid = CreateCellGrid(box, inputGridVectors, box.isOrthogonal)
    return Simulator(Vector{Atom}(), box, cellGrid, periodic, box.isOrthogonal, 0, 0)
end 

import Base: push!

function push!(simulator::Simulator, atom::Atom)
    atom.index = simulator.maxAtomID + 1
    simulator.maxAtomID += 1
    push!(simulator.atoms, atom)
    simulator.numberOfAtoms += 1
    
    cellIndex = zeros(Int64, 3)
    for d in 1:3
        cellIndex[d] = Int64(floor(atom.coordinate[d] / simulator.cellGrid.vectors[d,d])) + 1
        if cellIndex[d] < 1 
            cellIndex[d] = 1
        elseif cellIndex[d] > simulator.cellGrid.sizes[d]
            cellIndex[d] = simulator.cellGrid.sizes[d]
        end
    end
    
    atom.cellIndex = cellIndex
    push!(simulator.cellGrid.cells[cellIndex[1], cellIndex[2], cellIndex[3]].atoms, atom)
end 



function delete!(simulator::Simulator, atom::Atom)
    filter!(a -> a.index != atom.index, simulator.atoms )
    filter!(a -> a.index != atom.index, simulator.cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]].atoms)
end


function Simulator(primaryVectors::Matrix{Float64}, boxSizes::Vector{Int64}, 
                   inputGridVectors::Matrix{Float64},
                   periodic::Vector{Bool}, 
                   latticeRanges::Matrix{Int64}, basis::Matrix{Float64}, basisTypes::Vector{Int64})
    box = CreateBoxByPrimaryVectors(primaryVectors, boxSizes)
    simulator = Simulator(box, inputGridVectors, periodic)
    for x in latticeRanges[1,1]:latticeRanges[1,2]-1
        for y in latticeRanges[2,1]:latticeRanges[2,2]-1    
            for z in latticeRanges[3,1]:latticeRanges[3,2]-1
                for i in 1:length(basisTypes)
                    reducedCoordinate = Float64[x,y,z] + basis[i, :]
                    coordinate = primaryVectors' * reducedCoordinate
                    atom = Atom(basisTypes[i], coordinate)
                    push!(simulator, atom)
                end
            end
        end
    end
    return simulator
end 


function ComputeDistance_squared(atom1::Atom, atom2::Atom, crossFlag::Vector{Int64}, box::Box)
    dv = VectorDifference(atom1.coordinate, atom2.coordinate, crossFlag, box)
    distance_squared = dv[1]* dv[1] + dv[2]*dv[2] + dv[3]  * dv[3]
    return distance_squared
end


function ComputeDistance(atom1::Atom, atom2::Atom, crossFlag::Vector{Int64}, box::Box)
    return sqrt(ComputeDistance_squared(atom1, atom2, crossFlag, box))
end


function ComputeVDistance(atom1::Atom, atom2::Atom, crossFlag::Vector{Int64}, box::Box)
    # v for atom1
    dv = VectorDifference(atom1.coordinate, atom2.coordinate, crossFlag, box)
    return dv' * atom1.velocityDirection
end


function VectorDifference(v1::Vector{Float64}, v2::Vector{Float64}, crossFlag::Vector{Int64}, box::Box)
    # return v2 - v1
    if crossFlag == Vector{Int64}([0,0,0])
        return v2 - v1
    end 
    result = Vector{Float64}(undef, 3)
    for d in 1:3
        dv = v2[d] - v1[d] + crossFlag[d] * box.vectors[d,d]
        result[d] = dv
    end
    return result
end


function ComputeDistanceToRay_squared(atom1::Atom, atom2::Atom, crossFlag::Vector{Int64}, box::Box)
    # atom2 to atom1's ray
    dv = VectorDifference(atom1.coordinate, atom2.coordinate, crossFlag, box)
    area = cross(dv, atom1.velocityDirection)
    distance_squared = (area[1]*area[1] + area[2]*area[2] + area[3]*area[3]) / (atom1.velocityDirection[1]*atom1.velocityDirection[1] + atom1.velocityDirection[2]*atom1.velocityDirection[2] + atom1.velocityDirection[3]*atom1.velocityDirection[3])
    return distance_squared
end


function GetTargetFromNeighbor(atom::Atom, gridCell::GridCell, cellGrid::CellGrid, box::Box)
    minDistance = Inf 
    target = atom
    for (_, neighborCellInfo) in gridCell.neighborCellsInfo
        index = neighborCellInfo.index
        neighborCell = cellGrid.cells[index[1], index[2], index[3]]
        if !neighborCell.isExplored
            cross = neighborCellInfo.cross
            neighborCell.isExplored = true
            for neighborAtom in neighborCell.atoms
                if ComputeVDistance(atom, neighborAtom, cross, box) > 0
                    distance = ComputeDistanceToRay_squared(atom, neighborAtom, cross, box)
                    if distance < minDistance
                        target = neighborAtom
                        minDistance = distance
                    end
                end
            end
        end
    end
    return target
end

function ChangeAtomCell(atom::Atom, cellGrid::CellGrid, nextCellIndex::Vector{Float64})
    filter!(a -> a.index != atom.index, cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]].atoms)
    atom.cellIndex = nextCellIndex
    nextCell = cellGrid.cells[cellIndex[1], cellIndex[2], cellIndex[3]]
    push!(nextCell.atoms, atom)
end


function AtomOutFaceDimension(atom::Atom, cell::GridCell)
    for d in 1:3
        if atom.velocityDirection[d] >= 0
            rangeIndex = 2
        else
            rangeIndex = 1
        end
        faceCoordinate = cell.ranges[d, rangeIndex]
        t = (faceCoordinate - atom.coordinate[d]) / atom.velocityDirection[d]
        elseDs = [ed for ed in 1:3 if ed != d]
        allInRange = true
        for elseD in elseDs
            crossCoord = atom.coordinate[elseD] + atom.velocityDirection[elseD] * t
            if !(cell.ranges[elseD, 1] <= crossCoord <= cell.ranges[elseD, 2])
                allInRange = false
                break
            end
        end
        if allInRange
            return d, rangeIndex
        end
    end
    error("Atom $(atom.index) out face not found")
end


function ShotTarget(atom::Atom, simulator::Simulator)
    cellGrid = simulator.cellGrid
    box = simulator.box
    periodic = simulator.periodic       
    cell = cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    while true
        target = GetTargetFromNeighbor(atom, cell, cellGrid, box)
        if target.index != atom.index
            for cell in cellGrid.cells
                cell.isExplored = false
            end 
            return target
        else
            dimension, direction = AtomOutFaceDimension(atom, cell)
            neighborIndex = Vector{Int64}([0,0,0])
            neighborIndex[dimension] = direction == 1 ? -1 : 1
            neighborInfo = cell.neighborCellsInfo[neighborIndex]
            if !periodic[dimension]
                if neighborInfo.cross[dimension] != 0
                    for cell in cellGrid.cells
                        cell.isExplored = false
                    end 
                    return atom # means find nothing 
                end
            end 
            index = neighborInfo.index
            cell = cellGrid.cells[index[1], index[2], index[3]]
            println(cell.index)
        end
    end
end


function TypeToRadiusMass(type::Int64)
    if type == 1
        radius = 1.0
        mass = 1.0
    elseif type == 2
        radius = 2.0
        mass = 2.0
    end
    return radius, mass
end 

