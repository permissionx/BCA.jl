import Base: push!
using .BCA.ConstantFunctions
using LinearAlgebra


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
    radius, mass, Z = TypeToProperties(type)
    pValue = -1.0
    pVector = Vector{Float64}()
    pPoint = Vector{Float64}()  
    return Atom(id, type, coordinate, cellIndex, radius, mass, velocityDirection, energy, Z, pValue, pVector, pPoint)
end


function TypeToProperties(type::Int64)
    if type == 1
        radius = 1.0
        mass = 1.0
        Z = 1.0
    elseif type == 2
        radius = 2.0
        mass = 2.0
        Z = 2.0
    end
    return radius, mass, Z  
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

function CreateCellGrid(box::Box, inputVectors::Matrix{Float64})
    if !box.isOrthogonal
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

function InitConstants()
    p_max = 1.0
    q_max = 1.0
    stopEnergy = 1.0
    return Constants(p_max, q_max, q_max*q_max, stopEnergy)
end


function InitConstantsByType(types::Vector{Int64}, constants::Constants) 
    a_U = Dict{Vector{Int64}, Float64}()
    E_m = Dict{Int64, Float64}()
    S_e_UpTerm = Dict{Vector{Int64}, Float64}()
    S_e_DownTerm = Dict{Vector{Int64}, Float64}()
    x_nl = Dict{Vector{Int64}, Float64}()
    a = Dict{Vector{Int64}, Float64}()
    Q_nl = Dict{Vector{Int64}, Float64}()  
    Q_loc = Dict{Vector{Int64}, Float64}()
    for p in types
        for t in types
            _, mass_p, Z_p = TypeToProperties(p)
            _, _, Z_t = TypeToProperties(t)
            a_U[[p,t]] = BCA.ConstantFunctions.a_U(Z_p, Z_t)
            E_m[p] = BCA.ConstantFunctions.E_m(Z_p, mass_p)
            S_e_UpTerm[[p,t]] = BCA.ConstantFunctions.S_e_UpTerm(p, Z_p, Z_t, mass_p)
            x_nl[[p,t]] = BCA.ConstantFunctions.x_nl(p, Z_p, Z_t)
            a[[p,t]] = BCA.ConstantFunctions.a(Z_p, Z_t)
            Q_nl[[p,t]] = BCA.ConstantFunctions.Q_nl(Z_p, Z_t, constants.p_max)
            Q_loc[[p,t]] = BCA.ConstantFunctions.Q_loc(Z_p, Z_t)
        end
    end
    return ConstantsByType(a_U, E_m, S_e_UpTerm, S_e_DownTerm, x_nl, a, Q_nl, Q_loc)
end


function Simulator(box::Box, inputGridVectors::Matrix{Float64}, periodic::Vector{Bool}, types::Vector{Int64})
    cellGrid = CreateCellGrid(box, inputGridVectors)
    constants = InitConstants()
    constantsByType = InitConstantsByType(types, constants)
    return Simulator(Vector{Atom}(), box, cellGrid, periodic, box.isOrthogonal, 0, 0, 
                     types, constantsByType, constants)
end 


function Simulator(primaryVectors::Matrix{Float64}, boxSizes::Vector{Int64}, 
                   inputGridVectors::Matrix{Float64},
                   periodic::Vector{Bool}, 
                   latticeRanges::Matrix{Int64}, basis::Matrix{Float64}, basisTypes::Vector{Int64})
    box = CreateBoxByPrimaryVectors(primaryVectors, boxSizes)
    simulator = Simulator(box, inputGridVectors, periodic, unique(basisTypes))
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

function WhichCell(coordinate::Vector{Float64}, cellGrid::CellGrid)
    cellIndex = zeros(Int64, 3)
    for d in 1:3
        cellIndex[d] = Int64(floor(coordinate[d] / cellGrid.vectors[d,d])) + 1
        if cellIndex[d] < 1 
            cellIndex[d] = 1
        elseif cellIndex[d] > cellGrid.sizes[d]
            cellIndex[d] = cellGrid.sizes[d]
        end
    end
    return cellIndex
end

function push!(simulator::Simulator, atom::Atom)
    atom.index = simulator.maxAtomID + 1
    simulator.maxAtomID += 1
    push!(simulator.atoms, atom)
    simulator.numberOfAtoms += 1
    cellIndex = WhichCell(atom.coordinate, simulator.cellGrid)
    atom.cellIndex = cellIndex
    push!(simulator.cellGrid.cells[cellIndex[1], cellIndex[2], cellIndex[3]].atoms, atom)
end 



function delete!(simulator::Simulator, atom::Atom)
    filter!(a -> a.index != atom.index, simulator.atoms )
    filter!(a -> a.index != atom.index, simulator.cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]].atoms)
    simulator.numberOfAtoms -= 1
end

function DisplaceAtom(atom::Atom, newPosition::Vector{Float64}, simulator::Simulator)
    for d in 1:3
        if newPosition[d] < 0
            newPosition[d] += simulator.box.vectors[d,d]
        elseif newPosition[d] >= simulator.box.vectors[d,d]
            newPosition[d] -= simulator.box.vectors[d,d]
        end
    end
    atom.coordinate .= newPosition
    cellIndex = WhichCell(atom.coordinate, simulator.cellGrid)
    if cellIndex != atom.cellIndex
        ChangeAtomCell(atom, simulator.cellGrid, cellIndex)
    end
end


function ComputeDistance_squared(atom_p::Atom, atom_t::Atom, crossFlag::Vector{Int64}, box::Box)
    dv = VectorDifference(atom_p.coordinate, atom_t.coordinate, crossFlag, box)
    distance_squared = dv[1]* dv[1] + dv[2]*dv[2] + dv[3]  * dv[3]
    return distance_squared
end


function ComputeDistance(atom_p::Atom, atom_t::Atom, crossFlag::Vector{Int64}, box::Box)
    return sqrt(ComputeDistance_squared(atom_p, atom_t, crossFlag, box))
end


function ComputeVDistance(atom_p::Atom, atom_t::Atom, crossFlag::Vector{Int64}, box::Box)
    # v for atom_p
    dv = VectorDifference(atom_p.coordinate, atom_t.coordinate, crossFlag, box)
    return dv' * atom_p.velocityDirection
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


function ComputeP!(atom_p::Atom, atom_t::Atom, crossFlag::Vector{Int64}, box::Box)
    dv = VectorDifference(atom_p.coordinate, atom_t.coordinate, crossFlag, box)
    t = sum(dv .* atom_p.velocityDirection) / sum(atom_p.velocityDirection .* atom_p.velocityDirection)
    atom_t.pPoint = atom_p.coordinate + t * atom_p.velocityDirection
    atom_t.pVector = atom_t.pPoint - atom_t.coordinate
    atom_t.pValue = norm(atom_t.pVector)
end


function GetTargetsFromNeighbor(atom::Atom, gridCell::GridCell, simulator::Simulator)
    cellGrid = simulator.cellGrid
    box = simulator.box
    targets = Vector{Atom}()
    for (_, neighborCellInfo) in gridCell.neighborCellsInfo
        index = neighborCellInfo.index
        neighborCell = cellGrid.cells[index[1], index[2], index[3]]
        if neighborCell.isExplored
            continue
        end
        for neighborAtom in neighborCell.atoms
            if ComputeVDistance(atom, neighborAtom, neighborCellInfo.cross, box) > 0
                matchFlag = true
                for target in targets
                    if !SimultaneousCriteria(atom, neighborAtom, target, neighborCellInfo.cross, simulator)
                        matchFlag = false
                        break
                    end
                end
                if matchFlag
                    push!(targets, neighborAtom)
                end
            end
        end
        neighborCell.isExplored = true
    end
    return targets
end


function SimultaneousCriteria(atom::Atom, neighborAtom::Atom, addedTarget::Atom, crossFlag::Vector{Int64}, simulator::Simulator)
    box = simulator.box
    q_max_squared = simulator.constants.q_max_squared
    p_max = simulator.constants.p_max
    ComputeP!(atom, neighborAtom, crossFlag, box)
    flagP = abs(neighborAtom.pValue - addedTarget.pValue) < p_max
    flagQ = sum((neighborAtom.coordinate - addedTarget.coordinate) .* atom.velocityDirection) < q_max_squared
    return flagP && flagQ
end


function ChangeAtomCell(atom::Atom, cellGrid::CellGrid, nextCellIndex::Vector{Float64})
    filter!(a -> a.index != atom.index, cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]].atoms)
    atom.cellIndex = nextCellIndex
    nextCell = cellGrid.cells[cellIndex[1], cellIndex[2], cellIndex[3]]
    push!(nextCell.atoms, atom)
end



function SetvelocityDirection!(atom::Atom, velocity::Vector{Float64})
    atom.velocityDirection = velocity / norm(velocity)
end

function SetEnergy!(atom::Atom, energy::Float64)
    atom.energy = energy
end

function GetDTE(atom::Atom, simulator::Simulator)
    return 10.0
end