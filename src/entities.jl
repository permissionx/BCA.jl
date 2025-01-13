import Base: push!
using .BCA
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
    velocityNorm = Float64[0.0,0.0,0.0]  # length： 0 or one
    energy = 0.0
    radius, mass, Z = TypeToProperties(type)
    pValue = -1.0
    pVector = Vector{Float64}()
    pPoint = Vector{Float64}()  
    return Atom(id, type, coordinate, cellIndex, radius, mass, velocityNorm, energy, Z, pValue, pVector, pPoint)
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
    return Constants(p_max, q_max, q_max*q_max)
end


function InitConstantsByType(types::Vector{Int64}, constants::Constants) 
    using .BCA.Constants: a_U, E_m, S_e_UpTerm, x_nl, a, Q_nl, Q_loc
    for p in types
        for t in types
            _, mass_p, Z_p = TypeToProperties(p)
            _, _, Z_t = TypeToProperties(t)
            a_U[[p,t]] = a_U(Z_p, Z_t)
            E_m[p] = E_m(Z_p, mass_p)
            S_e_UpTerm[[p,t]] = S_e_UpTerm(p, Z_p, Z_t, mass_p)
            x_nl[[p,t]] = x_nl(p, Z_p, Z_t)
            a[p] = a(Z_p)
            Q_nl[p] = Q_nl(Z_p, constants.p_max)
            Q_loc[p] = Q_loc(Z_p, constants.p_max)
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
    return dv' * atom_p.velocityNorm
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
    t = sum(dv .* atom_p.velocityNorm) / sum(atom_p.velocityNorm .* atom_p.velocityNorm)
    atom.pPoint = atom_p.coordinate + t * atom_p.velocityNorm
    atom.pVector = atom.pPoint - atom_t.coordinate
    atom.pValue = norm(atom.pVector)
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
                    neighborAtom.pVector, neighborAtom.pPoint = ComputePVector(atom, neighborAtom, neighborCellInfo.cross, box)
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
    flagQ = sum((neighborAtom.coordinate - addedTarget.coordinate) .* atom.velocityNorm) < q_max_squared
    return flagP && flagQ
end


function ChangeAtomCell(atom::Atom, cellGrid::CellGrid, nextCellIndex::Vector{Float64})
    filter!(a -> a.index != atom.index, cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]].atoms)
    atom.cellIndex = nextCellIndex
    nextCell = cellGrid.cells[cellIndex[1], cellIndex[2], cellIndex[3]]
    push!(nextCell.atoms, atom)
end


function AtomOutFaceDimension(atom::Atom, cell::GridCell)
    for d in 1:3
        if atom.velocityNorm[d] >= 0
            rangeIndex = 2
        else
            rangeIndex = 1
        end
        faceCoordinate = cell.ranges[d, rangeIndex]
        t = (faceCoordinate - atom.coordinate[d]) / atom.velocityNorm[d]
        elseDs = [ed for ed in 1:3 if ed != d]
        allInRange = true
        for elseD in elseDs
            crossCoord = atom.coordinate[elseD] + atom.velocityNorm[elseD] * t
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
    periodic = simulator.periodic       
    cell = cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    while true
        targets = GetTargetsFromNeighbor(atom, cellGrid, simulator)
        if length(targets) > 0
            for cell in cellGrid.cells
                cell.isExplored = false
            end
            #for atom in simulator.atoms
            #    atom.pValue = -1.0
            #end
            return targets
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
                    #for atom in simulator.atoms
                    #    atom.pValue = -1.0
                    #end
                    return Vector{Atom}() # means find nothing 
                end
            end 
            index = neighborInfo.index
            cell = cellGrid.cells[index[1], index[2], index[3]]
        end
    end
end


function Collision!(atom_p::Atom, atoms_t::Vector{Atom}, simulator::Simulator)
    N_t = length(atoms_t)
    tanφList = Vector{Float64}(undef, N_t)
    tanψList = Vector{Float64}(undef, N_t)
    E_tList = Vector{Float64}(undef, N_t)
    E_pList = Vector{Float64}(undef, N_t)
    x_pList = Vector{Float64}(undef, N_t)
    x_tList = Vector{Float64}(undef, N_t)
    QList = Vector{Float64}(undef, N_t)
    for (i, atom_t) in enumerate(atoms_t)
        tanφList[i], tanψList[i], E_tList[i], E_pList[i], x_pList[i], x_tList[i], QList[i] = CollisionParams(atom_p, atom_t, atom_t.pValue, atom_t.pValue * atom_t.pValue, simulator)
    end
    sumE_t = sum(E_tList)
    η = N_t * atom_p.energy / (N_t * atom_p.energy + (N_t - 1) * sumE_t)

    # Update atoms_t (target atoms)
    avePPoint = Vector{Float64}([0,0,0])
    momentum = Vector{Float64}([0,0,0])
    for (i, atom_t) in enumerate(atoms_t)
        tCoordinate = atom_t.coordinate + x_tList[i] * η * atom_p.velocityNorm
        DisplaceAtom(atom_t, tCoordinate, simulator)
        velocityNorm = -atom_t.pVector / norm(atom_t.pVector) * tanψList[i] + atom_t.velocityNorm
        atom_t.velocityNorm = velocityNorm / norm(velocityNorm)
        atom_t.energy = E_tList[i] * η
        avePPoint += atom_t.pPoint
        momentum += sqrt(2 * atom_t.mass * atom_t.energy) * atom_t.velocityNorm
    end

    # Update atom_p
    avePPoint /= N_t
    x_p = η * sum(x_pList)
    pCoordinate = avePPoint + x_p * atom_p.velocityNorm
    DisplaceAtom(atom_p, pCoordinate, simulator)
    velocity = (sqrt(2 * atom_p.mass * atom_p.energy) * atom_p.velocityNorm - momentum) / atom_p.mass
    atom_p.velocityNorm = velocity / norm(velocity)
    atom_p.energy = atom_p.energy - sumE_t * η - sum(QList) * η
end 

function Cascade!(atom_p::Atom, simulator::Simulator)
    pAtoms = Vector{Atom}([atom_p])
    while true
        targetsList = Vector{Vector{Atom}}()
        for atom in pAtoms
            targets = ShotTarget(atom, simulator)
            for pAtom in pAtoms
                filter!(t->t.index != pAtom.index, targets)
            end
            push!(targetsList, targets)
        end
        UniqueTargets!(targetsList, pAtoms)
        nextPAtoms = Vector{Atom}()
        for (pAtom, targets) in zip(pAtoms, targetsList)
            Collision!(pAtom, targets, simulator)
            for target in targets
                if target.energy > target.dte
                    target.energy -= target.dte
                    if target.energy > simulator.constants.stopEnergy
                        push!(nextPAtoms, target)
                    else
                        Place!(target, simulator)
                    end
                else
                    Recover!(target, simulator)
                end
            end
            if pAtom.energy > simulator.constants.stopEnergy
                push!(nextPAtoms, pAtom)
            else
                Place!(pAtom, simulator)
            end
        end
        if length(nextPAtoms) == 0
            break
        end
        pAtoms = nextPAtoms
    end
end

function UniqueTargets!(targetsList::Vector{Vector{Atom}}, pAtoms::Vector{Atom})
    targetToListDict = Dict{Int64, Vector{Int64}}()
    for (i, targets) in enumerate(targetsList)
        for target in targets
            push!(targetToListDict[target.index], i)
        end
    end
    for (targetIndex, targetsListIndex) in targetToListDict
        if length(targetsListIndex) > 1
            minEnergy = Float(Inf) 
            minArg = 0
            for index in targetsListIndex
                energy = pAtoms[index].energy
                if energy < minEnergy
                    minEnergy = energy
                    minArg = index
                end
            end
            for index in targetsListIndex
                if index != minArg 
                    filter!(t -> t.index != targetIndex, targetsList[index])
                end
            end
        end
    end
end
