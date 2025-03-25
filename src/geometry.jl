import Base: push!
import Base: delete!
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


function Atom(type::Int64, coordinate::Vector{Float64}, parameters::Parameters)
    index = 0
    isAlive = true
    cellIndex = Vector{Int64}()
    velocityDirection = Float64[0.0,0.0,0.0]  # length： 0 or one
    energy = 0.0
    radius, mass, Z, dte, bde = TypeToProperties(type, parameters.typeDict)
    pValue = Dict{Int64, Float64}() # key is the index of the atom_p
    pVector = Dict{Int64, Vector{Float64}}()
    pPoint = Dict{Int64, Vector{Float64}}()  
    lastTargets = Vector{Int64}()
    pL = Dict{Int64, Float64}()
    latticePointIndex = -1
    return Atom(index, isAlive, type, coordinate, cellIndex, 
                radius, mass, velocityDirection, energy, Z, 
                dte, bde,
                pValue, pVector, pPoint, pL, lastTargets,
                latticePointIndex)
end


function TypeToProperties(type::Int64, 
    typeDict::Dict{Int, NamedTuple{(:radius, :mass, :Z, :dte, :bde), 
              Tuple{Float64, Float64, Float64, Float64, Float64}}})
    if haskey(typeDict, type)
        props = typeDict[type]
        return props.radius, props.mass, props.Z, props.dte, props.bde
    else
        error("Unknown atom type: $type")
    end 
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
                cells[x, y, z] = GridCell(Vector{Int64}([x,y,z]),Vector{Atom}(), Vector{LatticePoint}(), 
                                          ranges, centerCoordinate, 
                                          Dict{Vector{Int64}, Vector{Int64}}(), false, 0.0)
            end
        end    
    end
    cellVolume = vectors[1,1] * vectors[2,2] * vectors[3,3]
    cellGrid = CellGrid(cells, vectors, sizes, cellVolume) 
    for cell in cellGrid.cells
        IterPushCellNeighbors!(cellGrid, cell, Vector{Int64}(), Vector{Int64}(), Vector{Int64}(), 1)
    end
    return cellGrid
end


function InitConstants(parameters::Parameters)
    return Constants(parameters.pMax, parameters.stopEnergy,
                     parameters.vacancyRecoverDistance_squared, parameters.pLMax, 
                     parameters.dumpName, parameters.isDumpInCascade, parameters.isLog)
end


function InitConstantsByType(typeDict::Dict{Int, 
                             NamedTuple{(:radius, :mass, :Z, :dte, :bde), 
                                        Tuple{Float64, Float64, Float64, Float64, Float64}}},
                             constants::Constants) 
    V_upterm = Dict{Vector{Int64}, Float64}()
    a_U = Dict{Vector{Int64}, Float64}()
    E_m = Dict{Int64, Float64}()
    S_e_upTerm = Dict{Vector{Int64}, Float64}()
    S_e_downTerm = Dict{Vector{Int64}, Float64}()
    x_nl = Dict{Vector{Int64}, Float64}()
    a = Dict{Vector{Int64}, Float64}()
    Q_nl = Dict{Vector{Int64}, Float64}()  
    Q_loc = Dict{Vector{Int64}, Float64}()
    types = keys(typeDict)
    qMax_squared = Dict{Vector{Int64}, Float64}()
    for p in types
        for t in types
            radius_p, mass_p, Z_p, _, _ = TypeToProperties(p, typeDict)
            radius_t, _, Z_t, _, _ = TypeToProperties(t, typeDict)
            V_upterm[[p,t]] = BCA.ConstantFunctions.V_upterm(Z_p, Z_t)
            a_U[[p,t]] = BCA.ConstantFunctions.a_U(Z_p, Z_t)
            E_m[p] = BCA.ConstantFunctions.E_m(Z_p, mass_p)
            S_e_upTerm[[p,t]] = BCA.ConstantFunctions.S_e_upTerm(p, Z_p, Z_t, mass_p)
            x_nl[[p,t]] = BCA.ConstantFunctions.x_nl(p, Z_p, Z_t)
            a[[p,t]] = BCA.ConstantFunctions.a(Z_p, Z_t)
            Q_nl[[p,t]] = BCA.ConstantFunctions.Q_nl(Z_p, Z_t, constants.pMax)
            Q_loc[[p,t]] = BCA.ConstantFunctions.Q_loc(Z_p, Z_t)
            qMax_squared[[p,t]] = (radius_p + radius_t) * (radius_p + radius_t)
        end
    end
    return ConstantsByType(V_upterm, a_U, E_m, S_e_upTerm, S_e_downTerm, x_nl, a, Q_nl, Q_loc, qMax_squared)
end


function Simulator(box::Box, inputGridVectors::Matrix{Float64}, periodic::Vector{Bool}, parameters::Parameters)
    cellGrid = CreateCellGrid(box, inputGridVectors)
    constants = InitConstants(parameters)
    constantsByType = InitConstantsByType(parameters.typeDict, constants)
    return Simulator(Vector{Atom}(), Vector{LatticePoint}(), 
                     box, cellGrid, periodic, box.isOrthogonal, 0, 0, 
                     constantsByType, constants,
                     false, Vector{Int64}(), 0, 0, 
                     parameters.isDumpInCascade, parameters.dumpName, parameters.isLog)
end 


function Simulator(primaryVectors::Matrix{Float64}, boxSizes::Vector{Int64}, 
                   inputGridVectors::Matrix{Float64},
                   periodic::Vector{Bool}, 
                   latticeRanges::Matrix{Int64}, basis::Matrix{Float64}, basisTypes::Vector{Int64},
                   parameters::Parameters)
    box = CreateBoxByPrimaryVectors(primaryVectors, boxSizes)
    simulator = Simulator(box, inputGridVectors, periodic, parameters)
    for x in latticeRanges[1,1]:latticeRanges[1,2]-1
        for y in latticeRanges[2,1]:latticeRanges[2,2]-1    
            for z in latticeRanges[3,1]:latticeRanges[3,2]-1
                for i in 1:length(basisTypes)
                    reducedCoordinate = Float64[x,y,z] + basis[i, :]
                    coordinate = primaryVectors' * reducedCoordinate
                    atom = Atom(basisTypes[i], coordinate, parameters)
                    push!(simulator, atom)
                    latticePoint = LatticePoint(copy(atom.index), copy(atom.type), 
                                                copy(atom.coordinate), copy(atom.cellIndex),
                                                atom.index)
                    push!(simulator, latticePoint)
                end
            end
        end
    end
    for cell in simulator.cellGrid.cells
        cell.atomicDensity = length(cell.atoms) / simulator.cellGrid.cellVolume
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
    push!(simulator.cellGrid.cells[cellIndex[1], cellIndex[2], cellIndex[3]].atoms, atom.index)
end 


function push!(simulator::Simulator, latticePoint::LatticePoint)
    push!(simulator.latticePoints, latticePoint)
    push!(simulator.cellGrid.cells[latticePoint.cellIndex[1], latticePoint.cellIndex[2], latticePoint.cellIndex[3]].latticePoints, latticePoint.index)
    simulator.atoms[latticePoint.atomIndex].latticePointIndex = latticePoint.index
end 

function push!(cell::GridCell, atom::Atom, simulator::Simulator)
    atom.cellIndex = cell.index[:]
    push!(cell.atoms, atom.index)
    cell.atomicDensity = length(cell.atoms) / simulator.cellGrid.cellVolume
end 


function delete!(simulator::Simulator, atom::Atom)
    originalCell = simulator.cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    delete!(originalCell, atom.index, simulator)
    simulator.numberOfAtoms -= 1
    atom.isAlive = false 
    AddToStore!(atom, simulator)       
end


function delete!(cell::GridCell, atomIndex::Int64, simulator::Simulator)
    if !simulator.atoms[atomIndex].isAlive
        error("Atom $(atomIndex) is not alive when deleting")
    end
    deleteat!(cell.atoms, findfirst(==(atomIndex), cell.atoms)) 
    simulator.atoms[atomIndex].cellIndex = Vector{Int64}([-1, -1, -1])
    cell.atomicDensity = length(cell.atoms) / simulator.cellGrid.cellVolume  
end


function DisplaceAtom!(atom::Atom, newPosition::Vector{Float64}, simulator::Simulator)
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
        ChangeCell!(atom, cellIndex, simulator)
    end
end


function ComputeDistance_squared(coordinate1::Vector{Float64}, coordinate2::Vector{Float64}, crossFlag::Vector{Int64}, box::Box)
    dv = VectorDifference(coordinate1, coordinate2, crossFlag, box)
    distance_squared = dv[1]* dv[1] + dv[2]*dv[2] + dv[3]  * dv[3]
    return distance_squared
end


function ComputeDistance(coordinate1::Vector{Float64}, coordinate2::Vector{Float64}, crossFlag::Vector{Int64}, box::Box)
    return sqrt(ComputeDistance_squared(coordinate1, coordinate2, crossFlag, box))
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
    atom_t.pL[atom_p.index] = t
    atom_t.pPoint[atom_p.index] = atom_p.coordinate + t * atom_p.velocityDirection
    atom_t.pVector[atom_p.index] = atom_t.pPoint[atom_p.index] - atom_t.coordinate
    atom_t.pValue[atom_p.index] = norm(atom_t.pVector[atom_p.index])
end

function GetAllNeighbors(gridCell::GridCell, simulator::Simulator)
    cellGrid = simulator.cellGrid
    neighbors = Vector{Atom}()
    for (_, neighborCellInfo) in gridCell.neighborCellsInfo
        index = neighborCellInfo.index 
        neighborCell = cellGrid.cells[index[1], index[2], index[3]]
        neighborAtomIndex = neighborCell.atoms
        append!(neighbors, [simulator.atoms[index] for index in neighborAtomIndex])
    end
    return neighbors
end


function GetTargetsFromNeighbor(atom::Atom, gridCell::GridCell, simulator::Simulator)
    cellGrid = simulator.cellGrid
    box = simulator.box
    targets = Vector{Atom}()
    infiniteFlag = true
    for (_, neighborCellInfo) in gridCell.neighborCellsInfo
        index = neighborCellInfo.index
        neighborCell = cellGrid.cells[index[1], index[2], index[3]]
        if neighborCell.isExplored
            continue
        end
        infiniteFlag = false
        for neighborAtomIndex in neighborCell.atoms
            neighborAtom = simulator.atoms[neighborAtomIndex]
            #if simulator.nIrradiation == 332 && neighborAtom.index == 296
            #    println(atom.index, " ", neighborAtom.index)
            #    println(neighborAtom)
            #    println(neighborCell.index)
            #    println(simulator.cellGrid.cells[11,21,11].atoms)
            #end
            if ComputeVDistance(atom, neighborAtom, neighborCellInfo.cross, box) > 0
                ComputeP!(atom, neighborAtom, neighborCellInfo.cross, box)
                if neighborAtom.pValue[atom.index] >= simulator.constants.pMax
                    DeleteP!(neighborAtom, atom.index)
                    continue
                end
                matchFlag = true
                for target in targets
                    if !SimultaneousCriteria(atom, neighborAtom, target, neighborCellInfo.cross, simulator)
                        matchFlag = false
                        DeleteP!(neighborAtom, atom.index)
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
    if infiniteFlag
        error("Atom$(atom.index) fly inifinitely")
    end
    return targets
end

function DeleteP!(atom_t::Atom, atom_pIndex::Int64)
    delete!(atom_t.pValue, atom_pIndex)
    delete!(atom_t.pPoint, atom_pIndex)
    delete!(atom_t.pVector, atom_pIndex)
    delete!(atom_t.pL, atom_pIndex)
end

function EmptyP!(atom_t::Atom)
    empty!(atom_t.pValue)
    empty!(atom_t.pPoint)
    empty!(atom_t.pVector)
    empty!(atom_t.pL)
end


function SimultaneousCriteria(atom::Atom, neighborAtom::Atom, addedTarget::Atom, crossFlag::Vector{Int64}, simulator::Simulator)
    box = simulator.box
    qMax_squared = simulator.constantsByType.qMax_squared[[neighborAtom.type, addedTarget.type]]
    pMax = simulator.constants.pMax
    flagP = abs(neighborAtom.pValue[atom.index] - addedTarget.pValue[atom.index]) <= pMax
    flagQ = sum(VectorDifference(neighborAtom.coordinate, addedTarget.coordinate, crossFlag, box) .* atom.velocityDirection) <= sqrt(qMax_squared)
    return flagP && flagQ
end


function ChangeCell!(atom::Atom, nextCellIndex::Vector{Int64}, simulator::Simulator)
    originalCell = simulator.cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    delete!(originalCell, atom.index, simulator)
    nextCell = simulator.cellGrid.cells[nextCellIndex[1], nextCellIndex[2], nextCellIndex[3]]
    push!(nextCell, atom, simulator)
end


function SetVelocityDirection!(atom::Atom, velocity::Vector{Float64})
    atom.velocityDirection = velocity / norm(velocity)
end

function SetEnergy!(atom::Atom, energy::Float64)
    atom.energy = energy
end

function GetNeighborVacancy(atom::Atom, simulator::Simulator)
    cells = simulator.cellGrid.cells        
    cell = cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    nearestVacancyDistance_squared = Inf
    nearestVacancyIndex = -1
    for (_, neighborCellInfo) in cell.neighborCellsInfo
        index = neighborCellInfo.index
        cross = neighborCellInfo.cross
        neighborCell = cells[index[1], index[2], index[3]]
        for latticePointIndex in neighborCell.latticePoints
            latticePoint = simulator.latticePoints[latticePointIndex]
            if latticePoint.atomIndex == -1
                dr2 = ComputeDistance_squared(atom.coordinate, latticePoint.coordinate, cross, simulator.box)
                if dr2 < simulator.constants.vacancyRecoverDistance_squared && dr2 < nearestVacancyDistance_squared
                    nearestVacancyDistance_squared = dr2
                    nearestVacancyIndex = latticePoint.index
                end
            end 
        end
    end
    return nearestVacancyIndex
end


function Stop!(atom::Atom, simulator::Simulator)
    #SetEnergy!(atom, 0.0)
    Recover!(atom, simulator)
end


function Recover!(atom::Atom, simulator::Simulator)
    nearestVacancyIndex = GetNeighborVacancy(atom, simulator)
    if nearestVacancyIndex != -1
        SetOnLatticePoint!(atom, simulator.latticePoints[nearestVacancyIndex], simulator)
    else
        atom.latticePointIndex = -1
        AddToStore!(atom, simulator)
    end
end




function SetOnLatticePoint!(atom::Atom, latticePoint::LatticePoint, simulator::Simulator)
    if atom.latticePointIndex != latticePoint.index  # changed lattice point 
        if atom.index == latticePoint.index 
            DeleteFromStore!(atom, simulator)
        elseif atom.latticePointIndex == atom.index
            AddToStore!(atom, simulator)
        end
        latticePoint.atomIndex = atom.index
        atom.latticePointIndex = latticePoint.index
        atom.coordinate = latticePoint.coordinate[:]
        if atom.cellIndex != latticePoint.cellIndex
            ChangeCell!(atom, latticePoint.cellIndex, simulator)
        end
    else 
        atom.coordinate = latticePoint.coordinate[:]
        latticePoint.atomIndex = atom.index
    end
end


function AddToStore!(atom::Atom, simulator::Simulator)
    if simulator.isStore && atom.index <= simulator.atomNumberWhenStore
        push!(simulator.displacedAtoms, atom.index)
    end
end

function DeleteFromStore!(atom::Atom, simulator::Simulator)
    if simulator.isStore && atom.index <= simulator.atomNumberWhenStore
        deleteat!(simulator.displacedAtoms, findfirst(==(atom.index), simulator.displacedAtoms))
    end
end 

function Restore!(simulator::Simulator)
    if length(simulator.displacedAtoms) != length(Set(simulator.displacedAtoms))
        error("repeat atom index in displacedAtoms")
    end 
    for index in simulator.displacedAtoms
        atom = simulator.atoms[index]
        latticePoint = simulator.latticePoints[atom.index]
        atom.latticePointIndex = atom.index
        latticePoint.atomIndex = atom.index
        atom.coordinate = latticePoint.coordinate[:]
        if atom.isAlive
            ChangeCell!(atom, latticePoint.cellIndex, simulator)
        else
            atom.isAlive = true
            nextCell = simulator.cellGrid.cells[latticePoint.cellIndex[1], latticePoint.cellIndex[2], latticePoint.cellIndex[3]]
            push!(nextCell, atom, simulator)
        end
    end
    maxAtomID = simulator.atomNumberWhenStore
    simulator.atoms = simulator.atoms[1:maxAtomID]
    simulator.maxAtomID = maxAtomID
    simulator.numberOfAtoms  = maxAtomID
    empty!(simulator.displacedAtoms)
end

function Save!(simulator::Simulator)
    for atom in simulator.atoms
        if atom.latticePointIndex == -1
            error("Atom $(atom.index) is not on lattice when stored.")
        end
    end
    simulator.isStore = true
    simulator.atomNumberWhenStore = simulator.numberOfAtoms
end 
