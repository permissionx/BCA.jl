using .BCA


function ShotTarget(atom::Atom, simulator::Simulator)
    cellGrid = simulator.cellGrid
    periodic = simulator.periodic       
    cell = cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    accCrossFlag = Vector{Int64}([0,0,0])
    while true
        targets = GetTargetsFromNeighbor(atom, cell, simulator)
        if length(targets) > 0
            for cell in cellGrid.cells
                cell.isExplored = false
            end
            return targets
        else
            dimension, direction = AtomOutFaceDimension(atom, cell, accCrossFlag, simulator.box)
            neighborIndex = Vector{Int64}([0,0,0])
            neighborIndex[dimension] = direction == 1 ? -1 : 1
            neighborInfo = cell.neighborCellsInfo[neighborIndex]
            accCrossFlag += neighborInfo.cross
            if !periodic[dimension]
                if neighborInfo.cross[dimension] != 0
                    for cell in cellGrid.cells
                        cell.isExplored = false
                    end
                    return Vector{Atom}() # means find nothing  
                end
            end 
            index = neighborInfo.index
            cell = cellGrid.cells[index[1], index[2], index[3]]
        end
    end
end


function AtomOutFaceDimension(atom::Atom, cell::GridCell, crossFlag::Vector{Int64}, box::Box)
    coordinate = copy(atom.coordinate)
    for d in 1:3
        if crossFlag[d] != 0
            coordinate[d] = coordinate[d] - crossFlag[d] * box.vectors[d, d]
        end
    end 
    for d in 1:3
        if atom.velocityDirection[d] >= 0
            rangeIndex = 2
        else
            rangeIndex = 1
        end
        faceCoordinate = cell.ranges[d, rangeIndex]
        t = (faceCoordinate - coordinate[d]) / atom.velocityDirection[d]
        elseDs = [ed for ed in 1:3 if ed != d]
        allInRange = true
        for elseD in elseDs
            crossCoord = coordinate[elseD] + atom.velocityDirection[elseD] * t
            if !(cell.ranges[elseD, 1] <= crossCoord <= cell.ranges[elseD, 2])
                allInRange = false
                break
            end
        end
        if allInRange
            return d, rangeIndex
        end
    end
    error("Out face not found\n $(atom) \n $(cell)")
end


function Collision!(atom_p::Atom, atoms_t::Vector{Atom}, simulator::Simulator)
    N_t = length(atoms_t)
    cellGrid = simulator.cellGrid
    tanφList = Vector{Float64}(undef, N_t)
    tanψList = Vector{Float64}(undef, N_t)
    E_tList = Vector{Float64}(undef, N_t)
    E_pList = Vector{Float64}(undef, N_t)
    x_pList = Vector{Float64}(undef, N_t)
    x_tList = Vector{Float64}(undef, N_t)
    QList = Vector{Float64}(undef, N_t)
    pL = 0.0
    for atom_t in atoms_t
        l = atom_t.pL[atom_p.index]
        if l > simulator.constants.pLMax 
            l = simulator.constants.pLMax
        end
        pL += l
    end
    pL /= N_t   
    for (i, atom_t) in enumerate(atoms_t)
        p = atom_t.pValue[atom_p.index]
        N = cellGrid.cells[atom_t.cellIndex[1], atom_t.cellIndex[2], atom_t.cellIndex[3]].atomicDensity 
        tanφList[i], tanψList[i], E_tList[i], E_pList[i], x_pList[i], x_tList[i], QList[i] = CollisionParams(atom_p, atom_t, p, p * p, pL, N, simulator)
        #if atom_t.index == 441
            #println("DEBUG:\n tanφ: $(tanφList[i])\n tanψ: $(tanψList[i])\n E_t: $(E_tList[i])\n E_p: $(E_pList[i])\n x_p: $(x_pList[i])\n x_t: $(x_tList[i])\n Q: $(QList[i])")
        #end
    end
    sumE_t = sum(E_tList)
    η = N_t * atom_p.energy / (N_t * atom_p.energy + (N_t - 1) * sumE_t)
    # Update atoms_t (target atoms)         
    avePPoint = Vector{Float64}([0.0,0.0,0.0])
    momentum = Vector{Float64}([0.0,0.0,0.0])
    for (i, atom_t) in enumerate(atoms_t)
        tCoordinate = atom_t.coordinate + x_tList[i] * η * atom_p.velocityDirection
        DisplaceAtom(atom_t, tCoordinate, simulator)
        if atom_t.pValue[atom_p.index] != 0
            velocityDirectionTmp = -atom_t.pVector[atom_p.index] / atom_t.pValue[atom_p.index] * tanψList[i] + atom_p.velocityDirection
        else
            velocityDirectionTmp = atom_p.velocityDirection
        end   
        SetVelocityDirection!(atom_t, velocityDirectionTmp)
        SetEnergy!(atom_t, E_tList[i] * η)
        avePPoint += atom_t.pPoint[atom_p.index]
        momentum += sqrt(2 * atom_t.mass * atom_t.energy) * atom_t.velocityDirection
    end
    # Update atom_p
    avePPoint /= N_t
    x_p = η * sum(x_pList)
    pCoordinate = avePPoint + x_p * atom_p.velocityDirection
    DisplaceAtom(atom_p, pCoordinate, simulator)
    velocity = (sqrt(2 * atom_p.mass * atom_p.energy) * atom_p.velocityDirection - momentum) / atom_p.mass
    SetVelocityDirection!(atom_p, velocity)
    SetEnergy!(atom_p, atom_p.energy - (sumE_t + sum(QList)) * η)
end 


function Cascade!(atom_p::Atom, simulator::Simulator)
    pAtoms = Vector{Atom}([atom_p])
    nStep = 1
    if isDumpInCascade
        Dump(simulator, dumpName, nStep, false)
    end
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
            if length(targets) == 0
                delete!(simulator, pAtom)
                continue
            end
            Collision!(pAtom, targets, simulator)
            for target in targets
                if target.energy > GetDTE(target, simulator)        
                    SetEnergy!(target, target.energy - GetBDE(target, simulator)) # BDE: bonding energy
                    if target.latticePointIndex != -1
                        simulator.latticePoints[target.latticePointIndex].atomIndex = -1
                        target.latticePointIndex = -1
                    end
                    if target.energy > simulator.constants.stopEnergy
                        push!(nextPAtoms, target)
                        if simulator.isStore && target.index <= simulator.atomNumberWhenStore
                        # Store the displaced atom index
                            push!(simulator.displacedAtoms, target.index)
                        end
                    else
                        Stop!(target, simulator)
                    end
                else
                    if target.latticePointIndex != -1   
                        SetEnergy!(target, 0.0)
                        DisplaceAtom(target, simulator.latticePoints[target.latticePointIndex].coordinate, simulator)
                    end
                end 
            end
            if pAtom.energy > simulator.constants.stopEnergy
                push!(nextPAtoms, pAtom)
            else
                Stop!(pAtom, simulator)
            end
        end
        nStep += 1
        for targets in targetsList
            for target in targets
                EmptyP!(target)
            end
        end 
        if isDumpInCascade
            Dump(simulator, dumpName, nStep, true)
        end
        if isLog
            println(nStep)
        end
        if length(nextPAtoms) > 0
            pAtoms = nextPAtoms
        else
            break
        end
    end
end


function UniqueTargets!(targetsList::Vector{Vector{Atom}}, pAtoms::Vector{Atom})
    targetToListDict = Dict{Int64, Vector{Int64}}()
    for (i, targets) in enumerate(targetsList)
        for target in targets
            try 
                push!(targetToListDict[target.index], i)
            catch 
                targetToListDict[target.index] = Vector{Int64}([i])
            end
        end
    end
    for (targetIndex, targetsListIndex) in targetToListDict
        if length(targetsListIndex) > 1
            maxEnergy = -1.0
            maxArg = 0
            for index in targetsListIndex
                energy = pAtoms[index].energy
                if energy > maxEnergy
                    maxEnergy = energy
                    maxArg = index
                end
            end
            for index in targetsListIndex
                if index != maxArg 
                    filter!(t -> t.index != targetIndex, targetsList[index])
                end
            end
        end
    end
end
