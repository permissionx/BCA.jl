using .BCA

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
