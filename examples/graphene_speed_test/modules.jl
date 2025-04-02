

function RandomPointInCircle(radius::Float64=3.0)
    θ = 2π * rand()  
    r = sqrt(rand()) * radius  
    x = r * cos(θ)
    y = r * sin(θ)
    return [x, y, 0.0]
end  

function RandomInAnUnitGrapheneCell(a::Float64)
    X = a * 3
    Y = sqrt(3) * a
    x = rand() * X
    y = rand() * Y
    return [x, y, 0.0]
end  

function RandomVectorInUnitSphere(θrange::Float64)
    θ = θrange * rand()
    φ = 2π * rand()
    x = sin(θ) * cos(φ)
    y = sin(θ) * sin(φ)
    z = cos(θ)
    return [x, y, -z]
end



function CheckLatticePoint(simulator::Simulator)
    for i in 1:simulator.numberOfAtoms
        flag1 = simulator.atoms[i].latticePointIndex == simulator.latticePoints[i].atomIndex
        flag2 = i == simulator.atoms[i].latticePointIndex
        if !(flag1 && flag2)
            error("error: ", i, flag1, flag2)
        end
    end
end


function CheckCellDensity(simulator::Simulator)
    volume = simulator.cellGrid.cellVolume
    cells = simulator.cellGrid.cells
    for cell in cells
        density = length(cell.atoms) / volume
        if cell.atomicDensity != density
            error("error: ", cell.atomicDensity, density)
        end
    end
end

function CompareCellAtoms(simulator::Simulator, reference::Simulator, i::Int64)
    cells = simulator.cellGrid.cells
    referenceCells = reference.cellGrid.cells
    error_cell = Vector{Vector{Int64}}()
    error_ref = Vector{Vector{Int64}}()
    flag = true
    for (cell,refCell) in zip(cells, referenceCells)
        flag1 = Set(cell.atoms) == Set(refCell.atoms)
        if !flag1
            push!(error_cell, cell.atoms)
            push!(error_ref, refCell.atoms) 
            flag = false
        end
    end
    if !flag
        error("error: $(error_cell) $(error_ref) $(i)")
    end
end

