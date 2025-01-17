using LinearAlgebra
import Base: push!
include("types.jl")
include("bca.jl")  # In namespace BCA: BCA->(QLoss, Constants)
include("io.jl")
include("entities.jl")
include("dynamics.jl")


function main_1()
    primaryVectors = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    boxSizes = [10, 10, 15]
    inputGridVectors = [3.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 3.0]
    periodic = [true, true, false] 
    latticeRanges = [0 10;0 10;0 10]   
    basis = [0.0 0.0 0.0; 0.5 0.5 0.5]
    basisTypes = [1, 1]
    simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, periodic, latticeRanges, basis, basisTypes)
    Dump(simulator, "../output/test_bcc.dump", 0, false)
    return simulator
end


function main_2()
    # 2D Graphene
    a = 1.39667
    primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 2.0]
    boxSizes = [10, 20, 10]
    inputGridVectors = [a/2 0.0 0.0; 0.0 a/2 0.0; 0.0 0.0 a/2]
    periodic = [true, true, false] 
    latticeRanges = [0 10; 0 20; 5 6]   
    basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0]
    basisTypes = [1, 1, 1, 1]
    simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, periodic, latticeRanges, basis, basisTypes)    
    ion = Atom(1, [20.9, 24.1, 15.0])
    SetvelocityDirection!(ion, [6.0, 5.0, -5.0])
    SetEnergy!(ion, 10.0)
    push!(simulator, ion)
    Dump(simulator, "../output/test_graphene.dump", 0, false)
    target = ShotTarget(ion, simulator)
    return simulator, target
end


simulator, target = main_2()

