include("../../src/bca.jl")


using Random
using Profile
include("modules.jl")
# Box
a = 1.39667
b = 3.35
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
boxSizes = [10, 20, 10]
inputGridVectors = [a*1.1 0.0 0.0; 0.0 a*1.1 0.0; 0.0 0.0 a*1.1]  # never be same as primaryVectors 
periodic = [true, true, false]

latticeRanges = [0 10; 0 20; 5 6]   
basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0]
basisTypes = [1, 1, 1, 1]

# Parameters
pMax = 1.45
stopEnergy = 1
vacancyRecoverDistance_squared = 1.3
pLMax = 2.0
dumpName = "graphene.dump"
isDumpInCascade = false
isLog = false


typeDict = Dict(
    1 => (radius = 0.67, mass = 12.0, Z = 6.0, dte =22, bde = 22),  # Carbon
    2 => (radius = 0.38, mass = 20.0, Z = 10.0, dte = 22, bde = 22) # Neon
)


# Initialize
parameters = Parameters(pMax, stopEnergy, vacancyRecoverDistance_squared, 
                        pLMax, 
                        dumpName, isDumpInCascade, isLog, 
                        typeDict)

simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, periodic, latticeRanges, basis, basisTypes, parameters)  
Save!(simulator)

Random.seed!(42)
#Dump(simulator, "run.dump",0,false)
open("p.debug.log", "w") do f
    write(f, "nrun,ntargets,px0,py0,pz0,px1,py1,pz1\n")
end

Dump(simulator, "test", 1, false)

Profile.clear()
@profile begin
for i in 1:1000
    if i % 100 == 0
        println("Irradiation time: ", i)
    end
    simulator.nIrradiation += 1
    ionPosition = RandomPointInCircle(10.) + [20., 24., 20.]
    ion = Atom(2, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.0,0.0,-1.0])
    SetEnergy!(ion, 1000.0)
    push!(simulator, ion)
    Cascade!(ion, simulator)
    Restore!(simulator)
end
end
Profile.print(format=:flat) 


