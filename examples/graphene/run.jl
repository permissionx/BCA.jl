include("../../src/main.jl")

include("modules.jl")
# Box
a = 1.39667
b = 3.35
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
boxSizes = [10, 20, 10]
inputGridVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a] 
periodic = [true, true, false]

latticeRanges = [0 10; 0 20; 5 6]   
basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0]
basisTypes = [1, 1, 1, 1]

# Parameters
pMax = a/2
qMax = 3.0
stopEnergy = 0.01
vacancyRecoverDistance_squared = 1.3
pLMax = 2.0
dumpName = "graphene.dump"
isDumpInCascade = false
isLog = false

typeDict = Dict(
    1 => (radius = 1.0, mass = 12.0, Z = 6.0, dte = 0.10, bde = 0.10),  # Carbon
    2 => (radius = 1.0, mass = 132.0, Z = 54.0, dte = 0.10, bde = 0.10) # Xenon
)

# Initialize
parameters = Parameters(pMax, qMax, stopEnergy, vacancyRecoverDistance_squared, 
                        pLMax, 
                        dumpName, isDumpInCascade, isLog, 
                        typeDict)


simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, periodic, latticeRanges, basis, basisTypes, parameters)  
#checkpoint = Save!(simulator)
ion = Atom(2, [20.9, 24.1, 20.0], parameters)
SetVelocityDirection!(ion, [0.0,0.0,-1.0])
SetEnergy!(ion, 1000.0)
push!(simulator, ion)
Dump(simulator, "run.dump",0,false)
Cascade!(ion, simulator)
#Restore!(simulator, checkpoint)
Dump(simulator, "run.dump",1,true)


