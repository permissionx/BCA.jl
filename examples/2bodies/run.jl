include("../../src/main.jl")


include("modules.jl")
# Parameters
pMax = 3.0 
qMax = 3.0
stopEnergy = 1.0
vacancyRecoverDistance_squared = 1.3
pLMax = 2.0
dumpName = "2bodies.dump"
isDumpInCascade = false
isLog = false

typeDict = Dict(
    1 => (radius = 1.0, mass = 12.0, Z = 6.0, dte = 12.1, bde = 12.1),  # Carbon
    2 => (radius = 1.0, mass = 132.0, Z = 54.0, dte = 12.1, bde = 12.1) # Xenon
)

# Box
a = 1.39667
b = 3.35
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
boxSizes = [10, 20, 10]
inputGridVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]
periodic = [true, true, false]


# Initialize
parameters = Parameters(pMax, qMax, stopEnergy, vacancyRecoverDistance_squared, 
                        pLMax, 
                        dumpName, isDumpInCascade, isLog, 
                        typeDict)

box = CreateBoxByPrimaryVectors(primaryVectors, boxSizes)
simulator = Simulator(box, inputGridVectors, periodic, parameters)

# process
atom_t1 = Atom(1, [9.9, 9.9, 10.0], parameters)
push!(simulator, atom_t1)
atom_t2 = Atom(1, [9.0, 9.0, 10.0], parameters)
push!(simulator, atom_t2)
atom_p = Atom(1, [10.0, 10.0, 15.0], parameters)
SetVelocityDirection!(atom_p, [0.0, 0.0, -1.0])
SetEnergy!(atom_p, 10.0)
push!(simulator, atom_p)
targets = ShotTarget(atom_p, simulator) 
Collision!(atom_p, targets, simulator)

