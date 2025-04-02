include("../../src/main.jl")

using Random
include("modules.jl")
# Parameters
pMax = 3.0 
stopEnergy = 1.0
vacancyRecoverDistance_squared = 1.3
pLMax = 2.0
dumpName = "2bodies.dump"
isDumpInCascade = false
isLog = false

typeDict = Dict(
    1 => (radius = 1.0, mass = 4.0, Z = 2.0, dte = 12.1, bde = 12.1),  # Helium
    2 => (radius = 1.0, mass = 12.0, Z = 6.0, dte = 12.1, bde = 12.1) # Carbon
)

# Box
a = 1.39667
b = 3.35
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
boxSizes = [10, 20, 10]
inputGridVectors = [3.1 0.0 0.0; 0.0 3.1 0.0; 0.0 0.0 3.1]
periodic = [true, true, false]


# Initialize
parameters = Parameters(pMax,stopEnergy, vacancyRecoverDistance_squared, 
                        pLMax, 
                        dumpName, isDumpInCascade, isLog, 
                        typeDict)

box = CreateBoxByPrimaryVectors(primaryVectors, boxSizes)
simulator = Simulator(box, inputGridVectors, periodic, parameters)

#=
atom_t = Atom(1, [10.0, 10.0, 10.0], parameters)
atom_p = Atom(1, [10.0, 10.0, 15.0], parameters)
push!(simulator, atom_t)
push!(simulator, atom_p)
SetVelocityDirection!(atom_p, [0.0, 0.0, -1.0])
SetEnergy!(atom_p, 10.0)
targets = ShotTarget(atom_p, simulator) 
Collision!(atom_p, targets, simulator)
=#


# process
open("log.csv","w") do f
    write(f, "i init_energy atom_p.energy atom_t1.energy atom_t2.energy atom_t1_p atom_t1_px atom_t1_py atom_t1_pz atom_t2_p atom_t2_px atom_t2_py atom_t2_pz atom_p_p atom_p_px atom_p_py atom_p_pz\n")
end
Random.seed!(42)

positions = [[13.530000000000001, -47.57941480245123, 0.0], [11.07, -46.15913376357209, 0.0], [12.3, -45.448993244132524, 0.0], [13.530000000000001, -46.15913376357209, 0.0], [14.760000000000002, -45.448993244132524, 0.0], [12.3, -44.02871220525338, 0.0], [13.530000000000001, -43.31857168581381, 0.0], [14.760000000000002, -44.02871220525338, 0.0]]
shift = [10, 60, 10]
for i in 1:8
    position = positions[i] + shift
    atom_t = Atom(2, position, parameters)
    push!(simulator, atom_t)
end

atom_p = Atom(1, [13.366703575891517, -45.536729716941934, 10.0] + shift, parameters)
SetEnergy!(atom_p, 3000.0)
SetVelocityDirection!(atom_p, [0.0, 0.0, -1.0])
push!(simulator, atom_p)
Cascade!(atom_p, simulator)

