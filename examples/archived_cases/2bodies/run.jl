include("../../src/main.jl")

using Random
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
parameters = Parameters(pMax, qMax, stopEnergy, vacancyRecoverDistance_squared, 
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

#for i in 1:1
i=1
    atom_t1 = Atom(2, [10.5, 10.5, 10.0], parameters)
    push!(simulator, atom_t1)
    atom_t2 = Atom(2, [9.5, 9.5, 10.0], parameters)
    push!(simulator, atom_t2)
    x, y = random_point_in_circle(2.0)
    atom_p = Atom(1, [x+10.0, y+10.0, 15.0], parameters)
    SetVelocityDirection!(atom_p, [0.0, 0.0, -1.0])
    initEnergy = rand() * 10.0
    SetEnergy!(atom_p, initEnergy)
    push!(simulator, atom_p)
    targets = ShotTarget(atom_p, simulator) 
    Collision!(atom_p, targets, simulator)
    atom_t1_p = sqrt(atom_t1.energy * atom_t1.mass * 2) 
    atom_t1_px = atom_t1_p * atom_t1.velocityDirection[1]
    atom_t1_py = atom_t1_p * atom_t1.velocityDirection[2]
    atom_t1_pz = atom_t1_p * atom_t1.velocityDirection[3]
    atom_t2_p = sqrt(atom_t2.energy * atom_t2.mass * 2) 
    atom_t2_px = atom_t2_p * atom_t2.velocityDirection[1]
    atom_t2_py = atom_t2_p * atom_t2.velocityDirection[2]
    atom_t2_pz = atom_t2_p * atom_t2.velocityDirection[3]
    atom_p_p = sqrt(atom_p.energy * atom_p.mass * 2) 
    atom_p_px = atom_p_p * atom_p.velocityDirection[1]
    atom_p_py = atom_p_p * atom_p.velocityDirection[2]
    atom_p_pz = atom_p_p * atom_p.velocityDirection[3]
    open("log.csv","a") do f
        write(f, "$i $(initEnergy) $(atom_p.energy) $(atom_t1.energy) $(atom_t2.energy) $(atom_t1_p) $(atom_t1_px) $(atom_t1_py) $(atom_t1_pz) $(atom_t2_p) $(atom_t2_px) $(atom_t2_py) $(atom_t2_pz) $(atom_p_p) $(atom_p_px) $(atom_p_py) $(atom_p_pz)\n")
    end

    if atom_t1.isAlive == true
        delete!(simulator, atom_t1)
    end
    if atom_t2.isAlive == true
        delete!(simulator, atom_t2)
    end     
    if atom_p.isAlive == true
        delete!(simulator, atom_p)
    end
#end 

