# using BCA.jl
include("../../src/main.jl")
include("modules.jl")


using Random
using ProgressMeter

# Box
a = 1.42
b = 3.35
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
boxSizes = [10, 20, 10]
inputGridVectors = [a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]  # never be same as primaryVectors 
periodic = [true, true, false]

latticeRanges = [0 10; 0 20; 5 6]   
basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0]
basisTypes = [1, 1, 1, 1]

# Parameters
pMax = 1.45
stopEnergy = 0.1
vacancyRecoverDistance_squared = 1.3
pLMax = 2.0
isDumpInCascade = false
isLog = false


# for Ne ion 
typeDict = Dict(
    1 => Element("C", 22.0, 22.0),  
    2 => Element("Ne", 22.0, 22.0)  
)


# Initialize
println("Initializing simulator...")
parameters = Parameters(pMax, stopEnergy, vacancyRecoverDistance_squared, 
                        pLMax, 
                        isDumpInCascade, isLog, 
                        typeDict)

simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, periodic, latticeRanges, basis, basisTypes, parameters)  
Save!(simulator)

Random.seed!(42)
#Dump(simulator, "run.dump",0,false)




function Irradiation(simulator::Simulator, energy::Float64)
    #if i % 100 == 0
        #println("Irradiation time: ", i)
    #end
    Restore!(simulator)
    simulator.nIrradiation += 1
    ionPosition = RandomInAnUnitGrapheneCell(1.39667) + [18.855045, 22.981482313368623, 20]
    ion = Atom(2, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)
    nVacancies = CountVacancies(simulator)
    return nVacancies
end


println("Comiling...")
#nVacancies = Irradiation(simulator, 1000.)