# using BCA.jl
include("../../src/main.jl")
include("modules.jl")


using Random


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
θτRepository = "../../thetatau_repository/"
pMax = 1.45
vacancyRecoverDistance = 1.3
# for Xe ion 
typeDict = Dict(
    1 => Element("C", 22.0, 22.0),  
    2 => Element("Xe", 22.0, 22.0)  
)
parameters = Parameters(θτRepository, pMax, vacancyRecoverDistance, typeDict)


# Initialize
simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, periodic, latticeRanges, basis, basisTypes, parameters)  
Save!(simulator)

#Dump(simulator, "run.dump",0,false)




function Irradiation(simulator::Simulator, energy::Float64)
    Restore!(simulator)
    simulator.nIrradiation += 1
    ionPosition = RandomInAnUnitGrapheneCell(1.39667) + [18.855045, 22.981482313368623, 20]
    ion = Atom(2, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.,0.,-1.])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)
end


Random.seed!(42)
Irradiation(simulator, 1.0) # warm up
open("time.csv","w") do f
    write(f,"energy,time\n")
end
for energy_power in 1:6
    energy = 10.0^energy_power
    time_start = time()
    @showprogress desc="Irradiation of energy power: $(energy_power)" for i in 1:1000
        Irradiation(simulator, energy)
    end
    time_end = time()
    open("time.csv","a") do f
        write(f,"$(energy),$(time_end-time_start)\n")
    end
end


