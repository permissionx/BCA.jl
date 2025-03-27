include("../../src/bca.jl")


using Random

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
    1 => (radius = 0.67, mass = 12.0, Z = 6.0, dte =22, bde = 22, alpha = 1.0, beta = 0.44),  # Carbon
    2 => (radius = 0.38, mass = 20.0, Z = 10.0, dte = 22, bde = 22, alpha = 1.0, beta = 0.44) # Neon
)


# Initialize
parameters = Parameters(pMax, stopEnergy, vacancyRecoverDistance_squared, 
                        pLMax, 
                        isDumpInCascade, isLog, 
                        typeDict)

simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, periodic, latticeRanges, basis, basisTypes, parameters)  
Save!(simulator)

Random.seed!(42)
#Dump(simulator, "run.dump",0,false)

open("p.debug.log", "w") do f
    write(f, "nrun,ntargets,px0,py0,pz0,px1,py1,pz1\n")
end


function Irradiation(simulator::Simulator, energy::Float64)
    #if i % 100 == 0
        #println("Irradiation time: ", i)
    #end
    simulator.nIrradiation += 1
    ionPosition = RandomInAnUnitGrapheneCell(1.39667) + [18.855045, 22.981482313368623, 16.75]
    ion = Atom(2, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)
    Restore!(simulator)
end



vacancy_data = Dict{Int64, Vector{Int64}}()
for energy_order in 0:20
    energy = 2.0^energy_order
    vacancy_data[energy_order] = Vector{Int64}(undef, 10000)
    for i in 1:10000
        if i%1000 == 0
            println(energy_order, " ", i)
        end
        Irradiation(simulator, energy)
        nVacancies = CountVacancies(simulator)
        vacancy_data[energy_order][i] = nVacancies
    end
end

println("Outputting data...")
open("vacancy.csv", "w") do f
    write(f, "energy_order,energy,nVacancies\n")
    for energy_order in keys(vacancy_data)
        for nVacancies in vacancy_data[energy_order]
            write(f, "$(energy_order),$(2^energy_order),$(nVacancies)\n")
        end
    end
end

