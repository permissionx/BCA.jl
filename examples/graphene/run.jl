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
stopEnergy = 10
vacancyRecoverDistance_squared = 1.3
pLMax = 2.0
dumpName = "graphene.dump"
isDumpInCascade = false
isLog = false


typeDict = Dict(
    1 => Element("C", 22.0, 22.0),  # Carbon
    2 => Element("Ne", 20.0, 20.0) # Neon
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




function Irradiation(simulator::Simulator, energy::Float64, vacancyNumbers::Vector{Int64}, i::Int64)
    #if i % 100 == 0
        #println("Irradiation time: ", i)
    #end
    simulator.nIrradiation += 1
    ionPosition = RandomInAnUnitGrapheneCell(1.39667) + [18.855045, 22.981482313368623, 20]
    ion = Atom(2, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)
    nVacancies = CountVacancies(simulator)
    vacancyNumbers[i] = nVacancies
    Restore!(simulator)
end


computerNumberPerEnergy = 100000
vacancy_data = Dict{Int64, Vector{Int64}}()
for energy_order in 0:65
    energy = 10*1.3^energy_order
    vacancyNumbers = Vector{Int64}(undef, computerNumberPerEnergy)
    for i in 1:computerNumberPerEnergy
        if i%10000 == 0
            println(energy_order, " ", i)
        end
        Irradiation(simulator, energy, vacancyNumbers, i)
    end
    vacancy_data[energy_orde]
end


println("Outputting data...")
open("vacancy.csv", "w") do f
    write(f, "energy_order,energy,nVacancies\n")
    for energy_order in sort(collect(keys(vacancy_data)))
        for nVacancies in vacancy_data[energy_order]
            write(f, "$(energy_order),$(2^energy_order),$(nVacancies)\n")
        end
    end
end

