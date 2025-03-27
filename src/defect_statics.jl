function CountVacancies(simulator::Simulator)
    nVacancies = 0
    for latticePoint in simulator.latticePoints
        if latticePoint.atomIndex == -1
            nVacancies += 1
        end
    end
    return nVacancies
end

