function GetDTE(atom::Atom, simulator::Simulator)
    return atom.dte
end

function GetBDE(atom::Atom, simulator::Simulator)  # BDE: binding energy
    return atom.bde
end

function RandomPointInCircle(radius::Float64=3.0)
    θ = 2π * rand()  
    r = sqrt(rand()) * radius  
    x = r * cos(θ)
    y = r * sin(θ)
    return [x, y, 0.0]
end  