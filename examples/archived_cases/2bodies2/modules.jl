function GetDTE(atom::Atom, simulator::Simulator)
    return atom.dte
end

function GetBDE(atom::Atom, simulator::Simulator)  # BDE: binding energy
    return atom.bde
end

function random_point_in_circle(radius::Float64=3.0)
    θ = 2π * rand()  # 随机角度 (0 到 2π)
    r = sqrt(rand()) * radius  # 随机半径，满足面积均匀分布
    x = r * cos(θ)
    y = r * sin(θ)
    return (x, y)
end