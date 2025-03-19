function GetDTE(atom::Atom, simulator::Simulator)
    return atom.dte
end

function GetBDE(atom::Atom, simulator::Simulator)  # BDE: binding energy
    return atom.bde
end