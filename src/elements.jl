function Element(name::String, dte::Float64, bde::Float64)
    if name == "H"
        return (name=name, radius = 0.53, mass = 1.0, Z = 1.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.05)
    elseif name == "He"
        return (name=name, radius = 0.31, mass = 4.0, Z = 2.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.05)
    elseif name == "Li"
        return (name=name, radius = 1.67, mass = 7.0, Z = 3.0, dte = dte, bde = bde, alpha = 1.5, beta = 0.05)
    elseif name == "Be"
        return (name=name, radius = 1.12, mass = 9.0, Z = 4.0, dte = dte, bde = bde, alpha = 1.5, beta = 0.44)
    elseif name == "B"
        return (name=name, radius = 0.87, mass = 11.0, Z = 5.0, dte = dte, bde = bde, alpha = 1.5, beta = 0.44)
    elseif name == "C"
        return (name=name, radius = 0.67, mass = 12.0, Z = 6.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "N"
        return (name=name, radius = 0.56, mass = 14.0, Z = 7.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "O"
        return (name=name, radius = 0.48, mass = 16.0, Z = 8.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "F"
        return (name=name, radius = 0.42, mass = 19.0, Z = 9.0, dte = dte, bde = bde, alpha = 1.1, beta = 0.15)
    elseif name == "Ne"
        return (name=name, radius = 0.38, mass = 20.0, Z = 10.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Na"
        return (name=name, radius = 1.90, mass = 23.0, Z = 11.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Mg"
        return (name=name, radius = 1.45, mass = 24.0, Z = 12.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Al"
        return (name=name, radius = 1.18, mass = 27.0, Z = 13.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.12)
    elseif name == "Si"
        return (name=name, radius = 1.11, mass = 28.0, Z = 14.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "P"
        return (name=name, radius = 0.98, mass = 31.0, Z = 15.0, dte = dte, bde = bde, alpha = 1.25, beta = 0.44)
    elseif name == "S"
        return (name=name, radius = 0.88, mass = 32.0, Z = 16.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Cl"
        return (name=name, radius = 0.79, mass = 35.0, Z = 17.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ar"
        return (name=name, radius = 0.71, mass = 40.0, Z = 18.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "K"
        return (name=name, radius = 2.43, mass = 39.0, Z = 19.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ca"
        return (name=name, radius = 1.94, mass = 40.0, Z = 20.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Sc"
        return (name=name, radius = 1.84, mass = 45.0, Z = 21.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ti"
        return (name=name, radius = 1.76, mass = 48.0, Z = 22.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "V"
        return (name=name, radius = 1.71, mass = 51.0, Z = 23.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Cr"
        return (name=name, radius = 1.66, mass = 52.0, Z = 24.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Mn"
        return (name=name, radius = 1.61, mass = 55.0, Z = 25.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Fe"
        return (name=name, radius = 1.56, mass = 56.0, Z = 26.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Co"
        return (name=name, radius = 1.52, mass = 59.0, Z = 27.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ni"
        return (name=name, radius = 1.49, mass = 58.0, Z = 28.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Cu"
        return (name=name, radius = 1.45, mass = 63.0, Z = 29.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Zn"
        return (name=name, radius = 1.42, mass = 64.0, Z = 30.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ga"
        return (name=name, radius = 1.36, mass = 69.0, Z = 31.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ge"
        return (name=name, radius = 1.25, mass = 74.0, Z = 32.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "As"
        return (name=name, radius = 1.14, mass = 75.0, Z = 33.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Se"
        return (name=name, radius = 1.03, mass = 80.0, Z = 34.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Br"
        return (name=name, radius = 0.94, mass = 79.0, Z = 35.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Kr"
        return (name=name, radius = 0.88, mass = 84.0, Z = 36.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Rb"
        return (name=name, radius = 2.65, mass = 85.0, Z = 37.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Sr"
        return (name=name, radius = 2.19, mass = 88.0, Z = 38.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Y"
        return (name=name, radius = 2.12, mass = 89.0, Z = 39.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Zr"
        return (name=name, radius = 2.06, mass = 90.0, Z = 40.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Nb"
        return (name=name, radius = 1.98, mass = 93.0, Z = 41.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Mo"
        return (name=name, radius = 1.90, mass = 98.0, Z = 42.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Tc"
        return (name=name, radius = 1.83, mass = 98.0, Z = 43.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ru"
        return (name=name, radius = 1.78, mass = 102.0, Z = 44.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Rh"
        return (name=name, radius = 1.73, mass = 103.0, Z = 45.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Pd"
        return (name=name, radius = 1.69, mass = 106.0, Z = 46.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ag"
        return (name=name, radius = 1.65, mass = 107.0, Z = 47.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Cd"
        return (name=name, radius = 1.61, mass = 114.0, Z = 48.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "In"
        return (name=name, radius = 1.56, mass = 115.0, Z = 49.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Sn"
        return (name=name, radius = 1.45, mass = 120.0, Z = 50.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Sb"
        return (name=name, radius = 1.33, mass = 121.0, Z = 51.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Te"
        return (name=name, radius = 1.23, mass = 130.0, Z = 52.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "I"
        return (name=name, radius = 1.15, mass = 129.0, Z = 53.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Xe"
        return (name=name, radius = 1.08, mass = 132.0, Z = 54.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Cs"
        return (name=name, radius = 2.98, mass = 133.0, Z = 55.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ba"
        return (name=name, radius = 2.53, mass = 138.0, Z = 56.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "La"
        return (name=name, radius = 2.26, mass = 139.0, Z = 57.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ce"
        return (name=name, radius = 2.10, mass = 140.0, Z = 58.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Pr"
        return (name=name, radius = 2.47, mass = 141.0, Z = 59.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Nd"
        return (name=name, radius = 2.06, mass = 142.0, Z = 60.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Pm"
        return (name=name, radius = 2.05, mass = 145.0, Z = 61.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Sm"
        return (name=name, radius = 2.38, mass = 152.0, Z = 62.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Eu"
        return (name=name, radius = 2.31, mass = 153.0, Z = 63.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Gd"
        return (name=name, radius = 2.33, mass = 158.0, Z = 64.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Tb"
        return (name=name, radius = 2.25, mass = 159.0, Z = 65.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Dy"
        return (name=name, radius = 2.28, mass = 164.0, Z = 66.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ho"
        return (name=name, radius = 2.26, mass = 165.0, Z = 67.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Er"
        return (name=name, radius = 2.26, mass = 170.0, Z = 68.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Tm"
        return (name=name, radius = 2.22, mass = 169.0, Z = 69.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Yb"
        return (name=name, radius = 2.22, mass = 174.0, Z = 70.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Lu"
        return (name=name, radius = 2.17, mass = 176.0, Z = 71.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Hf"
        return (name=name, radius = 2.08, mass = 180.0, Z = 72.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ta"
        return (name=name, radius = 2.00, mass = 181.0, Z = 73.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "W"
        return (name=name, radius = 1.93, mass = 184.0, Z = 74.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Re"
        return (name=name, radius = 1.88, mass = 187.0, Z = 75.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Os"
        return (name=name, radius = 1.85, mass = 192.0, Z = 76.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ir"
        return (name=name, radius = 1.80, mass = 191.0, Z = 77.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Pt"
        return (name=name, radius = 1.77, mass = 195.0, Z = 78.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Au"
        return (name=name, radius = 1.74, mass = 197.0, Z = 79.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Hg"
        return (name=name, radius = 1.71, mass = 202.0, Z = 80.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Tl"
        return (name=name, radius = 1.56, mass = 205.0, Z = 81.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Pb"
        return (name=name, radius = 1.54, mass = 208.0, Z = 82.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Bi"
        return (name=name, radius = 1.43, mass = 209.0, Z = 83.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Po"
        return (name=name, radius = 1.35, mass = 209.0, Z = 84.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "At"
        return (name=name, radius = 1.27, mass = 210.0, Z = 85.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Rn"
        return (name=name, radius = 1.20, mass = 222.0, Z = 86.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    else
        error("Unknown element: $name")
    end
end
