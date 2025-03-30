function Element(name::String, dte::Float64, bde::Float64)
    if name == "H"
        return (radius = 0.53, mass = 1.0, Z = 1.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.05)
    elseif name == "He"
        return (radius = 0.31, mass = 4.0, Z = 2.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.05)
    elseif name == "Li"
        return (radius = 1.67, mass = 7.0, Z = 3.0, dte = dte, bde = bde, alpha = 1.5, beta = 0.05)
    elseif name == "Be"
        return (radius = 1.12, mass = 9.0, Z = 4.0, dte = dte, bde = bde, alpha = 1.5, beta = 0.44)
    elseif name == "B"
        return (radius = 0.87, mass = 11.0, Z = 5.0, dte = dte, bde = bde, alpha = 1.5, beta = 0.44)
    elseif name == "C"
        return (radius = 0.67, mass = 12.0, Z = 6.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "N"
        return (radius = 0.56, mass = 14.0, Z = 7.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "O"
        return (radius = 0.48, mass = 16.0, Z = 8.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "F"
        return (radius = 0.42, mass = 19.0, Z = 9.0, dte = dte, bde = bde, alpha = 1.1, beta = 0.15)
    elseif name == "Ne"
        return (radius = 0.38, mass = 20.0, Z = 10.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Na"
        return (radius = 1.90, mass = 23.0, Z = 11.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Mg"
        return (radius = 1.45, mass = 24.0, Z = 12.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Al"
        return (radius = 1.18, mass = 27.0, Z = 13.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.12)
    elseif name == "Si"
        return (radius = 1.11, mass = 28.0, Z = 14.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "P"
        return (radius = 0.98, mass = 31.0, Z = 15.0, dte = dte, bde = bde, alpha = 1.25, beta = 0.44)
    elseif name == "S"
        return (radius = 0.88, mass = 32.0, Z = 16.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Cl"
        return (radius = 0.79, mass = 35.0, Z = 17.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ar"
        return (radius = 0.71, mass = 40.0, Z = 18.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "K"
        return (radius = 2.43, mass = 39.0, Z = 19.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ca"
        return (radius = 1.94, mass = 40.0, Z = 20.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Sc"
        return (radius = 1.84, mass = 45.0, Z = 21.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ti"
        return (radius = 1.76, mass = 48.0, Z = 22.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "V"
        return (radius = 1.71, mass = 51.0, Z = 23.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Cr"
        return (radius = 1.66, mass = 52.0, Z = 24.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Mn"
        return (radius = 1.61, mass = 55.0, Z = 25.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Fe"
        return (radius = 1.56, mass = 56.0, Z = 26.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Co"
        return (radius = 1.52, mass = 59.0, Z = 27.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ni"
        return (radius = 1.49, mass = 58.0, Z = 28.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Cu"
        return (radius = 1.45, mass = 63.0, Z = 29.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Zn"
        return (radius = 1.42, mass = 64.0, Z = 30.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ga"
        return (radius = 1.36, mass = 69.0, Z = 31.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ge"
        return (radius = 1.25, mass = 74.0, Z = 32.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "As"
        return (radius = 1.14, mass = 75.0, Z = 33.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Se"
        return (radius = 1.03, mass = 80.0, Z = 34.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Br"
        return (radius = 0.94, mass = 79.0, Z = 35.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Kr"
        return (radius = 0.88, mass = 84.0, Z = 36.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Rb"
        return (radius = 2.65, mass = 85.0, Z = 37.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Sr"
        return (radius = 2.19, mass = 88.0, Z = 38.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Y"
        return (radius = 2.12, mass = 89.0, Z = 39.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Zr"
        return (radius = 2.06, mass = 90.0, Z = 40.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Nb"
        return (radius = 1.98, mass = 93.0, Z = 41.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Mo"
        return (radius = 1.90, mass = 98.0, Z = 42.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Tc"
        return (radius = 1.83, mass = 98.0, Z = 43.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ru"
        return (radius = 1.78, mass = 102.0, Z = 44.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Rh"
        return (radius = 1.73, mass = 103.0, Z = 45.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Pd"
        return (radius = 1.69, mass = 106.0, Z = 46.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ag"
        return (radius = 1.65, mass = 107.0, Z = 47.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Cd"
        return (radius = 1.61, mass = 114.0, Z = 48.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "In"
        return (radius = 1.56, mass = 115.0, Z = 49.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Sn"
        return (radius = 1.45, mass = 120.0, Z = 50.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Sb"
        return (radius = 1.33, mass = 121.0, Z = 51.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Te"
        return (radius = 1.23, mass = 130.0, Z = 52.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "I"
        return (radius = 1.15, mass = 129.0, Z = 53.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Xe"
        return (radius = 1.08, mass = 132.0, Z = 54.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Cs"
        return (radius = 2.98, mass = 133.0, Z = 55.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ba"
        return (radius = 2.53, mass = 138.0, Z = 56.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "La"
        return (radius = 2.26, mass = 139.0, Z = 57.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ce"
        return (radius = 2.10, mass = 140.0, Z = 58.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Pr"
        return (radius = 2.47, mass = 141.0, Z = 59.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Nd"
        return (radius = 2.06, mass = 142.0, Z = 60.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Pm"
        return (radius = 2.05, mass = 145.0, Z = 61.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Sm"
        return (radius = 2.38, mass = 152.0, Z = 62.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Eu"
        return (radius = 2.31, mass = 153.0, Z = 63.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Gd"
        return (radius = 2.33, mass = 158.0, Z = 64.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Tb"
        return (radius = 2.25, mass = 159.0, Z = 65.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Dy"
        return (radius = 2.28, mass = 164.0, Z = 66.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ho"
        return (radius = 2.26, mass = 165.0, Z = 67.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Er"
        return (radius = 2.26, mass = 170.0, Z = 68.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Tm"
        return (radius = 2.22, mass = 169.0, Z = 69.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Yb"
        return (radius = 2.22, mass = 174.0, Z = 70.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Lu"
        return (radius = 2.17, mass = 176.0, Z = 71.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Hf"
        return (radius = 2.08, mass = 180.0, Z = 72.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ta"
        return (radius = 2.00, mass = 181.0, Z = 73.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "W"
        return (radius = 1.93, mass = 184.0, Z = 74.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Re"
        return (radius = 1.88, mass = 187.0, Z = 75.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Os"
        return (radius = 1.85, mass = 192.0, Z = 76.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Ir"
        return (radius = 1.80, mass = 191.0, Z = 77.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Pt"
        return (radius = 1.77, mass = 195.0, Z = 78.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Au"
        return (radius = 1.74, mass = 197.0, Z = 79.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Hg"
        return (radius = 1.71, mass = 202.0, Z = 80.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Tl"
        return (radius = 1.56, mass = 205.0, Z = 81.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Pb"
        return (radius = 1.54, mass = 208.0, Z = 82.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Bi"
        return (radius = 1.43, mass = 209.0, Z = 83.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Po"
        return (radius = 1.35, mass = 209.0, Z = 84.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "At"
        return (radius = 1.27, mass = 210.0, Z = 85.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    elseif name == "Rn"
        return (radius = 1.20, mass = 222.0, Z = 86.0, dte = dte, bde = bde, alpha = 1.0, beta = 0.44)
    else
        error("Unknown element: $name")
    end
end
