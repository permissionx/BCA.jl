struct Element
    name::String
    radius::Float64
    mass::Float64
    Z::Float64
    dte::Float64
    bde::Float64
    alpha::Float64
    beta::Float64
end


function Element(name::String, dte::Float64, bde::Float64)
    if name == "H"
        return Element(name, 0.53, 1.0, 1.0, dte, bde, 1.0, 0.05)
    elseif name == "He"
        return Element(name, 0.31, 4.0, 2.0, dte, bde, 1.0, 0.05)
    elseif name == "Li"
        return Element(name, 1.67, 7.0, 3.0, dte, bde, 1.5, 0.05)
    elseif name == "Be"
        return Element(name, 1.12, 9.0, 4.0, dte, bde, 1.5, 0.44)
    elseif name == "B"
        return Element(name, 0.87, 11.0, 5.0, dte, bde, 1.5, 0.44)
    elseif name == "C"
        return Element(name, 0.67, 12.0, 6.0, dte, bde, 1.0, 0.44)
    elseif name == "N"
        return Element(name, 0.56, 14.0, 7.0, dte, bde, 1.0, 0.44)
    elseif name == "O"
        return Element(name, 0.48, 16.0, 8.0, dte, bde, 1.0, 0.44)
    elseif name == "F"
        return Element(name, 0.42, 19.0, 9.0, dte, bde, 1.1, 0.15)
    elseif name == "Ne"
        return Element(name, 0.38, 20.0, 10.0, dte, bde, 1.0, 0.44)
    elseif name == "Na"
        return Element(name, 1.90, 23.0, 11.0, dte, bde, 1.0, 0.44)
    elseif name == "Mg"
        return Element(name, 1.45, 24.0, 12.0, dte, bde, 1.0, 0.44)
    elseif name == "Al"
        return Element(name, 1.18, 27.0, 13.0, dte, bde, 1.0, 0.12)
    elseif name == "Si"
        return Element(name, 1.11, 28.0, 14.0, dte, bde, 1.0, 0.44)
    elseif name == "P"
        return Element(name, 0.98, 31.0, 15.0, dte, bde, 1.25, 0.44)
    elseif name == "S"
        return Element(name, 0.88, 32.0, 16.0, dte, bde, 1.0, 0.44)
    elseif name == "Cl"
        return Element(name, 0.79, 35.0, 17.0, dte, bde, 1.0, 0.44)
    elseif name == "Ar"
        return Element(name, 0.71, 40.0, 18.0, dte, bde, 1.0, 0.44)
    elseif name == "K"
        return Element(name, 2.43, 39.0, 19.0, dte, bde, 1.0, 0.44)
    elseif name == "Ca"
        return Element(name, 1.94, 40.0, 20.0, dte, bde, 1.0, 0.44)
    elseif name == "Sc"
        return Element(name, 1.84, 45.0, 21.0, dte, bde, 1.0, 0.44)
    elseif name == "Ti"
        return Element(name, 1.76, 48.0, 22.0, dte, bde, 1.0, 0.44)
    elseif name == "V"
        return Element(name, 1.71, 51.0, 23.0, dte, bde, 1.0, 0.44)
    elseif name == "Cr"
        return Element(name, 1.66, 52.0, 24.0, dte, bde, 1.0, 0.44)
    elseif name == "Mn"
        return Element(name, 1.61, 55.0, 25.0, dte, bde, 1.0, 0.44)
    elseif name == "Fe"
        return Element(name, 1.56, 56.0, 26.0, dte, bde, 1.0, 0.44)
    elseif name == "Co"
        return Element(name, 1.52, 59.0, 27.0, dte, bde, 1.0, 0.44)
    elseif name == "Ni"
        return Element(name, 1.49, 58.0, 28.0, dte, bde, 1.0, 0.44)
    elseif name == "Cu"
        return Element(name, 1.45, 63.0, 29.0, dte, bde, 1.0, 0.44)
    elseif name == "Zn"
        return Element(name, 1.42, 64.0, 30.0, dte, bde, 1.0, 0.44)
    elseif name == "Ga"
        return Element(name, 1.36, 69.0, 31.0, dte, bde, 1.0, 0.44)
    elseif name == "Ge"
        return Element(name, 1.25, 74.0, 32.0, dte, bde, 1.0, 0.44)
    elseif name == "As"
        return Element(name, 1.14, 75.0, 33.0, dte, bde, 1.0, 0.44)
    elseif name == "Se"
        return Element(name, 1.03, 80.0, 34.0, dte, bde, 1.0, 0.44)
    elseif name == "Br"
        return Element(name, 0.94, 79.0, 35.0, dte, bde, 1.0, 0.44)
    elseif name == "Kr"
        return Element(name, 0.88, 84.0, 36.0, dte, bde, 1.0, 0.44)
    elseif name == "Rb"
        return Element(name, 2.65, 85.0, 37.0, dte, bde, 1.0, 0.44)
    elseif name == "Sr"
        return Element(name, 2.19, 88.0, 38.0, dte, bde, 1.0, 0.44)
    elseif name == "Y"
        return Element(name, 2.12, 89.0, 39.0, dte, bde, 1.0, 0.44)
    elseif name == "Zr"
        return Element(name, 2.06, 90.0, 40.0, dte, bde, 1.0, 0.44)
    elseif name == "Nb"
        return Element(name, 1.98, 93.0, 41.0, dte, bde, 1.0, 0.44)
    elseif name == "Mo"
        return Element(name, 1.90, 98.0, 42.0, dte, bde, 1.0, 0.44)
    elseif name == "Tc"
        return Element(name, 1.83, 98.0, 43.0, dte, bde, 1.0, 0.44)
    elseif name == "Ru"
        return Element(name, 1.78, 102.0, 44.0, dte, bde, 1.0, 0.44)
    elseif name == "Rh"
        return Element(name, 1.73, 103.0, 45.0, dte, bde, 1.0, 0.44)
    elseif name == "Pd"
        return Element(name, 1.69, 106.0, 46.0, dte, bde, 1.0, 0.44)
    elseif name == "Ag"
        return Element(name, 1.65, 107.0, 47.0, dte, bde, 1.0, 0.44)
    elseif name == "Cd"
        return Element(name, 1.61, 114.0, 48.0, dte, bde, 1.0, 0.44)
    elseif name == "In"
        return Element(name, 1.56, 115.0, 49.0, dte, bde, 1.0, 0.44)
    elseif name == "Sn"
        return Element(name, 1.45, 120.0, 50.0, dte, bde, 1.0, 0.44)
    elseif name == "Sb"
        return Element(name, 1.33, 121.0, 51.0, dte, bde, 1.0, 0.44)
    elseif name == "Te"
        return Element(name, 1.23, 130.0, 52.0, dte, bde, 1.0, 0.44)
    elseif name == "I"
        return Element(name, 1.15, 129.0, 53.0, dte, bde, 1.0, 0.44)
    elseif name == "Xe"
        return Element(name, 1.08, 132.0, 54.0, dte, bde, 1.0, 0.44)
    elseif name == "Cs"
        return Element(name, 2.98, 133.0, 55.0, dte, bde, 1.0, 0.44)
    elseif name == "Ba"
        return Element(name, 2.53, 138.0, 56.0, dte, bde, 1.0, 0.44)
    elseif name == "La"
        return Element(name, 2.26, 139.0, 57.0, dte, bde, 1.0, 0.44)
    elseif name == "Ce"
        return Element(name, 2.10, 140.0, 58.0, dte, bde, 1.0, 0.44)
    elseif name == "Pr"
        return Element(name, 2.47, 141.0, 59.0, dte, bde, 1.0, 0.44)
    elseif name == "Nd"
        return Element(name, 2.06, 142.0, 60.0, dte, bde, 1.0, 0.44)
    elseif name == "Pm"
        return Element(name, 2.05, 145.0, 61.0, dte, bde, 1.0, 0.44)
    elseif name == "Sm"
        return Element(name, 2.38, 152.0, 62.0, dte, bde, 1.0, 0.44)
    elseif name == "Eu"
        return Element(name, 2.31, 153.0, 63.0, dte, bde, 1.0, 0.44)
    elseif name == "Gd"
        return Element(name, 2.33, 158.0, 64.0, dte, bde, 1.0, 0.44)
    elseif name == "Tb"
        return Element(name, 2.25, 159.0, 65.0, dte, bde, 1.0, 0.44)
    elseif name == "Dy"
        return Element(name, 2.28, 164.0, 66.0, dte, bde, 1.0, 0.44)
    elseif name == "Ho"
        return Element(name, 2.26, 165.0, 67.0, dte, bde, 1.0, 0.44)
    elseif name == "Er"
        return Element(name, 2.26, 170.0, 68.0, dte, bde, 1.0, 0.44)
    elseif name == "Tm"
        return Element(name, 2.22, 169.0, 69.0, dte, bde, 1.0, 0.44)
    elseif name == "Yb"
        return Element(name, 2.22, 174.0, 70.0, dte, bde, 1.0, 0.44)
    elseif name == "Lu"
        return Element(name, 2.17, 176.0, 71.0, dte, bde, 1.0, 0.44)
    elseif name == "Hf"
        return Element(name, 2.08, 180.0, 72.0, dte, bde, 1.0, 0.44)
    elseif name == "Ta"
        return Element(name, 2.00, 181.0, 73.0, dte, bde, 1.0, 0.44)
    elseif name == "W"
        return Element(name, 1.93, 184.0, 74.0, dte, bde, 1.0, 0.44)
    elseif name == "Re"
        return Element(name, 1.88, 187.0, 75.0, dte, bde, 1.0, 0.44)
    elseif name == "Os"
        return Element(name, 1.85, 192.0, 76.0, dte, bde, 1.0, 0.44)
    elseif name == "Ir"
        return Element(name, 1.80, 191.0, 77.0, dte, bde, 1.0, 0.44)
    elseif name == "Pt"
        return Element(name, 1.77, 195.0, 78.0, dte, bde, 1.0, 0.44)
    elseif name == "Au"
        return Element(name, 1.74, 197.0, 79.0, dte, bde, 1.0, 0.44)
    elseif name == "Hg"
        return Element(name, 1.71, 202.0, 80.0, dte, bde, 1.0, 0.44)
    elseif name == "Tl"
        return Element(name, 1.56, 205.0, 81.0, dte, bde, 1.0, 0.44)
    elseif name == "Pb"
        return Element(name, 1.54, 208.0, 82.0, dte, bde, 1.0, 0.44)
    elseif name == "Bi"
        return Element(name, 1.43, 209.0, 83.0, dte, bde, 1.0, 0.44)
    elseif name == "Po"
        return Element(name, 1.35, 209.0, 84.0, dte, bde, 1.0, 0.44)
    elseif name == "At"
        return Element(name, 1.27, 210.0, 85.0, dte, bde, 1.0, 0.44)
    elseif name == "Rn"
        return Element(name, 1.20, 222.0, 86.0, dte, bde, 1.0, 0.44)
    else
        error("Unknown element: $name")
    end
end
