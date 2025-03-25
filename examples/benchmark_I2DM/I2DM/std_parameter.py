#Constants
CONSTANT_E       = 1.6021892e-19                   #C,the unit of  charge
CONSTANT_MU      = 1.661e-27                       #kg,mass of u
CONSTANT_A       = 1e-10                           #m, angstrom
CONSTANT_EV      = 1.6021892e-19                   #J, the energy of micro world
CONSTANT_K       = 8.987551788e9                   #N*m*2/C^2, the Coulumb constant
CONSTANT_AB      = 5.2917721067e-11                #m,Born radius
CONSTANT_VB      = 2.19e6                          #m/s,Born speed
#CONSTANT_ME      = 9.10956e-31                     #kg,mass of e
#CONSTANT_C       = 299792458.0                     #m/s, the velocity of light
#CONSTANT_EPSILON = 8.854187817e-12                 #F/m，the vacuum dielectric constant
#CONSTANT_ER      = 13.6056                         #eV,Reyberg energy，eV
#CONSTANT_H       = 6.62607015e-34                  #J*s,Planck constant
#Element parameters
class Element:
    def __init__(self,name,charge,mass,radius,alpha,beta) -> None:
        self.name=name
        self.charge=charge      #unit：e
        self.mass=mass          #unit：Mu
        self.radius = radius    #unit：A
        self.alpha = alpha      #factor for Se
        self.beta = beta        #factor for x_nl
ELE_H  = Element('H',  1,  1,  0.53,1,0.05)
ELE_HE = Element('He', 2,  4,  0.31,1,0.05)
ELE_LI = Element('Li', 3,  7,  1.67,1.5,0.05)
ELE_BE = Element('Be', 4,  9,  1.12,1.5,0.44)
ELE_B  = Element('B',  5,  11, 0.87,1.5,0.44)
ELE_C  = Element('C',  6,  12, 0.67,1.0,0.44)
ELE_N  = Element('N',  7,  14, 0.56,1.0,0.44)
ELE_O  = Element('O',  8,  16, 0.48,1.0,0.44)
ELE_F  = Element('F',  9,  19, 0.42,1.1,0.15)
ELE_NE = Element('Ne', 10, 20, 0.38,1.0,0.44)
ELE_NA = Element('Na', 11, 23, 1.90,1.0,0.44)
ELE_MG = Element('Mg', 12, 24, 1.45,1.0,0.44)
ELE_AL = Element('Al', 13, 27, 1.18,1.0,0.12)
ELE_SI = Element('Si', 14, 28, 1.11,1.0,0.44)
ELE_P  = Element('P',  15, 31, 0.98,1.25,0.44)
ELE_S  = Element('S',  16, 32, 0.88,1.0,0.44)
ELE_CL = Element('Cl', 17, 35, 0.79,1.0,0.44)
ELE_AR = Element('Ar', 18, 40, 0.71,1.0,0.44)
ELE_K  = Element('K',  19, 39, 2.43,1.0,0.44)
ELE_CA = Element('Ca', 20, 40, 1.94,1.0,0.44)
ELE_SC = Element('Sc', 21, 45, 1.84,1.0,0.44)
ELE_TI = Element('Ti', 22, 48, 1.76,1.0,0.44)
ELE_V  = Element('V',  23, 51, 1.71,1.0,0.44)
ELE_CR = Element('Cr', 24, 52, 1.66,1.0,0.44)
ELE_MN = Element('Mn', 25, 55, 1.61,1.0,0.44)
ELE_FE = Element('Fe', 26, 56, 1.56,1.0,0.44)
ELE_CO = Element('Co', 27, 59, 1.52,1.0,0.44)
ELE_NI = Element('Ni', 28, 58, 1.49,1.0,0.44)
ELE_CU = Element('Cu', 29, 63, 1.45,1.0,0.44)
ELE_ZN = Element('Zn', 30, 64, 1.42,1.0,0.44)
ELE_GA = Element('Ga', 31, 69, 1.36,1.0,0.44)
ELE_GE = Element('Ge', 32, 74, 1.25,1.0,0.44)
ELE_AS = Element('As', 33, 75, 1.14,1.0,0.44)
ELE_SE = Element('Se', 34, 80, 1.03,1.0,0.44)
ELE_BR = Element('Br', 35, 79, 0.94,1.0,0.44)
ELE_KR = Element('Kr', 36, 84, 0.88,1.0,0.44)
ELE_RB = Element('Rb', 37, 85, 2.65,1.0,0.44)
ELE_SR = Element('Sr', 38, 88, 2.19,1.0,0.44)
ELE_Y  = Element('Y',  39, 89, 2.12,1.0,0.44)
ELE_ZR = Element('Zr', 40 ,90, 2.06,1.0,0.44)
ELE_NB = Element('Nb', 41, 93, 1.98,1.0,0.44)
ELE_MO = Element('Mo', 42, 98, 1.90,1.0,0.44)
ELE_TC = Element('Tc', 43, 98, 1.83,1.0,0.44)
ELE_RU = Element('Ru', 44, 102,1.78,1.0,0.44)
ELE_RH = Element('Rh', 45, 103,1.73,1.0,0.44)
ELE_PB = Element('Pb', 46, 106,1.69,1.0,0.44)
ELE_AG = Element('Ag', 47, 107,1.65,1.0,0.44)
ELE_CD = Element('Cd', 48, 114,1.61,1.0,0.44)
ELE_IN = Element('In', 49, 115,1.56,1.0,0.44)
ELE_SN = Element('Sn', 50, 120,1.45,1.0,0.44)
ELE_SB = Element('Sb', 51, 121,1.33,1.0,0.44)
ELE_TE = Element('Te', 52, 130,1.23,1.0,0.44)
ELE_I  = Element('I',  53, 129,1.15,1.0,0.44)
ELE_XE = Element('Xe', 54, 132,1.08,1.0,0.44)
ELE_CS = Element('Cs', 55, 133,2.98,1.0,0.44)
ELE_BA = Element('Ba', 56, 138,2.53,1.0,0.44)
ELE_LA = Element('La', 57, 139,2.26,1.0,0.44)
ELE_CE = Element('Ce', 58, 140,2.10,1.0,0.44)
ELE_PR = Element('Pr', 59, 141,2.47,1.0,0.44)
ELE_ND = Element('Nd', 60, 142,2.06,1.0,0.44)
ELE_PM = Element('Pm', 61, 145,2.05,1.0,0.44)
ELE_SM = Element('Sm', 62, 152,2.38,1.0,0.44)
ELE_EU = Element('Eu', 63, 153,2.31,1.0,0.44)
ELE_GD = Element('Gd', 64, 158,2.33,1.0,0.44)
ELE_TB = Element('Tb', 65, 159,2.25,1.0,0.44)
ELE_DY = Element('Dy', 66, 164,2.28,1.0,0.44)
ELE_HO = Element('Ho', 67, 165,2.26,1.0,0.44)
ELE_ER = Element('Er', 68, 170,2.26,1.0,0.44)
ELE_TM = Element('Tm', 69, 169,2.22,1.0,0.44)
ELE_YB = Element('Yb', 70, 174,2.22,1.0,0.44)
ELE_LU = Element('Lu', 71, 176,2.17,1.0,0.44)
ELE_HF = Element('Hf', 72, 180,2.08,1.0,0.44)
ELE_TA = Element('Ta', 73, 181,2.00,1.0,0.44)
ELE_W  = Element('W',  74, 184,1.93,1.0,0.44)
ELE_RE = Element('Re', 75, 187,1.88,1.0,0.44)
ELE_OS = Element('Os', 76, 192,1.85,1.0,0.44)
ELE_IR = Element('Ir', 77, 191,1.80,1.0,0.44)
ELE_PT = Element('Pt', 78, 195,1.77,1.0,0.44)
ELE_AU = Element('Au', 79, 197,1.74,1.0,0.44)
ELE_HG = Element('Hg', 80, 202,1.71,1.0,0.44)
ELE_TL = Element('Tl', 81, 205,1.56,1.0,0.44)
ELE_PB = Element('Pb', 82, 208,1.54,1.0,0.44)
ELE_BI = Element('Bi', 83, 209,1.43,1.0,0.44)
ELE_PO = Element('Po', 84, 209,1.35,1.0,0.44)
ELE_AT = Element('At', 85, 210,1.27,1.0,0.44)
ELE_RN = Element('Rn', 86, 222,1.20,1.0,0.44)

ALL_ELEMENTS = [ELE_H,ELE_HE,ELE_LI,ELE_BE,ELE_B,ELE_C,ELE_N,ELE_O,ELE_F,
                ELE_NE,ELE_NA,ELE_MG,ELE_AL,ELE_SI,ELE_P,ELE_S,ELE_CL,ELE_AR,
                ELE_K,ELE_CA,ELE_SC,ELE_TI,ELE_V,ELE_CR,ELE_MN,ELE_FE,ELE_CO,
                ELE_NI,ELE_CU,ELE_ZN,ELE_GA,ELE_GE,ELE_AS,ELE_SE,ELE_BR,ELE_KR,
                ELE_RB,ELE_SR,ELE_Y,ELE_ZR,ELE_NB,ELE_MO,ELE_TC,ELE_RU,ELE_RH,
                ELE_PB,ELE_AG,ELE_CD,ELE_IN,ELE_SN,ELE_SB,ELE_TE,ELE_I,ELE_XE,
                ELE_CS,ELE_BA,ELE_LA,ELE_CE,ELE_PR,ELE_ND,ELE_PM,ELE_SM,ELE_EU,
                ELE_GD,ELE_TB,ELE_DY,ELE_HO,ELE_ER,ELE_TM,ELE_YB,ELE_LU,ELE_HF,
                ELE_TA,ELE_W,ELE_RE,ELE_OS,ELE_IR,ELE_PT,ELE_AU,ELE_HG,ELE_TL,
                ELE_PB,ELE_BI,ELE_PO,ELE_AT,ELE_RN]
ALL_ELEMENTS_NAME   = tuple(atom.name for atom in ALL_ELEMENTS)
ALL_ELEMENTS_CHARGE = tuple(atom.charge for atom in ALL_ELEMENTS)
ALL_ELEMENTS_MASS   = tuple(atom.mass for atom in ALL_ELEMENTS)
ALL_ELEMENTS_RADIUS = tuple(atom.radius for atom in ALL_ELEMENTS)
ALL_ELEMENTS_ALPHA  = tuple(atom.alpha for atom in ALL_ELEMENTS)
ALL_ELEMENTS_BETA   = tuple(atom.beta for atom in ALL_ELEMENTS)
