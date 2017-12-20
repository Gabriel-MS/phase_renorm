#Thermodynamical correlations
import data
import numpy as np

R = 8.314462175e-6 #m3.MPa/K/mol

#Saturation pressure calculation, Antoine Equation--------------------------------------
def Psat_antoine(IDs,T):
    A = np.array(data.A_antoine(IDs))
    B = np.array(data.B_antoine(IDs))
    C = np.array(data.C_antoine(IDs))
    
    Psat = 10 ** (A-B/(T+C))
    
    #Convert bar to MPa
    Psat = Psat/10
    return Psat
#=======================================================================================

#ideal Cp calculation, SVN correlation--------------------------------------------------
def ideal_cp(IDs,T):
    #Cp/R = A + B + C*T^2 + D*T^(-2)
    coef = np.array(data.cp_ig_coeff(IDs))
    A = coef[0]
    B = coef[1]*1e-3
    C = coef[2]*1e-6
    D = coef[3]*1e-5
    
    Cp = R * (A + B*T + C*(T**2) + D*(T**(-2)))
    
    print A,B,C,D,Cp,T
    
    return Cp
#=======================================================================================
