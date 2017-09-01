#Thermodynamical correlations
import data
import numpy as np

#Saturation pressure calculation, Antoine Equation--------------------------------------
def Psat_antoine(IDs,T):
    A = np.array(data.A_antoine(IDs))
    B = np.array(data.B_antoine(IDs))
    C = np.array(data.C_antoine(IDs))
    
    Psat = 10 ** (A-B/(T+C))
    
    #Convert bar to MPa
    #Psat = Psat/10
    return Psat
#=======================================================================================