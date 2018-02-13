import menus
import data
import eos
import envelope
import association
import correlations
import renormalization
import numpy as np

print ('    ============================================')
print ('    --------Autor: Gabriel Moraes Silva---------')
print ('    ------------------V: 1.0--------------------')
print ('    ----------------Ano: 2017-------------------')
print ('    ---------------ATOMS/EQ/UFRJ----------------')
print ('    ============================================\n\n\n')

#Start user definitions------------------------------------------------------------------------------------
#User menu to define custom configuration
#user_options = menus.user()
user_options = []
user_options.append(1) #nc
user_options.append([1,1]) #IDs
user_options.append(3) #EoS
user_options.append(1) #MR
user_options.append([0.3,0.7]) #z
user_options.append([3,3]) #CPA AM
user_options.append(1) #CPA CR
print_options = []
print_options = menus.print_options(user_options)

#Show user defined menu
print ('    Defined Configuration:----------------------')
print ('    Number of components: %i ' %print_options[0])
print ('    Equation of State:    %s ' %print_options[2])
print ('    Mixing Rule:          %s ' %print_options[3])
print ('    Components:           %s ' %', '.join(print_options[1]))
print ('    Feed Composition:     %s ' %print_options[4])
print ('    Association Rule:     %s ' %', '.join(print_options[5]))
print ('    Combining Rule:       %s ' %print_options[6])
print ('    ============================================\n')
#End user definitions======================================================================================

#Start calculations----------------------------------------------------------------------------------------
#Define Variables
nc = user_options[0] #Defining number of components
IDs = user_options[1] #Defining IDs of components
EoS = user_options[2] #Define Equation of State
MR = user_options[3] #Define Mixing Rule
z = user_options[4] #Define Feed Composition
AR = user_options[5] #Define association rules - CPA
CR = user_options[6] #Define combining rule - CPA
env_type = 5
P = 0.05
T = 300.0
kij = np.zeros((nc,nc))

#CPA auto-association configurations
auto = []
auto = association.CPA_auto(AR,nc,IDs)
en_auto = auto[0]
beta_auto = auto[1]
#CPA cross-association configurations
SM = np.array([[1, 1],
               [1, 1]])
R = 8.314462175e-6 #m3.MPa/K/mol

envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

#End calculations==========================================================================================

print '\nEnd Program'
