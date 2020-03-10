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

"""
user_options = []
user_options.append(8) #nc
user_options.append([21,22,23,24,25,26,27,28]) #IDs
user_options.append(3) #EoS
user_options.append(1) #MR
user_options.append([0.3985,0.0512,0.0286,0.0042,0.0079,0.0018,0.0025,0.5016]) #z
user_options.append([3,3,3,3,3,3,3,3]) #CPA AM
user_options.append(1) #CPA CR
print_options = []
print_options = menus.print_options(user_options)
"""

user_options = []
user_options.append(2) #nc
user_options.append([31,29])#29,30,31]) #IDs #C20,CH4,CO2
user_options.append(1) #EoS
user_options.append(1) #MR
user_options.append([0.9953,0.0047])#0.2,0.3,0.5]) #z
user_options.append([8,8]) #CPA AM
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
env_type = 1
P = 0.05
T = 200.00
kij = np.zeros((nc,nc))
#kij = np.array([[0.0, 0.0, 0.02],
#                [0.0, 0.0, 0.0292],
#                [0.02, 0.0292, 0.0]])

#CPA auto-association configurations
auto = []
auto = association.CPA_auto(AR,nc,IDs)
en_auto = auto[0]
beta_auto = auto[1]
#CPA cross-association configurations
SM = np.array([[1.0, 0.0],
               [0.0, 1.0]])
R = 8.314462175e-6 #m3.MPa/K/mol

a=['{:.4f}'.format(x) for x in print_options[4]]
print print_options[1]+[print_options[2]]+a


#CO2+C20
#SRK
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#PR
user_options[2] = 3 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#CPA
user_options[2] = 5 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

#CH4+CO2+nC20
user_options[0] = 3 #nc
user_options[5] = [8,8,8] #CPA AM
user_options[1] = [30,31,29]#29,30,31]) #IDs #C20,CH4,CO2
user_options[2] = 1 #EoS
nc = user_options[0] #Defining number of components
IDs = user_options[1] #Defining IDs of components
EoS = user_options[2] #Define Equation of State
MR = user_options[3] #Define Mixing Rule
z = user_options[4] #Define Feed Composition
AR = user_options[5] #Define association rules - CPA
CR = user_options[6] #Define combining rule - CPA

#CPA auto-association configurations
auto = []
auto = association.CPA_auto(AR,nc,IDs)
en_auto = auto[0]
beta_auto = auto[1]
SM = np.array([[1.0, 0.0, 0.0],
               [0.0, 1.0, 0.0],
               [0.0, 0.0, 1.0]])

#20% A----------------------------------------------
print ("Begin 20%A+80%C")
user_options[4] = [0.2/3,0.2/3*2,0.8]
#SRK
user_options[2] = 1 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#PR
user_options[2] = 3 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#CPA
user_options[2] = 5 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

user_options[4] = [0.2/2,0.2/2,0.8]
#SRK
user_options[2] = 1 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#PR
user_options[2] = 3 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#CPA
user_options[2] = 5 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

user_options[4] = [0.2/3*2,0.2/3,0.8]
#SRK
user_options[2] = 1 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#PR
user_options[2] = 3 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#CPA
user_options[2] = 5 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

#40% A----------------------------------------------
print ("Begin 40%A+60%C")
user_options[4] = [0.4/3,0.4/3*2,0.6]
#SRK
user_options[2] = 1 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#PR
user_options[2] = 3 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#CPA
user_options[2] = 5 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

user_options[4] = [0.4/2,0.4/2,0.6]
#SRK
user_options[2] = 1 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#PR
user_options[2] = 3 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#CPA
user_options[2] = 5 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

user_options[4] = [0.4/3*2,0.4/3,0.6]
#SRK
user_options[2] = 1 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#PR
user_options[2] = 3 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#CPA
user_options[2] = 5 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

#60% A----------------------------------------------
print ("Begin 60%A+40%C")
user_options[4] = [0.6/3,0.6/3*2,0.4]
#SRK
user_options[2] = 1 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#PR
user_options[2] = 3 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#CPA
user_options[2] = 5 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

user_options[4] = [0.6/2,0.6/2,0.4]
#SRK
user_options[2] = 1 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#PR
user_options[2] = 3 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#CPA
user_options[2] = 5 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

user_options[4] = [0.6/3*2,0.6/3,0.4]
#SRK
user_options[2] = 1 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#PR
user_options[2] = 3 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#CPA
user_options[2] = 5 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

#80% A----------------------------------------------
print ("Begin 80%A+20%C")
user_options[4] = [0.8/3,0.8/3*2,0.2]
#SRK
user_options[2] = 1 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#PR
user_options[2] = 3 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#CPA
user_options[2] = 5 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

user_options[4] = [0.8/2,0.8/2,0.2]
#SRK
user_options[2] = 1 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#PR
user_options[2] = 3 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#CPA
user_options[2] = 5 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

user_options[4] = [0.8/3*2,0.8/3,0.2]
#SRK
user_options[2] = 1 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#PR
user_options[2] = 3 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)
#CPA
user_options[2] = 5 #EoS
print_options = []
print_options = menus.print_options(user_options)
EoS = user_options[2] #Define Equation of State
z = user_options[4] #Define Feed Composition
envelope.calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type)

#End calculations==========================================================================================

print '\nEnd Program'
