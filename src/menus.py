#This file is dedicated to user interface menus
import itertools
import data

#Brings out user menu--------------------------------------
def user():
    #Define number of compounds
    nc = input('Number of compounds = ')
    print '\n'
    
    #Define ID of compounds
    print 'Components ID: '
    compID = [int(x) for x in raw_input().split()]
    print '\n'
    
    #Define Feed Composition
    print 'Feed Composition: '
    z = [float(x) for x in raw_input().split()]
    print '\n'
    
    #Define EoS
    print 'EoS options: \n 1.SRK  2.SRK+RG \n 3.PR   4.PR+RG \n 5.CPA  6.CPA+RG \n 7.PC-SAFT  8.PC-SAFT+RG'
    EoS = input('EoS = ')
    print '\n'
    
    #Define MR
    print 'MR options: \n 1.VdW1f'
    MR = input('MR = ')
    print '\n'
    
    if EoS==5:
        #Define AR
        print 'AR options: \n 1. 1A \n 2. 2A \n 3. 2B \n 4. 3A \n 5. 3B \n 6. 4A \n 7. 4B \n 8. 4C'
        AR = [int(x) for x in raw_input().split()]
        print '\n'
        
        #Define CR
        print 'CR options: \n 1. CR-1 \n 2. CR-2 \n 3. CR-3 \n 4. CR-4 \n 5. ECR'
        CR = input('CR = ')
        print '\n'
        
    else:
        AR = 'NULL'
        CR = 'NULL'
        
    user_options = [nc,compID,EoS,MR,z,AR,CR]
    return user_options   
#==========================================================

#Gives String of EoS---------------------------------------
def printeos(eos):
    return{
        1: 'SRK',
        2: 'SRK+RG',
        3: 'PR',
        4: 'PR+RG',
        5: 'CPA',
        6: 'CPA+RG',
        7: 'PC-SAFT',
        8: 'PC-SAFT+RG'
    }.get(eos,'NULL')
#==========================================================

#Gives String of MR---------------------------------------
def printmr(MR):
    return{
        1: 'VdW1f'
    }.get(MR,'NULL')
#==========================================================

#Gives String of AR---------------------------------------
def printar(AR):
    return{
        1: '1A',
        2: '2A',
        3: '2B',
        4: '3A',
        5: '3B',
        6: '4A',
        7: '4B',
        8: '4C'
    }.get(AR,'NULL')
#==========================================================

#Gives String of AR---------------------------------------
def printcr(CR):
    return{
        1: 'CR-1',
        2: 'CR-2',
        3: 'CR-3',
        4: 'CR-4',
        5: 'ECR'
    }.get(CR,'NULL')
#==========================================================

#Gives String of Components--------------------------------
def printcompounds(comp):
    names = list()
    names = data.names()
    name = names[comp-1]
    return name
#==========================================================

#Flatten a list of lists into single list------------------
def flatten(xs):
    result = []
    if isinstance(xs, (list, tuple)):
        for x in xs:
            result.extend(flatten(x))
    else:
        result.append(xs)
    return result
#==========================================================

#Organize list to print user options-----------------------
def print_options(options):
    options_to_print = []
    options_to_print.append(options[0])
    #unpack compounds in options[1]
    names = flatten(options)
    names_str = []
    for i in range (1,options[0]+1):
        names_str.append(printcompounds(names[i]))
    options_to_print.append(names_str)
    options_to_print.append(printeos(options[2]))
    options_to_print.append(printmr(options[3]))
    options_to_print.append(options[4])
    #unpack ARs in options[1]
    ARs = flatten(options[5])
    ARs_str = []
    for i in range(0,options[0]):
        ARs_str.append(printar(ARs[i]))
    options_to_print.append(ARs_str)
    options_to_print.append(printcr(options[6]))
    return options_to_print
#==========================================================