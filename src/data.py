import csv
import menus
import pandas as pd
import numpy as np

#Function to open csv databank and store its values in specific lists---
def properties():
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        names = []
        Tcs = []
        for row in readCSV:
            name = row[1]           #Names of Compounds
            Tc   = float(row[2])    #Critical Temperature
            #Organize read values to respective lists
            names.append(name)
            Tcs.append(Tc)
    #Printing out for apparent no reason (test)
    print(names)
    print(Tcs)
#=======================================================================

#Function to open csv databank and return names of compounds------------
def names():
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        names = []
        for row in readCSV:
            name = row[1]           #Names of Compounds
            #Organize read values to respective lists
            names.append(name)
    return names
#=======================================================================
   
#Function to open csv databank and return critical temperatures---------
def Tc(ID):
    comp = menus.flatten(ID)
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        Tcs = []
        for row in readCSV:
            Tc = float(row[2])           #Critical Temperature
            #Organize read values to respective lists
            Tcs.append(Tc)
    #Finding specific Tcs
    Tcs = menus.flatten(Tcs)
    Tc_ID = []
    #Organizing and fingind specific critical temperatures
    for i in comp:
        Tc_ID.append(Tcs[i-1])
    return Tc_ID
#=======================================================================
    
#Function to open csv databank and return critical pressures------------
def Pc(ID):
    comp = menus.flatten(ID)
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        Pcs = []
        for row in readCSV:
            Pc = float(row[3])           #Critical pressure
            #Organize read values to respective lists
            Pcs.append(Pc)
    #Finding specific Pcs
    Pcs = menus.flatten(Pcs)
    Pc_ID = []
    #Organizing and fingind specific critical pressures
    for i in comp:
        Pc_ID.append(Pcs[i-1])
    return Pc_ID
#=======================================================================
    
#Function to open csv databank and return acentric factors--------------
def omega(ID):
    comp = menus.flatten(ID)
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        omegas = []
        for row in readCSV:
            omega = float(row[4])           #acentric factor
            #Organize read values to respective lists
            omegas.append(omega)
    #Finding specific omegas
    omegas = menus.flatten(omegas)
    omega_ID = []
    #Organizing and fingind specific acentric factors
    for i in comp:
        omega_ID.append(omegas[i-1])
    return omega_ID
#=======================================================================
   
#Function to get component mass----------
def mass(ID):
    
    #Read values
    df = pd.read_csv('../input/Mw.csv',sep=';',header=None)

    out = []
    out.append(df[2][ID[0]-1])

    return out
#=======================================================================   
    
#Function to open csv databank and return A antoine parameters----------
def A_antoine(ID):
    comp = menus.flatten(ID)
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        As = []
        for row in readCSV:
            A = float(row[5])           #A antoine parameter
            #Organize read values to respective lists
            As.append(A)
    #Finding specific As
    As = menus.flatten(As)
    A_ID = []
    #Organizing and fingind specific A antoine parameters
    for i in comp:
        A_ID.append(As[i-1])
    return A_ID
#=======================================================================

#Function to open csv databank and return B antoine parameters----------
def B_antoine(ID):
    comp = menus.flatten(ID)
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        Bs = []
        for row in readCSV:
            B = float(row[6])           #B antoine parameter
            #Organize read values to respective lists
            Bs.append(B)
    #Finding specific Bs
    Bs = menus.flatten(Bs)
    B_ID = []
    #Organizing and fingind specific B antoine parameters
    for i in comp:
        B_ID.append(Bs[i-1])
    return B_ID
#=======================================================================

#Function to open csv databank and return C antoine parameters----------
def C_antoine(ID):
    comp = menus.flatten(ID)
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        Cs = []
        for row in readCSV:
            C = float(row[7])           #C antoine parameter
            #Organize read values to respective lists
            Cs.append(C)
    #Finding specific Cs
    Cs = menus.flatten(Cs)
    C_ID = []
    #Organizing and fingind specific C antoine parameters
    for i in comp:
        C_ID.append(Cs[i-1])
    return C_ID
#=======================================================================

#Function to get Cp ideal coefficients for correlation------------------
def cp_ig_coeff(ID):
    
    #Read values
    df = pd.read_csv('../input/ideal_Cp.csv',sep=';',header=None)

    out = []
    out.append(df[2][ID[0]-1])
    out.append(df[3][ID[0]-1])
    out.append(df[4][ID[0]-1])
    out.append(df[5][ID[0]-1])

    return out
#=======================================================================

#Function to open csv databank and return a0 CPA parameters-------------
def a0(ID):
    comp = menus.flatten(ID)
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        a0s = []
        for row in readCSV:
            a0 = float(row[8])           #a0 parameter
            #Organize read values to respective lists
            a0s.append(a0)
    #Finding specific a0s
    a0s = menus.flatten(a0s)
    a0_ID = []
    #Organizing and fingind specific a0 parameters
    for i in comp:
        a0_ID.append(a0s[i-1])
    return a0_ID
#=======================================================================

#Function to open csv databank and return b CPA parameters--------------
def bCPA(ID):
    comp = menus.flatten(ID)
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        bCPAs = []
        for row in readCSV:
            bCPA = float(row[9])           #bCPA parameter
            #Organize read values to respective lists
            bCPAs.append(bCPA)
    #Finding specific bCPAs
    bCPAs = menus.flatten(bCPAs)
    bCPA_ID = []
    #Organizing and fingind specific bCPA parameters
    for i in comp:
        bCPA_ID.append(bCPAs[i-1])
    return bCPA_ID
#=======================================================================

#Function to open csv databank and return c1 CPA parameters-------------
def c1(ID):
    comp = menus.flatten(ID)
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        c1s = []
        for row in readCSV:
            c1 = float(row[10])           #c1 parameter
            #Organize read values to respective lists
            c1s.append(c1)
    #Finding specific c1s
    c1s = menus.flatten(c1s)
    c1_ID = []
    #Organizing and fingind specific c1 parameters
    for i in comp:
        c1_ID.append(c1s[i-1])
    return c1_ID
#=======================================================================

#Function to open csv databank and return en CPA parameters-------------
def en(ID):
    comp = menus.flatten(ID)
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        ens = []
        for row in readCSV:
            en = float(row[11])           #en parameter
            #Organize read values to respective lists
            ens.append(en)
    #Finding specific ens
    ens = menus.flatten(ens)
    en_ID = []
    #Organizing and fingind specific en parameters
    for i in comp:
        en_ID.append(ens[i-1])
    return en_ID
#=======================================================================

#Function to opbeta csv databank and return beta CPA parameters---------
def beta(ID):
    comp = menus.flatten(ID)
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        betas = []
        for row in readCSV:
            beta = float(row[12])           #beta parameter
            #Organize read values to respective lists
            betas.append(beta)
    #Finding specific betas
    betas = menus.flatten(betas)
    beta_ID = []
    #Organizing and fingind specific beta parameters
    for i in comp:
        beta_ID.append(betas[i-1])
    return beta_ID
#=======================================================================

#Function to opL csv databank and return L renormalization parameters---
def L(ID):
    comp = menus.flatten(ID)
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        Ls = []
        for row in readCSV:
            L = float(row[13])           #L parameter
            #Organize read values to respective lists
            Ls.append(L)
    #Finding specific Ls
    Ls = menus.flatten(Ls)
    L_ID = []
    #Organizing and fingind specific L parameters
    for i in comp:
        L_ID.append(Ls[i-1])
    return L_ID
#=======================================================================

#Function to opphi csv databank and return phi renormalization parameters---
def phi(ID):
    comp = menus.flatten(ID)
    with open('properties.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        phis = []
        for row in readCSV:
            phi = float(row[14])           #phi parameter
            #Organize read values to respective lists
            phis.append(phi)
    #Finding specific phis
    phis = menus.flatten(phis)
    phi_ID = []
    #Organizing and fingind specific phi parameters
    for i in comp:
        phi_ID.append(phis[i-1])
    return phi_ID
#===========================================================================

#Funciton to truncate given value and return with chosen decimal digits-----
def truncate(f, digits):
    return float(("{:.30f}".format(f))[:-30+digits])
#===========================================================================

#Function to open csv file with experimental T-density curve and output data---
def loadexp4():
    with open('../input/Methanol.csv') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        next(readCSV, None)  # skip the header
        dens_vaps = []
        dens_liqs = []
        Ts = []
        flag = False
        for row in readCSV:
            try:
                if flag==False:
                    T = float(row[0])
                    Ts.append(T)
                    dens_liq = float(row[2])
                    dens_liqs.append(dens_liq)
                else:
                    dens_vap = float(row[2])
                    dens_vaps.append(dens_vap)
            except ValueError:
                T = 999.0
                flag=True
    expdata = []
    expdata.append(dens_vaps)
    expdata.append(dens_liqs)
    expdata.append(Ts)
    return expdata
#=============================================================================

#Function to open csv file and extract one column of data---------------------
def loadexp(filename):
    filen = ('../input/%s' %filename)
    with open(filen) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=';')
        expdata = []
        flag = False
        for row in readCSV:
            dt = float(row[0])
            expdata.append(dt)
    return expdata
#=============================================================================

#Function to import dataframe-------------------------------------------------
def loadexp2(filename):

    #Read values
    df = pd.read_csv('../input/%s'%filename,sep=';',header=None)

    out = []
    out.append(df[:][0].values.tolist())
    out.append(df[:][1].values.tolist())
    out.append(df[:][2].values.tolist())
    out.append(df[:][3].values.tolist())
    out.append(df[:][4].values.tolist())

    return out
#=============================================================================

#Function to import dataframe-------------------------------------------------
def loadexp3(filename):

    #Read values
    df = pd.read_csv('../input/%s'%filename,sep=';',header=None)

    out = []
    out.append(df[:][0].values.tolist())
    out.append(df[:][1].values.tolist())
    out.append(df[:][2].values.tolist())
    out.append(df[:][3].values.tolist())

    return out
#=============================================================================

#Function to modify custom substance association parameters-------------------
def modify_assoc(ID,e,b):

    lin = ID[0]-1

    #Read values
    df = pd.read_csv('properties.csv',sep=';',header=None)

    #Modify association energy and association volume
    df.set_value(lin,11,e)
    df.set_value(lin,12,b)

    #Save new values
    df.to_csv('properties.csv',sep=';',index=False,index_label=False,header=False)

    return df
#=============================================================================

#Function to modify custom substance association parameters-------------------
def modify_CPA(ID,vec):

    lin = ID[0]-1

    #Read values
    df = pd.read_csv('properties.csv',sep=';',header=None)

    #Modify association energy and association volume
    a0 = vec[0]
    bCPA = vec[1]
    c1 = vec[2]
    e = vec[3]
    b = vec[4]

    df.set_value(lin,8,a0)
    df.set_value(lin,9,bCPA)
    df.set_value(lin,10,c1)
    df.set_value(lin,11,e)
    df.set_value(lin,12,b)

    #Save new values
    df.to_csv('properties.csv',sep=';',index=False,index_label=False,header=False)

    return df
#=============================================================================

#Function to modify custom substance association parameters-------------------
def modify_renorm(ID,vec):

    lin = ID[0]-1

    #Read values
    df = pd.read_csv('properties.csv',sep=';',header=None)

    #Modify association energy and association volume
    L = vec[0]
    phi = vec[1]

    df.set_value(lin,13,L)
    df.set_value(lin,14,phi)

    #Save new values
    df.to_csv('properties.csv',sep=';',index=False,index_label=False,header=False)

    return df
#=============================================================================

#Function to modify calculated isotherm---------------------------------------
def modify_isotherm(P,rho):

    #Read values
    df = pd.read_csv('../output/isotherm.csv',sep=';',header=None)

    #Change values
    for i in range(0,len(rho)):
        df.set_value(i,0,rho[i])
        df.set_value(i,1,P[i])

    #Save new values
    df.to_csv('../output/isotherm.csv',sep=';',index=False,index_label=False,header=False)
#=============================================================================

#Function to modify calculated isotherm---------------------------------------
def modify_isotherm2(rho,P,dP,d2P):

    #Read values
    df = pd.read_csv('../output/isotherm.csv',sep=';',header=None)

    #Change values
    for i in range(0,len(rho)):
        df.set_value(i,0,rho[i])
        df.set_value(i,1,P[i])
        df.set_value(i,2,dP[i])
        df.set_value(i,3,d2P[i])

    #Save new values
    df.to_csv('../output/isotherm.csv',sep=';',index=False,index_label=False,header=False)
#===============================================================================

#Function to output values to calculate surface for objetctive functions--------
def surface_out(L,phi,Fobj_dens_liq,Fobj_dens_vap,Fobj_P,Tc,Pc,rhoc):

    #Read values
    df = pd.read_csv('../output/surf.csv',sep=';',header='infer')
    colnames = []
    colnames = pd.read_csv('../output/surf.csv',sep=';',header='infer').columns.tolist()
    
    data2 = np.array([['',colnames[0],colnames[1],colnames[2],colnames[3],colnames[4],colnames[5],colnames[6],colnames[7]],
                ['0',L,phi,Fobj_dens_liq,Fobj_dens_vap,Fobj_P,Tc,Pc,rhoc]])
                
    df2 = pd.DataFrame(data=data2[1:,1:],
                  columns=data2[0,1:])
                  
    df3 = df.append(df2)

    #Save new values
    df3.to_csv('../output/surf.csv',sep=';',header=colnames,index=False)
#=============================================================================
