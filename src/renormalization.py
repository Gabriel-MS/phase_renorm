import numpy as np
import menus
import data
import eos
import association
import math
import numerical

R = 8.314462175e-6 #m3.MPa/K/mol
NA = 6.023e23      #Avogadro Number
kB = 1.3806503e-23 #Boltzmann constant, J/K

#Outputs renormalized csv file of given EoS at give T-------------------------------------
def renorm(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n):
    #nd    Size of density grid
    #nx    Size of mole fraction grid
    #n     Main loop iteration controller
    
    #If only 1 component is present mimic a binary mixture made of the same component
    if nc==1:
        IDs[1] = IDs[0]
    
    #Recover parameters
    L_rg = data.L(IDs)         #Vector with L parameters (cutoff length)
    phi_rg = data.phi(IDs)     #Vector with phi parameters
    Tc = data.Tc(IDs)
    
    #Components parameters
    a = eos.a_calc(IDs,EoS,T)
    b = eos.b_calc(IDs,EoS)
    Tr = T/np.array(Tc)
    
    #Main loop parameters
    x = np.array([0.001,0.999])
    stepx = float(1/float(nx)) #Step to calculate change
    k = 0               #Vector fill counter
    i = 1               #Main loop counter
    r = 0               #Report counter
    rho = np.empty((nd))            #Density vector
    rhov = []                       #Density vector to export
    x0v = []                        #Mole fraction vector to export
    f = np.empty((nd))              #Helmholtz energy density vector
    fv = []                         #Helmholtz energy density vector to export
    fresv = []                      #Residual Helmholtz energy density vector to export
    Tv = []                         #Temperature values to export
    df = np.empty((nd))             #Changes in helmholtz energy density vector
    f_orig = np.empty((nd))         #Unmodified Helmholtz energy density vector
    X = np.ones((4*nc))
    
    if nc==1:
        X = np.ones((8))
    
    #Main loop*************************************
    while x[0]<1.0:
        if x[0]==0.006: #after first step
            x[0]=0.005
            x[1]=1-x[0]
        
        if nc==1:
            x[0] = 0.99999
            x[1] = 0.00001
        
        #Mixture parameters
        bmix = eos.bmix_calc(MR,b,x)
        amix = eos.amix_calc(MR,a,x,kij)
        rhomax = 0.999999/bmix
        
        #Mixture Renormalization parameters
        L = np.dot(x,np.power(L_rg,3.0))
        L = np.power(L,1.0/3.0)
        phi = np.dot(x,phi_rg)
        
        while k<nd:
            rho[k] = np.array(float(k)/nd/bmix)
            if k==0:
                rho[0] = 1e-5
            if EoS==6:
                Xf = association.frac_nbs(nc,1/rho[k],CR,en_auto,beta_auto,b,bmix,X,0,x,0,T,SM)
                X = Xf
            f[k] = np.array(helm_rep(EoS,R,T,rho[k],amix,bmix,X,x,nc))   #Helmholtz energy density
            f_orig[k] = f[k]                                #Initial helmholtz energy density
            #Subtract attractive forces (due long range correlations)
            f[k] = f[k] + amix*(rho[k]**2)
            k = k+1 
            
        #Main loop****************************************************************
        i = 1
        while i<n:
            K = kB*T/((2**(3*i))*(L**3))
            
            #Long and Short Range forces
            fl = helm_long(EoS,rho,f)
            fs = helm_short(EoS,rho,f,phi,i)
            #Calculate df
            width = rhomax/nd
            w = 0
            while w<nd:
                df[w] = renorm_df(w,nd,fl,fs,K,rho,width)
                w = w+1
                
            #Update Helmholtz Energy Density
            f = f + df
            i = i+1
            
        #Add original attractive forces
        f = f - amix*np.power(rho,2)
        
        #Store residual value of f, calculate helmholtz energy density adding ideal gas energy
        fres = f
        f = f + rho*R*T*(np.log(rho)-1)
        
        fv.append(f)
        fresv.append(fres)
        x0v.append(x[0])
        
        if r==0:
            rhov.append(rho) #rho vector is always the same
        r=1
        
        x[0] = x[0]+stepx
        
    renorm_out = []
    renorm_out.append(fv)
    renorm_out.append(x0v)
    renorm_out.append(rhov)
    renorm_out.append(fresv)
    return renorm_out
#=========================================================================================

#Calculates Helmholtz repulsive forces----------------------------------------------------
def helm_rep(EoS,R,T,rho,amix,bmix,X,x,nc):
    
    f_CPA = 1
    if EoS==6:
        if nc==1:
            nc = 2
        one4c = np.ones((4*nc,nc))
        for i in range(0,4*nc):
            for j in range(0,nc):
                if math.trunc(i/4) == j:
                    one4c[i][j] = 1
                else:
                    one4c[i][j] = 0
        f_CPA1 = np.log(X)-0.5*X+0.5
        f_CPA = np.dot(np.dot(one4c,x),f_CPA1)
    
    f = {
        2: -rho*R*T*np.log(1-rho*bmix)-rho*amix/bmix*np.log(1+rho*bmix), #SRK+RG
        4: -rho*R*T*np.log(1-rho*bmix)-rho*amix/bmix*np.log(1+rho*bmix), #PR+RG
        6: -rho*R*T*np.log(1-rho*bmix)-rho*amix/bmix*np.log(1+rho*bmix)+rho*R*T*f_CPA #CPA+RG
    }.get(EoS,'NULL')
    
    return f
#=========================================================================================

#Calculates Helmholtz long range forces----------------------------------------------------
def helm_long(EoS,rho,f):
    
    f0a = {
        2: -0.5*rho**2, #SRK+RG
        4: -0.5*rho**2, #PR+RG
        6: -0.5*rho**2 #CPA+RG
    }.get(EoS,'NULL')
    
    flong = f - f0a 
    return flong
#=========================================================================================

#Calculates Helmholtz short range forces--------------------------------------------------
def helm_short(EoS,rho,f,phi,i):
    
    f0a = {
        2: -0.5*rho**2, #SRK+RG
        4: -0.5*rho**2, #PR+RG
        6: -0.5*rho**2 #CPA+RG
    }.get(EoS,'NULL')
    
    fshort = {
        2: f - f0a*phi/(2**i),      #SRK+RG
        4: f - f0a*phi/(2**i),      #PR+RG
        6: f - f0a*phi/(2**(2*i+1)) #CPA+RG
    }.get(EoS,'NULL')
    
    return fshort
#=========================================================================================

#Calculates Change in Helmholtz energy density--------------------------------------------
def renorm_df(w,n,fl,fs,K,rho,width):

    suml = 0
    sums = 0
    t = 0
    aminl = 0
    amins = 0
    
    Gl = np.ones((n))
    Gs = np.ones((n))
    Gl2 = np.ones((n))
    Gs2 = np.ones((n))
    argl = np.ones((n))
    args = np.ones((n))
    
    while t<min(w+1,n-w):
        Gl[t] = (fl[w+t] - 2*fl[w] + fl[w-t])/2
        Gs[t] = (fs[w+t] - 2*fs[w] + fs[w-t])/2
        Gl2[t] = np.exp(-Gl[t]/K)
        Gs2[t] = np.exp(-Gs[t]/K)
        argl[t] = Gl[t]/K
        args[t] = Gs[t]/K
        if argl[t]<aminl:
            aminl = argl[t]
        if args[t]<amins:
            aminl = args[t]
        t=t+1

    if w<(n/2):
        suml = 0.25*(np.exp(-argl[0]+aminl)+np.exp(-argl[w]+aminl))
        sums = 0.25*(np.exp(-args[0]+amins)+np.exp(-args[w]+amins))
    else:
        suml = 0.5*(np.exp(-argl[0]+aminl)+np.exp(-argl[n-w]+aminl))
        sums = 0.5*(np.exp(-args[0]+amins)+np.exp(-args[n-w]+amins))
    
    t=1
    while t<min(w,n-w):
        along = argl[t] - aminl
        ashort = args[t] - amins
        if along<30:
            suml = suml + 2.0*np.exp(-along)
        if ashort<30:
            sums = sums + 2.0*np.exp(-ashort)
        t=t+1

    Inl = np.log(width*suml)-aminl
    Ins = np.log(width*sums)-amins

    #using direct trapezoidal integration
    if w<n/2:
        Inl = numerical.trapezoidal(rho,Gl2,0,w)
        Ins = numerical.trapezoidal(rho,Gs2,0,w)
    else:
        Inl = numerical.trapezoidal(rho,Gl2,0,n-w)
        Ins = numerical.trapezoidal(rho,Gs2,0,n-w)
    
    Inl = np.log(Inl)
    Ins = np.log(Ins)
    df = -K*(Ins-Inl)

    if math.isnan(df) or math.isinf(df):
        df = 1e-15

    return df
#=========================================================================================

#Report f and rho value-------------------------------------------------------------------
def report_f_pure(data,options,title,print_options):
    n = len(data[0])
    f = np.array(menus.flatten(data[0]))
    x = np.array(menus.flatten(data[1]))
    rho = np.array(menus.flatten(data[2]))
    header = str(';%s\n' % ';'.join(rho))
    savedir = str('../output/%s' %title)
    with open(savedir,'w') as file:
        file.write(header)
        for i in range(0,n):
                lin1 = [str(round(rho[i],9)),str(round(f[i],9))]
                lin = ('%s\n' % ';'.join(lin1))
                file.write(lin)
#=========================================================================================