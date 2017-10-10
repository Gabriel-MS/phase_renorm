import numpy as np
import menus
import data
import eos
import association
import math
import numerical

import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline, splrep, splev, interp1d

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
    stepx = (1/float(nx)) #Step to calculate change
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
    fmat = []
    umat = []
    ures = np.empty((nd))
    uv = []
    
    if nc==1:
        X = np.ones((8))
    
    #Main loop*************************************
    while x[0]<1.0:
        print x[0]
        if x[0]==0.006: #after first step
            x[0]=0.005
            x[1]=1-x[0]
        
        if nc==1:
            x[0] = 0.999999
            x[1] = 0.000001
        
        #Mixture parameters
        bmix = eos.bmix_calc(MR,b,x)
        amix = eos.amix_calc(MR,a,x,kij)
        rhomax = 0.999999
        
        #Mixture Renormalization parameters
        L = np.dot(x,np.power(L_rg,3.0))
        L = np.power(L,1.0/3.0)
        phi = np.dot(x,phi_rg)
        
        for k in range(0,nd):
            rho[k] = np.array(float(k)/nd/bmix)
            if k==0:
                rho[0] = 1e-6
            if EoS==6:
                Xf = association.frac_nbs(nc,1/rho[k],CR,en_auto,beta_auto,b,bmix,X,0,x,0,T,SM)
                X = Xf
            f[k] = np.array(helm_rep(EoS,R,T,rho[k],amix,bmix,X,x,nc))   #Helmholtz energy density
            f_orig[k] = f[k]                                #Initial helmholtz energy density
            
            #Subtract attractive forces (due long range correlations)
            f[k] = f[k] + amix*(rho[k]**2)
            k = k+1
            
            
        #Adimensionalization
        rho = rho*bmix
        f = f*bmix*bmix/amix
        T = T*bmix*R/amix
    

        rho1 = rho.flatten()

        #Main loop****************************************************************
        i = 1
        while i<n:
            #K = kB*T/((2**(3*i))*(L**3))
            K = T/(2**(3*i))/((L**3)/bmix*6.023e23)
            
            #Long and Short Range forces
            fl = helm_long(EoS,rho,f)
            fs = helm_short(EoS,rho,f,phi,i)

            #Calculate df
            width = rhomax/nd
            w = 0
            for w in range(0,nd):
                df[w] = renorm_df(w,nd,fl,fs,K,rho,width)

            #Update Helmholtz Energy Density
            df = np.array(df)
            f = f + df
            i = i+1

        #Dimensionalization
        rho = rho/bmix
        f = f/bmix/bmix*amix
        T = T/bmix/R*amix
        
        #Add original attractive forces
        f = f - amix*(rho**2)
        
        #Store residual value of f, calculate helmholtz energy density adding ideal gas energy
        fres = f
        f = f + rho*R*T*(np.log(rho)-1)

        if(EoS==6):
			f = fres
        
        fv.append(f)
        fresv.append(fres)
        x0v.append(x[0])
        
        if r==0:
            rhov.append(rho) #rho vector is always the same
        r=1

        fmat.append(f)

        if nc>1:
            drho = rho[nd/2]-rho[nd/2-1]
            for i in range(1,nd-2):
                ures[i] = (fres[i+1]-fres[i-1])/(2*drho)
            ures[nd-1] = (fres[nd-1]-fres[nd-2])/drho
            ures[0] = (fres[1]-fres[0])/drho
            uv.append(ures)
            umat.append(uv)

        x[0] = x[0]+stepx
        
    renorm_out = []
    renorm_out.append(fv)
    renorm_out.append(x0v)
    renorm_out.append(rhov)
    renorm_out.append(fresv)
    renorm_out.append(fmat)
    renorm_out.append(umat)
    if nc>1: #If binary mixture, report calculated values
        print 'before report'
        report_renorm_bin(rhov,x0v,fmat,nx,nd,MR,IDs,EoS)
    return renorm_out
#=========================================================================================

#Output renorm values of binary mixtures--------------------------------------------------
def report_renorm_bin(rhov,x0v,fmat,nx,nd,MR,IDs,EoS):

    x = x0v
    x1 = np.array(menus.flatten(x))
    x2 = 1.0-x1
    rhob = rhov[0]
    f2d = np.empty((nx,nd))
    rho = np.empty((nx,nd))
    dfdrhob = np.empty((nx,nd))
    dfdx1 = np.empty((nx,nd))
    dfdx2 = np.empty((nx,nd))
    rho1 = np.empty((nx,nd))
    rho2 = np.empty((nx,nd))
    u1mat = np.empty((nx,nd))
    u2mat = np.empty((nx,nd))

    bmix = np.empty((nx))
    b = eos.b_calc(IDs,EoS)
    b1 = b[0]
    b2 = b[1]
    for i in range(0,nx):
        bmix[i] = eos.bmix_calc(MR,b,x[i])

    for i in range(0,nx):
        f = np.array(fmat[i])
        for j in range(0,nd):
            f2d[i][j] = f[j]
            rho[i][j] = rhob[j]/bmix[i]
    
    for i in range(0,nx):
        for j in range(0,nd):
            rho1[i][j] = x1[i]*rho[i][j]
            rho2[i][j] = x2[i]*rho[i][j]

    for i in range(0,nx):
        for j in range(0,nd-1):
            dfdrhob[i][j] = (f2d[i][j+1]-f2d[i][j])/(rhob[j+1]-rhob[j])
        dfdrhob[i][nd-1] = (f2d[i][nd-1]-f2d[i][nd-2])/(rhob[nd-1]-rhob[nd-2])

    for i in range(0,nx-1):
        for j in range(0,nd):
            dfdx1[i][j] = (f2d[i+1][j]-f2d[i][j])/(x1[i+1]-x1[i])
        dfdx1[nx-1][j] = (f2d[nx-1][j]-f2d[nx-2][j])/(x1[nx-1]-x1[nx-2])

    for i in range(0,nx-1):
        for j in range(0,nd):
            dfdx2[i][j] = (f2d[i+1][j]-f2d[i][j])/(x2[i+1]-x2[i])
        dfdx2[nx-1][j] = (f2d[nx-1][j]-f2d[nx-2][j])/(x2[nx-1]-x2[nx-2])

    for i in range(0,nx):
        for j in range(0,nd):
            u1mat[i][j] = dfdrhob[i][j]*b1+dfdx2[i][j]*(0-rho2[i][j])/(rho[i][j]*rho[i][j])
            u2mat[i][j] = dfdrhob[i][j]*b2+dfdx1[i][j]*(0-rho1[i][j])/(rho[i][j]*rho[i][j])

    #u1res_mat = RectBivariateSpline(x,rhob,u1mat)
    #u2res_mat = RectBivariateSpline(x,rhob,u2mat)

    u1mat_r = []
    u2mat_r = []
    for i in range(0,nx):
        u1mat_r.append(u1mat[i])
        u2mat_r.append(u2mat[i])

    title = 'ren_f.tmp'
    savedir = str('%s' %title)
    with open(savedir,'w') as file:
        dv = str("")
        for j in range(0,nd):
            d1 = str(round(rhob[j],9))
            dv = dv+';'+d1
        file.write(dv)
        file.write('\n')
        
        for i in range(0,nx):
            x = str(round(x0v[i],9))
            f = fmat[i]
            lin1 = x
            for j in range(0,nd):
                f1 = str(round(f[j],9))
                lin1 = lin1+';'+f1
            file.write(lin1)
            file.write('\n')

    title = 'u1.tmp'
    savedir = str('%s' %title)
    with open(savedir,'w') as file:
        dv = str("")
        for j in range(0,nd):
            d1 = str(round(rhob[j],9))
            dv = dv+';'+d1
        file.write(dv)
        file.write('\n')
        
        for i in range(0,nx):
            x = str(round(x0v[i],9))
            u = u1mat_r[i]
            lin1 = x
            for j in range(0,nd):
                u1 = str(round(u[j],9))
                lin1 = lin1+';'+u1
            file.write(lin1)
            file.write('\n')

    title = 'u2.tmp'
    savedir = str('%s' %title)
    with open(savedir,'w') as file:
        dv = str("")
        for j in range(0,nd):
            d1 = str(round(rhob[j],9))
            dv = dv+';'+d1
        file.write(dv)
        file.write('\n')
        
        for i in range(0,nx):
            x = str(round(x0v[i],9))
            u = u2mat_r[i]
            lin1 = x
            for j in range(0,nd):
                u2 = str(round(u[j],9))
                lin1 = lin1+';'+u2
            file.write(lin1)
            file.write('\n')


#=========================================================================================

#Calculates Fugacity coefficient of a component in a mixtures with renormalization-------
def fugacity_renormalized(phase, x, P, bmix, T, V, ures_Bspln):

    #Calculate ures
    ures = ures_Bspln(x,(1/V)*bmix) #bmix is here because in the spline, rho is adimensional

    #Calculate Z
    Z = P*V/(R*T)

    #Calculate phi
    lnphi = ures/R/T - np.log(Z)
    phi = np.exp(lnphi)

    return phi
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
        6: rho*R*T*(np.log(rho/(1-rho*bmix))-1)-rho*amix/bmix*np.log(1+rho*bmix)+rho*R*T*f_CPA #CPA+RG
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
    #print 'inside renorm_df',w,n
    suml = 0
    sums = 0
    t = 0
    aminl = 0
    amins = 0
    Gl2 = np.empty((n))
    Gs2 = np.empty((n))
    argl = np.zeros((n))
    args = np.zeros((n))
    
    while t<min(w+1,n-w):
        Gl2[t] = np.exp(-((fl[w+t] - 2*fl[w] + fl[w-t])/2)/K)
        Gs2[t] = np.exp(-((fs[w+t] - 2*fs[w] + fs[w-t])/2)/K)
        t=t+1

    Inl = numerical.trapezoidal(rho,Gl2,0,t) #era w
    Ins = numerical.trapezoidal(rho,Gs2,0,t) #era w
    df = -K*np.log(Ins/Inl)

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

#Outputs renormalized csv file of given EoS at give T-------------------------------------
def renorm_est(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,phi,L):
    #nd    Size of density grid
    #nx    Size of mole fraction grid
    #n     Main loop iteration controller
    
    #If only 1 component is present mimic a binary mixture made of the same component
    if nc==1:
        IDs[1] = IDs[0]
    
    #Recover parameters
    #L_rg = data.L(IDs)         #Vector with L parameters (cutoff length)
    #phi_rg = data.phi(IDs)     #Vector with phi parameters
    Tc = data.Tc(IDs)
    
    #Components parameters
    a = eos.a_calc(IDs,EoS,T)
    b = eos.b_calc(IDs,EoS)
    Tr = T/np.array(Tc)
    
    #Main loop parameters
    x = np.array([0.0001,0.9999])
    stepx = (1/float(nx)) #Step to calculate change
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
            x[0] = 0.999999
            x[1] = 0.000001
        
        #Mixture parameters
        bmix = eos.bmix_calc(MR,b,x)
        amix = eos.amix_calc(MR,a,x,kij)
        rhomax = 0.999999
        
        #Mixture Renormalization parameters
        #L = np.dot(x,np.power(L_rg,3.0))
        #L = np.power(L,1.0/3.0)
        #phi = np.dot(x,phi_rg)
        
        for k in range(0,nd):
            rho[k] = np.array(float(k)/nd/bmix)
            if k==0:
                rho[0] = 1e-6
            if EoS==6:
                Xf = association.frac_nbs(nc,1/rho[k],CR,en_auto,beta_auto,b,bmix,X,0,x,0,T,SM)
                X = Xf
            f[k] = np.array(helm_rep(EoS,R,T,rho[k],amix,bmix,X,x,nc))   #Helmholtz energy density
            f_orig[k] = f[k]                                #Initial helmholtz energy density
            
            #Subtract attractive forces (due long range correlations)
            f[k] = f[k] + amix*(rho[k]**2)
            k = k+1
            
            
        #Adimensionalization
        rho = rho*bmix
        f = f*bmix*bmix/amix
        T = T*bmix*R/amix
    

        rho1 = rho.flatten()

        #Main loop****************************************************************
        i = 1
        while i<n:
            #K = kB*T/((2**(3*i))*(L**3))
            K = T/(2**(3*i))/((L**3)/bmix*6.023e23)
            
            #Long and Short Range forces
            fl = helm_long(EoS,rho,f)
            fs = helm_short(EoS,rho,f,phi,i)

            #Calculate df
            width = rhomax/nd
            w = 0
            for w in range(0,nd):
                df[w] = renorm_df(w,nd,fl,fs,K,rho,width)

            #Update Helmholtz Energy Density
            df = np.array(df)
            f = f + df
            i = i+1

        #Dimensionalization
        rho = rho/bmix
        f = f/bmix/bmix*amix
        T = T/bmix/R*amix
        
        #Add original attractive forces
        f = f - amix*(rho**2)
        
        #Store residual value of f, calculate helmholtz energy density adding ideal gas energy
        fres = f
        f = f + rho*R*T*(np.log(rho)-1)

        if(EoS==6):
			f = fres
        
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
