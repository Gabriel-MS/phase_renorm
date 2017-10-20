import numpy as np
import menus
import data
import eos
import association
import math
import numerical
from scipy.interpolate import InterpolatedUnivariateSpline, splrep, splev, interp1d
from scipy.interpolate import RectBivariateSpline

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
    x = np.array([0.000001,0.999999])
    stepx = (1/float(nx)) #Step to calculate change
    k = 0               #Vector fill counter
    i = 1               #Main loop counter
    r = 0               #Report counter
    count = 0
    rho = np.empty((nd))            #Density vector
    rhov = []                       #Density vector to export
    x0v = []                        #Mole fraction vector to export
    f = np.empty((nd))              #Helmholtz energy density vector
    fv = []                         #Helmholtz energy density vector to export
    fresv = []                      #Residual Helmholtz energy density vector to export
    Tv = []                         #Temperature values to export
    df = np.empty((nd))             #Changes in helmholtz energy density vector
    f_orig = np.empty((nd))         #Unmodified Helmholtz energy density vector
    rhob = []                       #Adimensional density vector
    u = np.empty((nd))
    X = np.ones((4*nc))
    Pv = []
    fmat = []
    Pmatv = np.empty((nx,nd))
    fmatres = []
    umat = []
    ures = np.empty((nd))
    uv = []
    
    if nc==1:
        X = np.ones((8))
    
    #Main loop*************************************
    while x[0]<1.0:
        if nc>1:
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
                X = association.frac_nbs(nc,1/rho[k],CR,en_auto,beta_auto,b,bmix,X,0,x,0,T,SM)
            f[k] = np.array(helm_rep(EoS,R,T,rho[k],amix,bmix,X,x,nc))   #Helmholtz energy density
            f_orig[k] = f[k]                                #Initial helmholtz energy density
            
            #Subtract attractive forces (due long range correlations)
            f[k] = f[k] + 0.5*amix*(rho[k]**2)
            k = k+1
            
        #Adimensionalization
        rho = rho*bmix
        f = f*bmix*bmix/amix
        T = T*bmix*R/amix

        rho1 = rho.flatten()

        #Main loop****************************************************************
        i = 1
        while i<=n:
            #print i
            #K = kB*T/((2**(3*i))*(L**3))
            #K = R*T/((L**3)*(2**(3*i)))
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
            #print 'i=',i,K/bmix/bmix*amix,f[60]/bmix/bmix*amix,df[60]/bmix/bmix*amix,fl[60]/bmix/bmix*amix,fs[60]/bmix/bmix*amix
            i = i+1

        #Dimensionalization
        rho = rho/bmix
        f = f/bmix/bmix*amix
        T = T/bmix/R*amix
        
        #Add original attractive forces
        f = f - 0.5*amix*(rho**2)
        
        #Store residual value of f
        fres = f - rho*R*T*(np.log(rho)-1)
        #f = f + rho*R*T*(np.log(rho)-1) #Already accounting ideal gas energy

        #if(EoS==6):
        #    f = fres
        
        fv.append(f)
        fresv.append(fres)
        x0v.append(x[0])
        
        if r==0:
            rhob.append(rho*bmix) #rhob vector is always the same
            rhov.append(rho) #in case the calculation is done for one-component
        r=1

        drho = rho[int(nd/2)]-rho[int(nd/2)-1]
        for i in range(1,nd-2):
            u[i] = (f[i+1]-f[i-1])/(2*drho)
        u[nd-1] = (f[nd-1]-f[nd-2])/drho
        u[0] = (f[1]-f[0])/drho

        P = -f+rho*u
        Pv.append(P)
        for j in range(0,nd):
            Pmatv[count][j] = P[j]
            #print Pmatv[count][j],count,j,x[0]
        count = count+1

        fmat.append(f)
        fmatres.append(fres)

        x[0] = x[0]+stepx
        
    if nc>1:
        Pmat = RectBivariateSpline(x0v,rhob,Pmatv)
    else:
        Pmat = 'NULL'

    renorm_out = []
    renorm_out.append(fv)
    renorm_out.append(x0v)
    renorm_out.append(rhov)
    renorm_out.append(rhob)
    renorm_out.append(fmat)
    renorm_out.append(Pmat)
    if nc>1: #If binary mixture, report calculated values
        print 'before report'
        ren_u = report_renorm_bin(rhob,x0v,fmatres,nx,nd,MR,IDs,EoS)
        renorm_out.append(ren_u)
    else:
        renorm_out.append(0)
    renorm_out.append(fresv)
    renorm_out.append(Pv)
    return renorm_out
#=========================================================================================

#Output renorm values of binary mixtures--------------------------------------------------
def report_renorm_bin(rhov,x0v,fm,nx,nd,MR,IDs,EoS):

    rep = []
    x = x0v
    x1 = np.array(menus.flatten(x))
    x2 = 1.0-x1
    rhob = rhov[0]
    xb = np.empty((2))
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
        xb[0] = x1[i]
        xb[1] = x2[i]
        bmix[i] = eos.bmix_calc(MR,b,xb)

    for i in range(0,nx):
        f = np.array(fm[i])
        for j in range(0,nd):
            f2d[i][j] = f[j]
            rho[i][j] = rhob[j]/bmix[i]
    
    for i in range(0,nx):
        for j in range(0,nd):
            rho1[i][j] = x1[i]*rho[i][j]
            rho2[i][j] = x2[i]*rho[i][j]

    for i in range(0,nx):
        dfdrhob[i][0] = (f2d[i][1]-f2d[i][0])/(rhob[1]-rhob[0])
        for j in range(1,nd-1):
            dfdrhob[i][j] = (f2d[i][j+1]-f2d[i][j-1])/(rhob[j+1]-rhob[j-1])
        dfdrhob[i][nd-1] = (f2d[i][nd-1]-f2d[i][nd-2])/(rhob[nd-1]-rhob[nd-2])

    for i in range(0,nx-1):
        for j in range(0,nd):
            if i!=0:
                dfdx1[i][j] = (f2d[i+1][j]-f2d[i-1][j])/(x1[i+1]-x1[i-1])
            else:
                dfdx1[0][j] = (f2d[1][j]-f2d[0][j])/(x1[1]-x1[0])
        dfdx1[nx-1][j] = (f2d[nx-1][j]-f2d[nx-2][j])/(x1[nx-1]-x1[nx-2])

    for i in range(0,nx-1):
        for j in range(0,nd):
            if i!=0:
                dfdx2[i][j] = (f2d[i+1][j]-f2d[i-1][j])/(x2[i+1]-x2[i-1])
            else:
                dfdx2[0][j] = (f2d[1][j]-f2d[0][j])/(x2[1]-x2[0])
        dfdx2[nx-1][j] = (f2d[nx-1][j]-f2d[nx-2][j])/(x2[nx-1]-x2[nx-2])

    for i in range(0,nx):
        for j in range(0,nd):
            u1mat[i][j] = dfdrhob[i][j]*b1+dfdx2[i][j]*(0-rho2[i][j])/(rho[i][j]*rho[i][j])
            u2mat[i][j] = dfdrhob[i][j]*b2+dfdx1[i][j]*(0-rho1[i][j])/(rho[i][j]*rho[i][j])

    u1res_mat = RectBivariateSpline(x,rhob,u1mat)
    u2res_mat = RectBivariateSpline(x,rhob,u2mat)

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
            f = fm[i]
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

    rep.append(u1res_mat)
    rep.append(u2res_mat)
    return rep
#=========================================================================================

#Calculates Fugacity coefficient of a component in a mixtures with renormalization-------
def lnfugcoef_renorm(P,T,x,phase,V,r_data,bmix):

    rhob = (1/V)*bmix
    size = x.shape[0]
    ures = np.empty((size))

    #Calculate ures
    uspln = r_data[6]
    u1res_Bspln = uspln[0]
    u2res_Bspln = uspln[1]
    u1res = u1res_Bspln(x[0],rhob) #bmix is here because rho is adimensional inside the spline
    u2res = u2res_Bspln(x[0],rhob) #bmix is here because rho is adimensional inside the spline, IT REALLY IS x[0] HERE
    ures[0] = u1res
    ures[1] = u2res

    #Calculate Z
    Z = P*V/(R*T)

    #Calculate phi
    lnphi = ures/R/T - np.log(Z)
    #print 'ures=',ures,x[0],x[1],rhob
    #input('...')

    return lnphi
#=========================================================================================

#Calculates Volume of a component in phase of a mixture with renormalization--------------
def volume_renorm(phase, xint, Pint, bmix, R, T, r_data):

    Pspln = r_data[5]
    rho = r_data[3][0]
    x = np.array(menus.flatten(r_data[1]))

    nd = len(rho)
    nx = x.shape[0]

    Pvec = np.empty((nd))
    Pfvec = np.empty((nd))
    dPdrho = np.empty((nd))

    flag = False
    inflex = False

    while flag!=True:
        #Interpolate specific pressure
        Pvint = Pspln(xint,0.0001)
        Pvec[0] = Pvint
        Pfvec[0] = Pvec[0] - Pint

        for i in range(1,nd-1):
            rhoint = float(i)/nd
            Pvint = Pspln(xint,rhoint)
            #print i,rhoint,Pvint
            Pvec[i] = Pvint
            Pfvec[i] = Pvec[i] - Pint
            dPdrho[i] = (Pvec[i+1] - Pvec[i-1])/(float(i+1)/nd-float(i-1)/nd)
            if inflex==False and dPdrho[i]<0:
                inflex=True
        Pvint = Pspln(xint,int(nd-1))
        Pvec[nd-1] = Pvint
        Pfvec[nd-1] = Pvec[nd-1] - Pint
        dPdrho[0] = (Pvec[1] - Pvec[0])/(float(1)/nd-float(0)/nd)
        dPdrho[nd-1] = (Pvec[nd-1] - Pvec[nd-2])/(float(nd-1)/nd-float(nd-2)/nd)

        #Bracketing the real densities at given P
        #print Pvec
        #print Pfvec
        #plt.plot(rho,Pvec)
        #plt.ylim(-15,15)
        #plt.show()
        max1 = 2
        min1 = int(0.90*nd)
        max2 = max1+2
        min2 = min1-2
        #print rho[max1],rho[max2]
        while Pfvec[max1]*Pfvec[max2]>0:
            #max2 = max2+int(nd/200)
            max2 = max2+2
            #print 'max',max2
        if max2-int(nd/100)<0:
            max1 = 0
            #print 'max1',max1
        else:
            #max1 = max2-int(nd/100)
            max1 = max2-2
            #print 'else',max1

        while Pfvec[min1]*Pfvec[min2]>0:
            #min2 = min2-int(nd/200)
            min2 = min2-2
        #min1 = min2+int(nd/100)
        min1 = min2+2

        #print 'falsi_spline',rho[max1],rho[max2],rho[min1],rho[min2]
        #Calculate coexistence densities in interpolated isotherm for given P
        rho_vap = numerical.falsi_spline(rho, Pfvec, rho[max1], rho[max2], 1e-5)
        rho_liq = numerical.falsi_spline(rho, Pfvec, rho[min2], rho[min1], 1e-5)

        if inflex==True and abs(rho_vap-rho_liq)<1e-5:
            Pint=Pint/2

        if inflex==True and abs(rho_vap-rho_liq)>1e-5:
            flag=True

        if xint>0.05:
            flag=True

    #Select desired density
    if phase<0:
        rho_out = rho_liq
    if phase>0:
        rho_out = rho_vap

    V = 1/(rho_out/bmix)
    V_out = []
    V_out.append(V)
    V_out.append(0)
    return V_out
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

    #Considering ideal gas energy contribution
    f = {
        2: -rho*R*T*np.log(1-rho*bmix)-rho*amix/bmix*np.log(1+rho*bmix)+rho*R*T*(np.log(rho)-1), #SRK+RG
        4: -rho*R*T*np.log(1-rho*bmix)-rho*amix/bmix*np.log(1+rho*bmix)+rho*R*T*(np.log(rho)-1), #PR+RG
        6: rho*R*T*(np.log(rho/(1-rho*bmix))-1)-rho*amix/bmix*np.log(1+rho*bmix)+rho*R*T*f_CPA #CPA+RG
    }.get(EoS,'NULL')
    
    #if abs(rho*bmix-0.15)<1e-4:
    #    print f_CPA1,f_CPA,f
    
    return f
#=========================================================================================

#Calculates Helmholtz long range forces----------------------------------------------------
def helm_long(EoS,rho,f):
    
    f0a = {
        2: -0.5*rho*rho, #SRK+RG
        4: -0.5*rho*rho, #PR+RG
        6: -0.5*rho*rho #CPA+RG
    }.get(EoS,'NULL')
    
    flong = f - f0a
    return flong
#=========================================================================================

#Calculates Helmholtz short range forces--------------------------------------------------
def helm_short(EoS,rho,f,phi,i):
    
    f0a = {
        2: -0.5*rho*rho, #SRK+RG
        4: -0.5*rho*rho, #PR+RG
        6: -0.5*rho*rho #CPA+RG
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
    t = 1
    Inl = 0
    Ins = 0
    Gl2 = np.zeros((n))
    Gs2 = np.zeros((n))
    
    #while t<min(w+1,n-w):
    while t<min(w,n-w):
        Gl2[t] = np.exp(-((fl[w+t] - 2*fl[w] + fl[w-t])/2)/K)
        Gs2[t] = np.exp(-((fs[w+t] - 2*fs[w] + fs[w-t])/2)/K)
        t=t+1

    Inl = numerical.trapezoidal(rho,Gl2,0,t-1) #era w
    Ins = numerical.trapezoidal(rho,Gs2,0,t-1) #era w
    
    logInl = np.log(Inl)
    logIns = np.log(Ins)
    df = -K*(logIns-logInl)
    #print Inl
    #if w==60:
    #    print w,t,rho[60],Ins,Inl,df
    #    raw_input('......')
        
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
