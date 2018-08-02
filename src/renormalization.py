import numpy as np
import menus
import data
import eos
import association
import math
import PSO
import numerical
import envelope
import time
from scipy.interpolate import InterpolatedUnivariateSpline, splrep, splev, interp1d
from scipy.interpolate import RectBivariateSpline
from scipy import stats
from scipy.linalg import inv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = 8.314462175e-6 #m3.MPa/K/mol
NA = 6.023e23      #Avogadro Number
kB = 1.3806503e-23 #Boltzmann constant, J/K

#Outputs renormalized csv file of given EoS at give T-------------------------------------
def renorm(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate,L_est,phi_est):
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
    x = np.array([0.0001,0.9999])
    stepx = (1/float(nx)) #Step to calculate change
    k = 0               #Vector fill counter
    i = 1               #Main loop counter
    r = 0               #Report counter
    count = 0
    rho = np.empty((nd))            #Density vector
    rhov = []                       #Density vector to export
    x0v = []                        #Mole fraction vector to export
    bmixv = []
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
    
    df_vec = []
    f_vec2 = []
    
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
        Nav = 6.023e23
        rhomax = 0.999999
        
        #Mixture Renormalization parameters
        L = np.dot(x,np.power(L_rg,3.0))
        L = np.power(L,1.0/3.0)
        phi = np.dot(x,phi_rg)
        
        #print L
        #print phi
        
        
        pi = math.pi
        #sig = np.power(6/pi*b/Nav,1.0/3.0)[0]
        #sig = np.power(b/Nav,1.0/3.0)[0]
        #c1 = data.c1(IDs)[0]
        #en = data.en(IDs)[0]
        #sig = np.power(1.15798*b/Nav,1.0/3.0)[0]
        #L = 1.25*sig
        #L = 1.5*sig
        #L = 1/c1*sig
        #print L,phi
        #L = 0.5/c1*sig
        #PHI = 4*(pi**2.0)
        
        #PHI = 1.0/pi/4.0
        #lamda = 1.5
        #w_LJ = (9.0*sig/7.0) #lennard-jones
        #print 'LJ=',w_LJ
        #w_SW = np.sqrt((1./5.)*(sig**2.0)*(lamda**5.0-1)/(lamda**3.0-1)) #square-well potential
        #print 'SW=',w_SW
        #phi = PHI*(w_LJ**2)/2/(L**2)
        #phi = PHI*(w_SW**2)/2/(L**2)
        
        #om = data.omega(IDs)
        #phi = 2/np.power(np.exp(om),4)[0]
        #w = 0.575*sig*en/T/kB/b[0]*1e6
        #print 'w=',w
        #phi = 2/np.power(np.exp(c1),4)[0]
        
        #w = 100.0*1e-9/100 #van der waals wavelength 100nm
        #phi = PHI*(w**2)/2/(L**2)
        
        #print L
        #print phi
        #print '---------'
        
        #New parameters
        #L = 1.5*np.power(b/Nav,1.0/3.0)
        #h = 6.626e-34
        #kkB = 1.38e-23
        #MM = 0.034
        #deBroglie = h/np.sqrt(3*kkB*T*MM/Nav)
        #phi = (deBroglie**2.0)/(L**2.0)*150*3.14
        #L = L[0]
        #phi = phi[0]
        #print 'L=',L
        #print 'phi=',phi
        

        if estimate==True:
            L = L_est
            phi = phi_est

        for k in range(0,nd):
            rho[k] = np.array(float(k)/nd/bmix)
            if k==0:
                rho[0] = 1e-6
            if EoS==6:
                if k==0:
                    X = association.frac_nbs(nc,1/rho[k],CR,en_auto,beta_auto,b,bmix,X,0,x,0,T,SM)
                else:
                    X = association.frac_nbs(nc,1/rho[k],CR,en_auto,beta_auto,b,bmix,X,1,x,0,T,SM)
            #print X,k
            #raw_input('...')
            f[k] = np.array(helm_rep(EoS,R,T,rho[k],amix,bmix,X,x,nc))   #Helmholtz energy density
            k = k+1
            
        f_orig = f                                #Initial helmholtz energy density
    
        #Subtract attractive forces (due long range correlations)
        f = f + 0.5*amix*(rho**2)
        
        df_vec.append(rho)

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
            #df = np.array(df) #used to evaluate each step
            #df_vec.append(df)
            f = f + df
            #f_vec2.append(f)
            #print 'i=',i,K/bmix*amix,f[60]/bmix*amix,df[60]/bmix*amix,T
            i = i+1

        #Dimensionalization
        rho = rho/bmix
        f = f/bmix/bmix*amix
        T = T/bmix/R*amix
        
        #df_total = 
        #df = np.array(df)
        #df_vec.append(df)
        
        #Add original attractive forces
        f = f - 0.5*amix*(rho**2)
        
        #Store residual value of f
        #fres = f - rho*R*T*(np.log(rho)-1) #WRONG
        fres = f - rho*R*T*np.log(rho)
        #f = f + rho*R*T*(np.log(rho)-1) #Already accounting ideal gas energy
        
        #strT = str(T)
        #dfT = ('df_%s.csv' %strT)
        #envelope.report_df(df_vec,'df.csv')
        #envelope.report_df(f_vec2,'f.csv')

        #if(EoS==6):
        #    f = fres
        
        fv.append(f)
        fresv.append(fres)
        x0v.append(x[0])
        bmixv.append(bmix)
        
        if r==0:
            rhob.append(rho*bmix) #rhob vector is always the same
            rhov.append(rho) #in case the calculation is done for one-component
        r=1

        drho = rho[int(nd/2)]-rho[int(nd/2)-1]
        for i in range(1,nd-2):
            u[i] = (f[i+1]-f[i-1])/(2*drho)
        u[nd-1] = (f[nd-1]-f[nd-2])/drho
        u[0] = (f[1]-f[0])/drho
        
        fspl = splrep(rho,f,k=3)         #Cubic Spline Representation
        f = splev(rho,fspl,der=0)
        u = splev(rho,fspl,der=1)        #Evaluate Cubic Spline First derivative

        P = -f+rho*u
        Pv.append(P)
        for j in range(0,nd):
            Pmatv[count][j] = P[j]
            #print Pmatv[count][j],count,j,x[0]
        count = count+1

        fmat.append(f)
        fmatres.append(fres)

        x[0] = x[0]+stepx
        #if nc>1:
        #    if abs(x[0]-1.0)<1e-5:
        #        x[0] = 0.9999
        
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
    renorm_out.append(bmixv)
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
        #print rho
        #print Pvec
        #print Pfvec
        
        #NEW
        P_fvec_spl = splrep(rho,Pfvec,k=3)         #Cubic Spline Representation
        Pfvec = splev(rho,P_fvec_spl,der=0)
        #NEW
        
        #plt.plot(rho,Pvec)
        #plt.ylim(-15,15)
        #plt.show()
        max1 = 2
        min1 = int(0.90*nd) #it was 0.90 before
        max2 = max1+2
        min2 = min1-2
        #raw_input('before')
        while Pfvec[max1]*Pfvec[max2]>0:# and max2<len(Pfvec):
            #max2 = max2+int(nd/200)
            max2 = max2+1
            #print 'max',max2
            #raw_input('max')
        if max2-int(nd/100)<0:
            max1 = 0
            #print 'max1',max1
        else:
            #max1 = max2-int(nd/100)
            max1 = max2-4
            #print 'else',max1

        while Pfvec[min1]*Pfvec[min2]>0 and min2>0:
            #min2 = min2-int(nd/200)
            min2 = min2-1
            #print 'min',min2
            #raw_input('min')
        #min1 = min2+int(nd/100)
        min1 = min2+4

        #print 'int',Pint,xint,phase
        #print 'falsi_spline',rho[max1],rho[max2],rho[min1],rho[min2]
        #print 'falsi_pressures',Pfvec[max1],Pfvec[max2],Pfvec[min1],Pfvec[min2]
        #Calculate coexistence densities in interpolated isotherm for given P
        rho_vap = numerical.falsi_spline(rho, Pfvec, rho[max1], rho[max2], 1e-5)
        #print 'rho_vap',rho_vap
        rho_liq = numerical.falsi_spline(rho, Pfvec, rho[min2], rho[min1], 1e-5)
        #print 'rho_liq',rho_liq
        #raw_input('...')

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
    V = 1/rho
    f = {
        2: rho*R*T*np.log(rho/(1-rho*bmix))-rho*amix/bmix*np.log(1+rho*bmix), #SRK+RG
        #2: -rho*R*T*np.log(1-rho*bmix)-rho*amix/bmix*np.log(1+rho*bmix)+rho*R*T*(np.log(rho)-1), #SRK+RG
        #4: rho*R*T*np.log(rho/(1-rho*bmix))-rho*amix/bmix*np.log(1+rho*bmix), #PR+RG
        #4: -rho*(R*T*np.log(V-bmix)-(amix/(2*np.sqrt(2)*bmix)*(np.log(1-(V+bmix)/np.sqrt(2)/bmix)-np.log(1+(V+bmix)/np.sqrt(2)/bmix)))), #PR+RG
        4: rho*R*T*np.log(rho/(1-rho*bmix))-rho*amix/bmix/np.sqrt(8)*np.log((1+rho*bmix*(1+np.sqrt(2)))/(1+rho*bmix*(1-np.sqrt(2)))), #PR+RG
        #4: rho*R*T*np.log(rho)+rho*R*T*np.log(1-amix*(V-bmix)/(R*T*(V**2+2*bmix*V-bmix**2)))+rho*amix/(2*np.sqrt(2)*bmix)*np.log((V+(1-np.sqrt(2)*bmix))/(V+(1+np.sqrt(2)*bmix))), #PR+RG
        #4: -rho*R*T*np.log(1-rho*bmix)-rho*amix/bmix*np.log(1+rho*bmix)+rho*R*T*(np.log(rho)-1), #PR+RG
        #6: rho*R*T*(np.log(rho/(1-rho*bmix))-1)-rho*amix/bmix*np.log(1+rho*bmix)+rho*R*T*f_CPA #CPA+RG
        6: rho*R*T*np.log(rho/(1-rho*bmix))-rho*amix/bmix*np.log(1+rho*bmix)+rho*R*T*f_CPA #CPA+RG
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
        2: f - f0a*phi/(2**(2*i)),      #SRK+RG
        4: f - f0a*phi/(2**(2*i)),      #PR+RG
        6: f - f0a*phi/(2**(2*i)) #CPA+RG
        #6: f - f0a*phi/(2**(2*i)) #CPA+RG
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

#Calculates critical exponents------------------------------------------------------------
def critical_exponents(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,L__est,phi__est,Tc,rhoc):

    beta = beta_exponent(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,L__est,phi__est,Tc,rhoc)
    #delta = delta_exponent()

    exponents = []
    exponents.append(beta)
    #exponents.append(delta)
    return exponents
#=========================================================================================

#Calculates  beta critical exponent-------------------------------------------------------
def beta_exponent(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,L__est,phi__est,Tc,rhoc):

    print 'calculating beta exponent'
    Tv = []
    Pv = []
    rhov = []
    rhol = []
    y = np.empty((10))
    x = np.empty((10))

    Tvec = np.linspace(0.99*Tc,0.999*Tc,10)

    #Calculate coexisting densities
    for i in range(0,10):
        T = Tvec[i]
        ren = renorm(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,L__est,phi__est)
        dens = envelope.coexistence_dens(ren[2],ren[0])
        Tv.append(T)
        rhov.append(dens[0])
        rhol.append(dens[1])
        Pv.append(dens[2])
        print T,dens[0],dens[1],dens[2]

    #Calculate critical exponent
    for i in range(0,10):
        y[i] = np.log(abs(rhov[i]-rhol[i])/rhoc)
        x[i] = np.log(abs(Tv[i]-Tc)/Tc)

    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    beta = slope
    print 'beta',beta

    return beta
#=========================================================================================

#Calculates delta critical exponent-------------------------------------------------------
def delta_exponent(env,Tc,Pc,rhoc):



    return delta
#=========================================================================================

#Calculate Objective function based on Critical Temperature and Critical Pressure---------
def objFunc_Tc_Pc(par,argss):

    #Parameters
    EoS     = argss[0]
    IDs     = argss[1]
    MR     = argss[2]
    T      = argss[3]
    Tfinal  = argss[4]
    stepT   = argss[5]
    nd      = argss[6]
    nx      = argss[7]
    kij     = argss[8]
    nc      = argss[9]
    CR      = argss[10]
    en_auto = argss[11]
    beta_auto = argss[12]
    SM      = argss[13]
    n       = argss[14]
    estimate_bool = argss[15]
    crit_bool = argss[16]
    expfile = argss[17]

    #Renormalization Parameters
    L__est   = par[0]
    phi__est = par[1]

    #Calculate Critical Point for given parameters
    env       = envelope.PV_findTc3_envelope(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,crit_bool,L__est,phi__est)
    size      = len(env[0])
    Tc_calc   = env[0][size-1]
    rhoc_calc = env[1][size-1]
    Pc_calc   = env[3][size-1]

    #Recover Experimental data for critical temperature
    expT = data.loadexp(expfile[0])
    Tc_exp = np.mean(expT)
    var_Tc_exp = np.var(expT)

    #Recover Experimental data for critical pressure
    expP = data.loadexp(expfile[1])
    Pc_exp = np.mean(expP)
    var_Pc_exp = np.var(expP)

    #Calculate Objective Function
    Fobj_T = (Tc_calc-Tc_exp)**2/var_Tc_exp
    Fobj_P = (Pc_calc-Pc_exp)**2/var_Pc_exp
    Fobj = Fobj_T + Fobj_P

    """
    print '--------------------------------'
    print 'Parameters:',L__est,phi__est
    print 'Critical Temperature:',Tc_calc,Tc_exp
    print 'Critical Pressure:',Pc_calc,Pc_exp
    print 'Objective Function:',Fobj,(Tc_calc-Tc_exp)**2/var_Tc_exp,(Pc_calc-Pc_exp)**2/var_Pc_exp
    print '--------------------------------\n'
    """

    out = []
    out.append(Fobj)
    out.append(Tc_calc)
    out.append(Pc_calc)
    out.append(rhoc_calc)
    out.append(Fobj_T)
    out.append(Fobj_P)
    return out
#=========================================================================================

#Given initial T, using renormalization method, estimate L and phi parameters----------
def Estimate_Parameters(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,expfile,estimate_bool,crit_bool):

    #Parameters for PSO
    nswarm = 5
    nparameter = 2
    ndata = 1

    #Create Particles
    p = np.empty((nswarm,nparameter))
    for i in range(0,nswarm):
        p[i][0] = np.random.uniform(1e-10,9e-10)
        p[i][1] = np.random.uniform(0.01,10)
    print 'particles'
    print p

    #Organize Parameters
    argss = []
    argss.append(EoS)
    argss.append(IDs)
    argss.append(MR)
    argss.append(T)
    argss.append(Tfinal)
    argss.append(stepT)
    argss.append(nd)
    argss.append(nx)
    argss.append(kij)
    argss.append(nc)
    argss.append(CR)
    argss.append(en_auto)
    argss.append(beta_auto)
    argss.append(SM)
    argss.append(n)
    argss.append(estimate_bool)
    argss.append(crit_bool)
    argss.append(expfile)

    #Initialize PSO method
    #PSO.LJPSO(nparameter,ndata,nswarm,objFunc_Tc_Pc,argss,p)

    
    #Initialize Newton method
    phi = 1.0
    L = 5.6e-10
    L_orig = L
    tol = 1e-5
    k = 0
    err = 1.0
    par = np.empty((2))
    par[0] = L
    par[1] = phi
    Fobj = 1.0
    Lold = 1.0
    while err>tol:

        #Calculate objective function
        F = objFunc_Tc_Pc(par,argss)
        Fold = Fobj
        Fobj = F[0]     #Objective Function value
        argss[3] = F[1]-5 #New temperature estimate is near the last critical temperature
        Tc = F[1]
        Pc = F[2]
        rhoc = F[3]

        #Apply correction to parameter L
        derF = (Fobj-Fold)/(L-Lold)
        Lold = L
        L = L - 0.5*Fobj/derF #damping to help
        if k==0:
            L = L_orig-0.1e-10
        par[0] = L
        errL = L-Lold
        err = Fobj

        k = k+1
        print Tc,Pc,rhoc,L+0.5*Fobj/derF,err
    
    """
    #Initialize Full Newton method
    phi = 5.0
    L = 5e-10
    L_orig = L
    phi_orig = phi
    tol = 1e-5
    k = 0
    err = 1.0
    par = np.empty((2))
    Fv = np.empty((2))
    dF = np.empty((2,2))
    invD = np.empty((2,2))
    D = np.empty((2))
    par[0] = L
    par[1] = phi
    Fobj = 1.0
    Fobj_T = 1.0
    Fobj_P = 1.0
    Lold = 1.0
    phiold = 1.0
    while err>tol:

        #Calculate objective function
        F = objFunc_Tc(par,argss)
        Fold_T = Fobj_T
        Fold_P = Fobj_P
        Fobj = F[0]     #Objective Function value
        Fobj_T = F[4]     #Objective Function value for critical temperature
        Fobj_P = F[5]     #Objective Function value for critical pressure
        argss[3] = F[1]-5 #New temperature estimate is near the last critical temperature
        Tc = F[1]
        Pc = F[2]
        rhoc = F[3]

        #Apply correction to parameter L
        derF_TL = (Fobj_T-Fold_T)/(L-Lold)
        derF_Tphi = (Fobj_T-Fold_T)/(phi-phiold)
        derF_PL = (Fobj_P-Fold_P)/(L-Lold)
        derF_Pphi = (Fobj_P-Fold_P)/(phi-phiold)

        #print L,Lold,phi,phiold
        #print Fobj_T,Fobj_P
        #print derF_TL,derF_Tphi,derF_PL,derF_Pphi

        Lold = L
        phiold = phi
        

        Fv[0] = -Fobj_T
        Fv[1] = -Fobj_P

        dF[0][0] = derF_TL
        dF[0][1] = derF_Tphi
        dF[1][0] = derF_PL
        dF[1][1] = derF_Pphi
        print '---dF---'
        print dF

        detD = derF_TL*derF_Pphi-derF_Tphi*derF_PL
        print '---det---'
        print detD
        invD[0][0] = derF_Pphi/detD
        invD[0][1] = -derF_Tphi/detD
        invD[1][0] = -derF_PL/detD
        invD[1][1] = derF_TL/detD
        D[0] = (invD[0][0]*Fv[0])+(invD[0][1]*Fv[1])
        D[1] = (invD[1][0]*Fv[0])+(invD[1][1]*Fv[1])
        print '---D---'
        print D
        
        #D = inv(dF)*Fv

        #L   = L   + D[0]
        #phi = phi + D[1]

        L   = L   + Fobj_T/derF_TL
        phi = phi + Fobj_P/derF_Pphi

        if k==0:
            L = L_orig+0.1e-10
            phi = phi_orig+0.2

        par[0] = L
        par[1] = phi
        errL = L-Lold
        errphi = phi-phiold
        err = Fobj

        k = k+1
        print Tc,Pc,rhoc,Lold,phiold,err
    """
    par = 1
    return par
#====================================================================================== 

#Given initial T, using renormalization method, estimate L and phi parameters----------
def Estimate_Parameters(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,expfile,estimate_bool,crit_bool):

    #Parameters for PSO
    nswarm = 5
    nparameter = 2
    ndata = 1

    #Create Particles
    p = np.empty((nswarm,nparameter))
    for i in range(0,nswarm):
        p[i][0] = np.random.uniform(1e-10,9e-10)
        p[i][1] = np.random.uniform(0.01,10)
    print 'particles'
    print p

    #Organize Parameters
    argss = []
    argss.append(EoS)
    argss.append(IDs)
    argss.append(MR)
    argss.append(T)
    argss.append(Tfinal)
    argss.append(stepT)
    argss.append(nd)
    argss.append(nx)
    argss.append(kij)
    argss.append(nc)
    argss.append(CR)
    argss.append(en_auto)
    argss.append(beta_auto)
    argss.append(SM)
    argss.append(n)
    argss.append(estimate_bool)
    argss.append(crit_bool)
    argss.append(expfile)

    #Initialize PSO method
    #PSO.LJPSO(nparameter,ndata,nswarm,objFunc_Tc_Pc,argss,p)

    
    #Initialize Newton method
    phi = 1.0
    L = 5.6e-10
    L_orig = L
    tol = 1e-5
    k = 0
    err = 1.0
    par = np.empty((2))
    par[0] = L
    par[1] = phi
    Fobj = 1.0
    Lold = 1.0
    while err>tol:

        #Calculate objective function
        F = objFunc_Tc_Pc(par,argss)
        Fold = Fobj
        Fobj = F[0]     #Objective Function value
        argss[3] = F[1]-5 #New temperature estimate is near the last critical temperature
        Tc = F[1]
        Pc = F[2]
        rhoc = F[3]

        #Apply correction to parameter L
        derF = (Fobj-Fold)/(L-Lold)
        Lold = L
        L = L - 0.5*Fobj/derF #damping to help
        if k==0:
            L = L_orig-0.1e-10
        par[0] = L
        errL = L-Lold
        err = Fobj

        k = k+1
        print Tc,Pc,rhoc,L+0.5*Fobj/derF,err
    
    """
    #Initialize Full Newton method
    phi = 5.0
    L = 5e-10
    L_orig = L
    phi_orig = phi
    tol = 1e-5
    k = 0
    err = 1.0
    par = np.empty((2))
    Fv = np.empty((2))
    dF = np.empty((2,2))
    invD = np.empty((2,2))
    D = np.empty((2))
    par[0] = L
    par[1] = phi
    Fobj = 1.0
    Fobj_T = 1.0
    Fobj_P = 1.0
    Lold = 1.0
    phiold = 1.0
    while err>tol:

        #Calculate objective function
        F = objFunc_Tc(par,argss)
        Fold_T = Fobj_T
        Fold_P = Fobj_P
        Fobj = F[0]     #Objective Function value
        Fobj_T = F[4]     #Objective Function value for critical temperature
        Fobj_P = F[5]     #Objective Function value for critical pressure
        argss[3] = F[1]-5 #New temperature estimate is near the last critical temperature
        Tc = F[1]
        Pc = F[2]
        rhoc = F[3]

        #Apply correction to parameter L
        derF_TL = (Fobj_T-Fold_T)/(L-Lold)
        derF_Tphi = (Fobj_T-Fold_T)/(phi-phiold)
        derF_PL = (Fobj_P-Fold_P)/(L-Lold)
        derF_Pphi = (Fobj_P-Fold_P)/(phi-phiold)

        #print L,Lold,phi,phiold
        #print Fobj_T,Fobj_P
        #print derF_TL,derF_Tphi,derF_PL,derF_Pphi

        Lold = L
        phiold = phi
        

        Fv[0] = -Fobj_T
        Fv[1] = -Fobj_P

        dF[0][0] = derF_TL
        dF[0][1] = derF_Tphi
        dF[1][0] = derF_PL
        dF[1][1] = derF_Pphi
        print '---dF---'
        print dF

        detD = derF_TL*derF_Pphi-derF_Tphi*derF_PL
        print '---det---'
        print detD
        invD[0][0] = derF_Pphi/detD
        invD[0][1] = -derF_Tphi/detD
        invD[1][0] = -derF_PL/detD
        invD[1][1] = derF_TL/detD
        D[0] = (invD[0][0]*Fv[0])+(invD[0][1]*Fv[1])
        D[1] = (invD[1][0]*Fv[0])+(invD[1][1]*Fv[1])
        print '---D---'
        print D
        
        #D = inv(dF)*Fv

        #L   = L   + D[0]
        #phi = phi + D[1]

        L   = L   + Fobj_T/derF_TL
        phi = phi + Fobj_P/derF_Pphi

        if k==0:
            L = L_orig+0.1e-10
            phi = phi_orig+0.2

        par[0] = L
        par[1] = phi
        errL = L-Lold
        errphi = phi-phiold
        err = Fobj

        k = k+1
        print Tc,Pc,rhoc,Lold,phiold,err
    """
    par = 1
    return par
#====================================================================================== 

#Calculates critical point for given mixture composition-------------------------------
def crit_mix(estimate_crit,Tci,rhoci,EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,x):

    #Initial estimate of critical point
    if estimate_crit==True:
        Tci = Tci
        Vci = rhoci
    else:
        Tcvec = data.Tc(IDs)
        Tci = 1.5*np.dot(x,Tcvec)
        b = eos.b_calc(IDs,EoS)
        bmix = eos.bmix_calc(MR,b,x)
        Vci = 4*bmix
    Tc = Tci
    P = 0.1
    
    print 'Initial Estimates:','Tc=',Tci,'Vci=',Vci
    raw_input('...')
        
    #Solve det(Q)=0 for initial condition
    nt = 1.0 #Total number of moles   
    delta = np.identity(nc) #Kronecker delta
    
    #Calculate thermodynamic properties at initial T
    nx = 4 #Renorm will calculate only given x
    nd = 400
    
    h = 1e-4
    x0 = np.array([0.1,0.9])
    nt = np.sum(x0)
    n0 = nt*x0
    
    tolC = 1e-20
    tolQ = 1e-20
    C = tolC+1
    detQ = tolQ+1
    Q = np.zeros((nc,nc))
    n1 = np.zeros((nc))
    n1 = n0
    
    r_dat0 = 1
    r_dat1 = 1
    
    
    #External V loop - solve triple sum = 0
    while C>tolC:
        #Internal T loop - solve det(Q)=0------------------------------------------------------------------------------------
        while detQ>tolQ:
            
            dT = 1e-5*Tc
            
            #Q at T
            Q = np.zeros((nc,nc))
            if EoS==2 or EoS==4 or EoS==6:
                r_dat0 = renorm(EoS,IDs,MR,Tc,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)    
            lnfugcoef0 = eos.lnfugcoef_func(IDs,EoS,MR,P,Tc,x0,kij, 1,Vci,en_auto,beta_auto,CR,SM,0,0,r_dat0)[0] #vapor
            for j in range(0,nc):
                n1[j] = n0[j]+h
                x1 = n1/(nt+h)
                if EoS==2 or EoS==4 or EoS==6:
                    r_dat1 = renorm(EoS,IDs,MR,Tc,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
                lnfugcoef1 = eos.lnfugcoef_func(IDs,EoS,MR,P,Tc,x1,kij, 1,Vci,en_auto,beta_auto,CR,SM,0,0,r_dat1)[0] #vapor
                print lnfugcoef1,x1
                for i in range(0,nc):
                    Q[i][j] = (lnfugcoef1[i]-lnfugcoef0[i])/h
                n1 = n0
                x1 = x0
            detQ = np.linalg.det(Q)
            print 'lnfug0',lnfugcoef0,x0
            print 'lnfug1',lnfugcoef1,x1
            print 'detQ=',detQ
            print 'Q=',Q
            raw_input('...')
            
            #Q at T+dT
            Q = np.zeros((nc,nc))
            if EoS==2 or EoS==4 or EoS==6:
                r_dat0 = renorm(EoS,IDs,MR,Tc+dT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
            lnfugcoef0 = eos.lnfugcoef_func(IDs,EoS,MR,P,Tc+dT,x0,kij, 1,Vci,en_auto,beta_auto,CR,SM,0,0,r_dat0)[0] #vapor
            for j in range(0,nc):
                n1[j] = n0[j]+h
                x1 = n1/(nt+h)
                if EoS==2 or EoS==4 or EoS==6:
                    r_dat1 = renorm(EoS,IDs,MR,Tc+dT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
                lnfugcoef1 = eos.lnfugcoef_func(IDs,EoS,MR,P,Tc+dT,x1,kij, 1,Vci,en_auto,beta_auto,CR,SM,0,0,r_dat1)[0] #vapor
                for i in range(0,nc):
                    Q[i][j] = (lnfugcoef1[i]-lnfugcoef0[i])/h
                n1 = n0
            detQ_dT = np.linalg.det(Q)
            
            dF = (detQ_dT-detQ)/(dT*Tc)
            Tc = Tc - detQ/dF
            print 'Tc=',Tc
            #End Internal T loop - solve det(Q)=0===============================================================================
            
        nt = 1.0
        nt0 = 1.0
        n0 = x0*nt0
        dn = 1.0
        """
        n1k = np.zeros((4))
        nt1k = np.zeros((4))
        x1k = np.zeros((4))
        k = -2
        for i in range(0,4)
            n1k[i] = (x0+k*1e-3*dn/nt0)*nt0
            nt1k[i] = np.sum(n1)
            x1k[i] = n1/nt1
            k = k+1
        """
        
        dV = 1e-6*Vci
            
        #C at V
        C = 0
        if EoS==2 or EoS==4 or EoS==6:
            r_dat0 = renorm(EoS,IDs,MR,Tc,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
        lnfugcoef0 = eos.lnfugcoef_func(IDs,EoS,MR,P,Tc,x0,kij, 1,Vci,en_auto,beta_auto,CR,SM,0,0,r_dat0)[0] #vapor
        for j in range(0,nc):
            eps = 1e-3
            n1[j] = (x0[j]+eps*dn/nt)*nt
            nt1 = np.sum(n1)
            x1 = n1/nt1
            if EoS==2 or EoS==4 or EoS==6:
                r_dat1 = renorm(EoS,IDs,MR,Tc,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
            lnfugcoef1 = eos.lnfugcoef_func(IDs,EoS,MR,P,Tc,x1,kij, 1,Vci,en_auto,beta_auto,CR,SM,0,0,r_dat1)[0] #vapor
            for i in range(0,nc):
                Q[i][j] = (lnfugcoef1[i]-lnfugcoef0[i])/eps
                C = C + Q[i][j]*n1[j]*n1[j]
            n1 = n0
        
        #C at T+dV
        C_dV = 0
        if EoS==2 or EoS==4 or EoS==6:
            r_dat0 = renorm(EoS,IDs,MR,Tc,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
        lnfugcoef0 = eos.lnfugcoef_func(IDs,EoS,MR,P,Tc,x0,kij, 1,Vci+dV,en_auto,beta_auto,CR,SM,0,0,r_dat0)[0] #vapor
        for j in range(0,nc):
            eps = 1e-3
            n1[j] = (x0[j]+eps*dn/nt)*nt
            nt1 = np.sum(n1)
            x1 = n1/nt1
            if EoS==2 or EoS==4 or EoS==6:
                r_dat1 = renorm(EoS,IDs,MR,Tc,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
            lnfugcoef1 = eos.lnfugcoef_func(IDs,EoS,MR,P,Tc,x1,kij, 1,Vci+dV,en_auto,beta_auto,CR,SM,0,0,r_dat1)[0] #vapor
            for i in range(0,nc):
                Q[i][j] = (lnfugcoef1[i]-lnfugcoef0[i])/eps
                C_dV = C_dV + Q[i][j]*n1[j]*n1[j]
            n1 = n0
        
        dF = (C_dV-C)/(dV*Vci)
        Vci = Vci - Vci/dF
    
    #Solve det(Q)=0 to find Tc and rhoc
        
    print Tc
    print Vci
    rhoc = 1/Vci
    crit = []
    crit.append(Tc)
    #crit.append(Pc)
    crit.append(rhoc)
    return crit
#======================================================================================
