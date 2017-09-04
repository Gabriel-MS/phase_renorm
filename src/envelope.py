#VLE algorithms to calculate and plot envelopes
import numpy as np
import scipy
import data
import menus
import eos
import math
import correlations
import association
import numerical
import renormalization
from scipy.interpolate import InterpolatedUnivariateSpline, splrep, splev, interp1d

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = 8.314462175e-6 #m3.MPa/K/mol

#Calculates initial dew point to given initial z, P, and T guess-----------------------
def initial_dewT(z,T,P,IDs,EoS,MR,kij):
    z = np.array(z)
    nc = z.shape[0]
    Tc = np.array(data.Tc(IDs))
    Pc = np.array(data.Pc(IDs))
    omega = np.array(data.omega(IDs))

    #Initial point calculation, find dew T and y, given P, beta
    beta = 1 #specified, considering bubble point
    errT = 1.0
    tolT = 1e-3
    iteri = 1
    
    while (errT>tolT):
        K = Pc/P*(np.exp(5.373*(1+omega)*(1-Tc/T))) #Initial estimate for K, Wilson
        K2 = Pc/P*5.373*(1+omega)*Tc*np.exp(5.373*(1+omega)*(1-Tc/T))/T/T
        sumKz = np.dot(K,z)
        sumKz2 = np.sum(K2)
        fT = P*(1-sumKz)
        dfT = P*(-sumKz2)
        if iteri<2:
            dfT = 1e4
        step = -fT/dfT
        T = T + step
        errT = abs(fT/dfT)
        print errT
        input('wow')
        iteri = iteri + 1
    
    K = Pc/P*np.exp(5.373*(1+omega)*(1-Tc/T))
    x = K*z
    y = z
    sumKz = np.dot(K,z)
	
	#Normalize x
    x = x/sumKz
    
    out = []
    out.append(T)
    out.append(x)
    
    return out
#======================================================================================

#Calculates initial dew point to given initial z, P, and T guess-----------------------
def dewT_guess(z,T,P,IDs,EoS,MR,kij):
    z = np.array(z)
    nc = z.shape[0]
    Tc = np.array(data.Tc(IDs))
    Pc = np.array(data.Pc(IDs))
    omega = np.array(data.omega(IDs))

    T_old = T - 1.0
    while T_old!=T:
        T_old = T
        lnfugcoef_v = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,z,kij, 1) #Vapor
        lnfugcoef_l = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,z,kij,-1) #Liquid
        K = np.exp(lnfugcoef_l - lnfugcoef_v)
        Gvap = 0.0
        Gliq = 0.0
        for i in range(0,nc):
            Gvap = Gvap + z[i]*lnfugcoef_v[i]
            Gliq = Gliq + z[i]*lnfugcoef_l[i]
        if (Gliq<Gvap) or (Gliq==Gvap):
            T = T + 10.0
    
    x = K*z
    
    out = []
    out.append(T)
    out.append(x)
    out.append(K)
    return out
#======================================================================================

#Calculates initial dew point to given initial z, x, P, and T guess--------------------
def dewT(z,x,T,P,IDs,EoS,MR,kij):
    z = np.array(z)
    x = np.array(x)
    nc = z.shape[0]
    Tc = np.array(data.Tc(IDs))
    Pc = np.array(data.Pc(IDs))
    omega = np.array(data.omega(IDs))

    #Initial point calculation, find dew T and y, given P, beta
    beta = 1 #specified, considering bubble point
    err = 1.0
    tol = 1e-5
    iteri = 1
    h = 1e-5
    
    while (err>tol or iteri<3):
        #ln fugacity coefficient
        lnfugcoef_v = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,z,kij, 1) #Vapor
        lnfugcoef_l = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,-1) #Liquid
        K = np.exp(lnfugcoef_l-lnfugcoef_v)
        
        #ln fugacity coefficient derivative with respect to Temperature
        lnfugcoef_v2 = eos.lnfugcoef_calc(IDs,EoS,MR,P,T+h,z,kij, 1) #Vapor
        lnfugcoef_l2 = eos.lnfugcoef_calc(IDs,EoS,MR,P,T+h,x,kij,-1) #Liquid
        lnfugcoef_v1 = eos.lnfugcoef_calc(IDs,EoS,MR,P,T-h,z,kij, 1) #Vapor
        lnfugcoef_l1 = eos.lnfugcoef_calc(IDs,EoS,MR,P,T-h,x,kij,-1) #Liquid
        
        F = np.dot(K,z)-1
        dF2 = np.sum(K*z*(lnfugcoef_l2-lnfugcoef_v2))
        dF1 = np.sum(K*z*(lnfugcoef_l1-lnfugcoef_v1))
        dF = (dF2-dF1)/(2*h)
        
        delT = -F/dF
        
        T = T + delT
        
        lnfugcoef_v = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,z,kij, 1) #Vapor
        lnfugcoef_l = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,-1) #Liquid
        K = np.exp(lnfugcoef_l-lnfugcoef_v)
        
        x = K*z
        err = abs(F/dF)
        iteri = iteri + 1
    
    K = np.exp(lnfugcoef_l-lnfugcoef_v)
    x = K*z
    y = z
    sumKz = np.dot(K,z)
	
	#Normalize x
    x = x/sumKz
    
    out = []
    out.append(T)
    out.append(x)
    out.append(K)
    return out
#======================================================================================

#Calculates initial dew point to given z, P, and T guess, using sucessive substitution-
def dewT_successive(z,x,K,T,Torig,P,IDs,EoS,MR,kij):
    tol = 1.0e-5
    z = np.array(z)
    nc = z.shape[0]
    h = 1.0e-6
    step = tol+1
    maxit = 100
    it = 0
    
    Tc = np.array(data.Tc(IDs))
    Pc = np.array(data.Pc(IDs))
    omega = np.array(data.omega(IDs))
    K = Pc/P*(np.exp(5.373*(1+omega)*(1-Tc/T))) #Initial estimate for K, Wilson
    K = 1.0/(np.exp(np.log(Pc/P) + 5.373*(1.0+omega)*(1.0 - Tc/T)))
  
    while (abs(step)>tol and it<maxit):
        it = it + 1
        T_old = T
      
        n = 0
        x = z*K
        n = np.sum(x)
      
        F = -1.0
        dF = 0.0
        x = x/n
        F = F + np.sum(z*K)
    
        T = T_old + h
        lnphi_z = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,z,kij, 1)
        lnphi_x = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,-1)
        dF = dF + np.sum(z*K*(lnphi_z-lnphi_x))
        
        T = T_old - h
        lnphi_z = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,z,kij, 1)
        lnphi_x = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,-1)
        dF = dF - np.sum(z*K*(lnphi_z-lnphi_x))
            
        dF = dF/(2*h)
        step = F/dF
        
        if abs(step>0.25*T_old):
            step = 0.25*T_old*step/abs(step)
        
        T = T_old - step
        
        lnphi_z = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,z,kij, 1)
        lnphi_x = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,-1)
        K = np.exp(lnphi_z-lnphi_x)
    
    x = K*z
    lnphi_z = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,z,kij, 1)
    lnphi_x = eos.lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,-1)
    K = np.exp(lnphi_z-lnphi_x)
    
    out = []
    out.append(T)
    out.append(x)
    out.append(K)
    return out
#======================================================================================

#Given initial dew point, calculate PT envelope using continuation method--------------
def PT_envelope(z,T,P,IDs,EoS,MR,kij,en_auto,beta_auto,CR,SM):
    
    #Reference for this algorithm:
    #Algorithm based on other algorithm developed by Rafael Pereira and Iuri Segtovich
    #Avaible at PyTherm
    #https://github.com/iurisegtovich/PyTherm-applied-thermodynamics/tree/master/contents/models-and-algorithms-laboratory/LVE%20algorithms/PxT_L-V_Phase_Envelope_given_z
    
    #Initial dew point estimation, find T guess
    Toriginal = T
    idewT = dewT_guess(z,T,P,IDs,EoS,MR,kij)
    T = idewT[0]
    x = idewT[1]
    K = idewT[2]

    #First Dew point calculation
    dew = dewT_successive(z,x,K,T,Toriginal,P,IDs,EoS,MR,kij)
    T = dew[0]
    x = dew[1]
    K = dew[2]
    
    #Initial configuration in first step
    y = np.array(z)
    K = np.array(K)
    nc = K.shape[0]
    lnfugcoef_v = eos.lnfugcoef_func(IDs,EoS,MR,P,T,x,kij, 1,0.0,en_auto,beta_auto,CR,SM,0,0)[0] #Vapor
    lnfugcoef_l = eos.lnfugcoef_func(IDs,EoS,MR,P,T,x,kij,-1,0.0,en_auto,beta_auto,CR,SM,0,0)[0] #Liquid
    K = np.exp(lnfugcoef_v-lnfugcoef_l)
    lnK = np.log(K)
    
    #Continuation method specifications--------------------------------
    F = np.empty((nc+2))
    alfa = np.empty((nc+2))
    J = np.empty((nc+2,nc+2))
    X = []
    wlist = []
    rlist = []
    Tlist = []
    Plist = []
    
    #Specified variable and position in vector
    ns = nc + 2 - 1 #-1 because first element in zero index
    S = np.log(P)
    delS = 0.1
    
    #Dependent variables
    lnT = np.log(T)
    for i in range(0,nc):
        X.append(lnK[i])
    X.append(lnT)
    X.append(S) #Specified variable
    X = np.array(X)
    
    #dFdS vector
    dFdS = np.empty((nc+2))
    for i in range(0,nc+1):
        dFdS[i] = 0
    dFdS[ns] = -1
    
    #Stop conditions
    Plow  = 0.025 #Stop condition for low pressure
    Phigh = 30.0 #Stop condition for high pressure
    Tmax_step = 5.0 #Max step for temperature
    lnKmax_step = 0.05 #Max step for lnK
    
    #Critical point detector
    critK = 0.05 #Reference value for approaching critical point
    
    h = 1e-6 #Value to calculate numerical derivatives
    
    tolN = 1e-7 #Tolerance
    maxitN = 50
    
    phase = -1 #Incipient phase is liquid
    w = x      #incipient phase is liquid
    r = y      #Reference phase
    ##################################################################
    
    #Continuation method----------------------------------------------
    print 'Begin Continuation Method'
    critKpass = False
    pt = 0 #Point counter
    pt1 = 0 #
    while P>Plow and P<Phigh: #Pressure between defined limits
        mstep = tolN+1
        itN = 0
        while (mstep>tolN and itN<maxitN): #Newton loop
        
            #Equations
            lnfugcoef_r = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij,-phase,0.0,en_auto,beta_auto,CR,SM,0,0)[0] #Reference
            lnfugcoef_w = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0)[0] #Incipient
            for i in range(0,nc):
                F[i] = (lnK[i]+lnfugcoef_w[i]-lnfugcoef_r[i])
            F[nc+1-1] = (np.sum(w-r))
            F[nc+2-1] = (X[ns]-S)
            #print 'F',F
            
            #Jacobian
            #NC equations derivatives with respect to mole number
            J = PT_jacobian(IDs,EoS,MR,P,T,w,r,h,kij,phase,X,ns,en_auto,beta_auto,CR,SM)
            #print 'J',J
            
            #Solve system
            step = scipy.linalg.solve(J,F)
            #print 'J',J
            #print 'F',F
            #print 'step',step
            if ((phase==1) and (EoS==5) and (itN<20) and pt1<20):
                step = step/2.5
            
            #PROBABLY J IS WRONG AFTER CRITICAL POINT=====================================
            #if phase==1:
            #    for i in range(0,nc+2):
            #        step[i] = -step[i]
            #print 'step',step
            #PROBABLY J IS WRONG AFTER CRITICAL POINT=====================================
            #input('--------------')
            
            #Update independent variables-------------------------------------------
            mstep = 0.0
            #print 'Xbevor',X
            for i in range(0,nc+2):
                X[i] = X[i]-step[i]
                if mstep<abs(step[i]):
                    mstep = abs(step[i])
            #print X
                    
            #Update variables-------------------------------------------------------
            for i in range(0,nc):
                lnK[i] = X[i]
                K[i] = np.exp(lnK[i])
            w = r*K #Adjust incipient phase
            T = np.exp(X[nc+1-1])
            P = np.exp(X[nc+2-1])
            
            if phase==1:
                pt1 = pt1 + 1
            #    print 'X',X
            #    print 'expX',np.exp(X)
            #    print 'T,P',T,P
            #    input ('...')
            
            itN = itN+1
        
        #Calculate Sensivities--------------------------------------------------
        ns_old = ns
        dXdS = scipy.linalg.solve(J,dFdS)
        dXdS = -dXdS
        #print 'dXdS',dXdS
        
        #Search highest sensivity and update specification variable position----
        ns = np.argmax(np.absolute(dXdS))
        #print 'dXdS',dXdS
        #print 'new ns:',ns
        
        #Translating variation in specification variable------------------------
        if ns!=ns_old:
            delS = dXdS[ns]*delS
            for i in range(0,nc+2):
                if i!=ns:
                    dXdS = dXdS/dXdS[ns]
            dXdS[ns] = 1.0
            S = X[ns] #Update specification variable
            #print 'SSS',ns,X[ns],S
        
        #Step in S for next point-----------------------------------------------
        #delmax = max(np.sqrt(abs(X[ns]))/10.0,0.1)
        #delmax = delmax*abs(delS)/delS
        #updel = delS*3/itN
        #if delS>0:
        #    delS = min(updel,delmax)
        #else:
        #    delS = max(updel,-delmax)
        #S = S + delS
        #print 'new S:',S
        
        delmax = ((abs(X[ns]))**0.5)/10.0
        if delmax < 0.1:
            delmax = 0.1
        delmax = abs(delmax)*(abs(delS)/delS)
        delS = delS*4.0/itN
        if abs(delmax) < abs(delS):
            delS = delmax
        S = S + delS
        #print 'SSS',ns,X[ns],S
        
        #Estimates for next point-----------------------------------------------
        Told = T
        lnKold = np.empty((nc))
        for i in range(0,nc):
            lnKold[i] = X[i]
        Xold = X
        X = Xold + dXdS*delS    #STEP
        #print 'SSS',ns,X[ns],S
        for i in range(0,nc):
            lnK[i] = X[i]
        T = np.exp(X[nc+1-1])
        #print 'Xold',Xold
        #print 'X',X
        #print 'dXdS',dXdS
        #print 'delS',delS
        #print 'New estimates:',np.exp(X)
        
        #Adjust large temperature stepsize--------------------------------------
        while abs(T-Told) > Tmax_step:
            delS = delS/2.0
            S = S - delS
            X = X - dXdS*delS
            T = np.exp(X[nc+1-1])
            #print 'SSS',ns,X[ns],S
        #Adjust large lnK stepsize
        #while np.amax(np.absolute(lnK-lnKold)) > lnKmax_step:
        #    delS = delS/2.0
        #    S = S - delS
        #    X = X - dXdS*delS
        #    for i in range(0,nc):
        #        lnK[i] = X[i]
      
        #Adjust large lnK stepsize------------------------------------------------
        lnKm = 0.0
        for i in range(0,nc):
            if abs(X[i])>lnKm:
                lnKm = abs(X[i])
              
        if lnKm<0.1:
            mstep = 0.0
            for i in range(0,nc):
                if abs(delS*dXdS[i])>mstep:
                    mstep = abs(delS*dXdS[i])

            while mstep>critK:
                #print 'mstep>critK'
                delS = delS/2.0
                S = S - delS
                X = X - dXdS*delS

                mstep = 0.0
                for i in range(0,nc):
                    if abs(delS*dXdS[i])>mstep:
                        mstep = abs(delS*dXdS[i])
        
        #Jump critical point----------------------------------------------------
        critKpass = False
        if lnKm<critK:
            while lnKm<critK:
                S = S + delS
                X = X + dXdS*delS
                lnKm = 0.0
                for i in range(0,nc):
                    if X[i]>lnKm:
                        lnKm = X[i]
            phase = -phase #Change incipient phase
            critKpass = True
        #print 'SSS',ns,X[ns],S
        
        #Update variables-------------------------------------------------------
        #if (ns<(nc+1-1)) and (phase==1):
        #    S = -S
        for i in range(0,nc):
            lnK[i] = X[i]
            K[i] = np.exp(lnK[i])
            #if phase==1:
            #    K[i] = 1/K[i]
            #    lnK[i] = np.log(K[i])
            #    X[i] = lnK[i]
        w = r*K #Adjust incipient phase
        T = np.exp(X[nc+1-1])
        P = np.exp(X[nc+2-1])
        #if phase==1:
        #    print 'X',X
        #    print 'expX',np.exp(X)
        #    print 'T,P',T,P
        #    input ('...at end')
        
        print 'Phase:',phase,'lnK:',lnK[0],lnK[1],'T:',T,'P:',P,'w1:',w[0],'r1:',r[0]
        wlist.append(w[0])
        rlist.append(r[0])
        Tlist.append(T)
        Plist.append(P)
        itN = 0
        
        pt = pt+1
                        
    out = []
    out.append(wlist)
    out.append(rlist)
    out.append(Tlist)
    out.append(Plist)
    return out
#======================================================================================

#Build Jacobian to be used in continuation method--------------------------------------
def PT_jacobian(IDs,EoS,MR,P,T,w,r,h,kij,phase,X,ns,en_auto,beta_auto,CR,SM):
    nc = w.shape[0]
    J = np.empty((nc+2,nc+2))
    
    #NC equations derivatives with respect to lnK
    for i in range(0,nc):
        for j in range(0,nc):
            w_orig = w[j]
            w[j] = w_orig + w_orig*h
            lnphi_2 = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0)[0]
            #print 'lnphi2',lnphi_2[i],lnphi_2[i],w[j],lnphi_2
            w[j] = w_orig - w_orig*h
            lnphi_1 = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0)[0]
            #print 'lnphi2-lnphi1',lnphi_2[i]-lnphi_1[i],lnphi_1[i],w[j]
            w[j] = w_orig
            J[i][j] = (lnphi_2[i]-lnphi_1[i])*w[j]/(2*w[j]*h)
            #print 'J',J[i][j],w[j],w[j]/(2*w[j]*h)
            #print 'P',P,'T',T,'kij',kij,'comp',w
            if i==j:
                J[i][j] = J[i][j]+1
    #print 'J ln K',J
    #input('-----------')
    #print '----------'
                    
    #NC equations derivatives with respect to lnT
    for i in range(0,nc):
        T_orig = X[nc+1-1]
        T = np.exp(T_orig*(1+h))
        lnphi_2w = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0)[0]
        lnphi_2r = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij,-phase,0.0,en_auto,beta_auto,CR,SM,0,0)[0]
        T = np.exp(T_orig*(1-h))
        lnphi_1w = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0)[0]
        lnphi_1r = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij,-phase,0.0,en_auto,beta_auto,CR,SM,0,0)[0]
            
        J[i][nc+1-1] = ((lnphi_2w[i]-lnphi_2r[i])-(lnphi_1w[i]-lnphi_1r[i]))/(2*T_orig*h)
        T = np.exp(X[nc+1-1])
            
    #NC equations derivatives with respect to lnP
    for i in range(0,nc):
        P_orig = X[nc+2-1]
        P = np.exp(P_orig*(1+h))
        lnphi_2w = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0)[0]
        lnphi_2r = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij,-phase,0.0,en_auto,beta_auto,CR,SM,0,0)[0]
        P = np.exp(P_orig*(1-h))
        lnphi_1w = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0)[0]
        lnphi_1r = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij,-phase,0.0,en_auto,beta_auto,CR,SM,0,0)[0]
                
        J[i][nc+2-1] = ((lnphi_2w[i]-lnphi_2r[i])-(lnphi_1w[i]-lnphi_1r[i]))/(2*P_orig*h)
        P = np.exp(X[nc+2-1])
                
    #NC+1 equation derivatives
    for i in range(0,nc):
        J[nc+1-1][i] = w[i] #lnKi derivative
    J[nc+1-1][nc+1-1] = 0.0 #lnT derivative
    J[nc+1-1][nc+2-1] = 0.0 #lnP derivative
            
    #NC+2 equation derivatives
    for i in range(0,nc+2):
        J[nc+2-1][i] = 0
    J[nc+2-1][ns] = 1
            
    return J
#======================================================================================

#Report PTenvelope-------------------------------------------------------------------
def report_PT(data,options,title,print_options):
    n = len(data[0])
    p1 = np.array(menus.flatten(data[0]))
    p2 = np.array(menus.flatten(data[1]))
    T = np.array(menus.flatten(data[2]))
    P = np.array(menus.flatten(data[3]))
    header = 'Phase 1;Phase 2;T(K);P(MPa)\n'
    savedir = str('../output/%s' %title)
    with open(savedir,'w') as file:
        file.write(' ')
        file.write('Defined Configuration:----------------------\n')
        file.write('Number of components: %i \n' %print_options[0])
        file.write('Equation of State:    %s \n' %print_options[2])
        file.write('Mixing Rule:          %s \n' %print_options[3])
        file.write('Components:           %s \n' %', '.join(print_options[1]))
        file.write('Feed Composition:     %s \n' %print_options[4])
        file.write('Association Rule:     %s \n' %', '.join(print_options[5]))
        file.write('Combining Rule:       %s \n' %print_options[6])
        file.write('============================================\n')
        file.write(header)
        file.write('\n') 
        for i in range(0,n):
                lin1 = [str(round(p1[i],9)),str(round(p2[i],9)),str(round(T[i],9)),str(round(P[i],9))]
                lin = ('%s\n' % ';'.join(lin1))
                file.write(lin)
#======================================================================================

#Plot PT envelope---------------------------------------------------------------------
def plot_PT(title,xtitle,ytitle,figname,boxtext,xdata,ydata):
    
    xdata = menus.flatten(xdata)
    ydata = menus.flatten(ydata)
    
    xdata = np.array(xdata)
    ydata = np.array(ydata)
    size = xdata.shape[0]
    
    #Create figure
    fig, ax = plt.subplots(1,1)
    
    #Data
    plt.plot(xdata,ydata)
    
    #Labels
    plt.xlabel(xtitle)      #x axis label
    plt.ylabel(ytitle)      #y axis label
    plt.title(title)        #title
    
    #plt.axis([x1data[0],x1data[size-1],ydata[0],ydata[size-1]])
    
    #Text box
    ax.text(0.05, 0.90, boxtext, transform=ax.transAxes, style='italic',
        bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    
    #Save figure with 1000 dpi resolution
    figsavedir = str('../output/%s' %figname)
    fig.savefig(figsavedir, dpi=1000)
#======================================================================================

#Given T, iterate over specified x to find y and P-------------------------------------
def Pxy_envelope(T,IDs,EoS,MR,kij,nc,AR,CR,SM):
    
    Pvec = []
    xvec = []
    yvec = []
    
    #Definitions--------------------------------------------------
    #Main iteration conditions
    x = np.array([0.001,0.999]) #x array
    xf = 1.0                    #Main stop condition
    stepx = 5e-3                #Main step
    it = 0                      #Iteration counter
    ity = 0                     #Iteration over y loop counter
    pt = 0                      #Point counter
    
    #tolerances
    tolK = 1e-5             #iteration over Kx
    toly = 1e-6             #iteration over y
    
    #CPA auto-association configurations
    auto = []
    auto = association.CPA_auto(AR,nc,IDs)
    en_auto = auto[0]
    beta_auto = auto[1]
    #=============================================================
    
    
    #Initial guess------------------------------------------------
    #Saturation pressure
    Psat = np.array(correlations.Psat_antoine(IDs,T))
    
    #Initial guess
    P = np.dot(Psat,x)
    y = Psat*x/P
    Vv = R*T/P #Not going to be used, just starting
    Vl = 0.99  #Not going to be used, just starting
    #=============================================================
    
    #Main iteration start, range x1------------------------------------------
    while x[0]<xf:
    
        #Iteration Kx start--------------------------------------------------
        errK = tolK+1 #Force enter iteration
        while errK>tolK:
            func_l = eos.lnfugcoef_func(IDs,EoS,MR,P,T,x,kij,-1,Vl,en_auto,beta_auto,CR,SM,it,pt) #Liquid
            lnfugcoef_l = func_l[0]
            Vl = func_l[1]
            func_v = eos.lnfugcoef_func(IDs,EoS,MR,P,T,y,kij, 1,Vv,en_auto,beta_auto,CR,SM,it,pt) #Vapor
            lnfugcoef_v = func_v[0]
            Vv = func_v[1]
            K = np.exp(lnfugcoef_l-lnfugcoef_v)
            Kx = K*x
            sumKx = np.sum(Kx)
            
            #print 'phil = ',np.exp(lnfugcoef_l),Vl
            #print 'phiv = ',np.exp(lnfugcoef_v),Vv
            #print 'sumKx = ',sumKx
            
            #Iteration y start-----------------------------------------------
            erry = toly+1
            ity = 0
            while erry>toly or ity<2:
                y = Kx/sumKx
                sumKxold = sumKx
                func_v = eos.lnfugcoef_func(IDs,EoS,MR,P,T,y,kij, 1,Vv,en_auto,beta_auto,CR,SM,it,pt) #Vapor
                lnfugcoef_v = func_v[0]
                Vv = func_v[1]
                K = np.exp(lnfugcoef_l-lnfugcoef_v)
                Kx = K*x
                sumKx = np.sum(Kx)
                erry = abs(sumKx-sumKxold)/sumKx
                ity = ity+1
                #print 'erry',erry
                #print 'y inside erry',y
            #Iteration y end=================================================
            
            errK = abs(sumKx-1)/sumKx
            y = Kx/sumKx
            P = P*sumKx #Applying break condition
            it = it+1
            #print 'y',y
            #print 'P',P,sumKx
            #print 'it',it
            #print 'errK---------------',errK,tolK
            #input('...')
        #Iteration Kx end====================================================
        
        #Save P, x1 and y1
        Pvec.append(P)
        xvec.append(x[0])
        yvec.append(y[0])
        
        #Update main condition
        #print'======='
        #print 'x',x
        #print 'P',P
        #print 'it',it
        #print 'y',y
        #input('End a point--------------')
        print x[0],y[0],P,it
        x[0] = x[0]+stepx
        if stepx==5e-3 and pt==0:
            x[0] = 5e-3
        x[1] = 1-x[0]
        pt = pt+1       #add counter pt
        it = 0          #zero counter it
        
            
    #Main iteration end, range x1============================================
    
    out = []
    out.append(Pvec)
    out.append(xvec)
    out.append(yvec)
    
    return out
#======================================================================================

#Plot Pxy envelope---------------------------------------------------------------------
def plot_Pxy(title,xtitle,ytitle,figname,boxtext,ydata,x1data,x2data):
    
    x1data = menus.flatten(x1data)
    x2data = menus.flatten(x2data)
    ydata = menus.flatten(ydata)
    
    x1data = np.array(x1data)
    x2data = np.array(x2data)
    ydata = np.array(ydata)
    size = x1data.shape[0]
    
    #Create figure
    fig, ax = plt.subplots(1,1)
    
    #Data
    plt.plot(x1data,ydata)
    plt.plot(x2data,ydata)
    
    #Labels
    plt.xlabel(xtitle)      #x axis label
    plt.ylabel(ytitle)      #y axis label
    plt.title(title)        #title
    
    #plt.axis([x1data[0],x1data[size-1],ydata[0],ydata[size-1]])
    
    #Text box
    ax.text(0.05, 0.90, boxtext, transform=ax.transAxes, style='italic',
        bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    
    #Save figure with 1000 dpi resolution
    figsavedir = str('../output/%s' %figname)
    fig.savefig(figsavedir, dpi=1000)
#======================================================================================

#Report Pxy envelope-------------------------------------------------------------------
def report_Pxy(data,options,title,print_options):
    n = len(data[0])
    P = np.array(menus.flatten(data[0]))
    x = np.array(menus.flatten(data[1]))
    y = np.array(menus.flatten(data[2]))
    header = 'P(MPa);x1;y1\n'
    savedir = str('../output/%s' %title)
    with open(savedir,'w') as file:
        file.write(' ')
        file.write('Defined Configuration:----------------------\n')
        file.write('Number of components: %i \n' %print_options[0])
        file.write('Equation of State:    %s \n' %print_options[2])
        file.write('Mixing Rule:          %s \n' %print_options[3])
        file.write('Components:           %s \n' %', '.join(print_options[1]))
        file.write('Feed Composition:     %s \n' %print_options[4])
        file.write('Association Rule:     %s \n' %', '.join(print_options[5]))
        file.write('Combining Rule:       %s \n' %print_options[6])
        file.write('============================================\n')
        file.write(header)
        file.write('\n')
        for i in range(0,n):
                lin1 = [str(round(P[i],9)),str(round(x[i],9)),str(round(y[i],9))]
                lin = ('%s\n' % ';'.join(lin1))
                file.write(lin)
#======================================================================================

#Given pure component isotherm, calculates phase coexistence densities-----------------
def coexistence_dens(rho1,f1):
    
    rho1 = np.array(rho1)
    f1 = np.array(f1)
    n = f1.shape[1]
    rho = rho1.flatten()
    f = f1.flatten()
    #Spline to get chemical potential
    fspl = splrep(rho,f,k=3)         #Cubic Spline Representation
    u = splev(rho,fspl,der=1)        #Evaluate Cubic Spline First derivative
     
    #Calcule pressure
    P = -f+rho*u
    
    #First derivative
    drho = rho[1]-rho[0]
    dPdrho = np.diff(P)/drho
    
    #Find max and min pressure of isotherm inside binodal curve
    max1 = int(numerical.bin_max(dPdrho))
    min1 = int(numerical.bin_min(dPdrho))
    rhomax = rho[max1]
    rhomin = rho[min1]
    Pmax = P[max1]
    Pmin = P[min1]
    min2 = min1+30
    max2 = max1-30

    if Pmin>Pmax:
        min1 = numerical.bin_min_seed(dPdrho,max1)
        rhomin = rho[min1]
        Pmin = P[min1]
        min2 = min1+10

    if Pmin<0:
        Pmin=1e-3
    
    Pf1 = np.empty((n))
    Pf2 = np.empty((n))
    i = 0
    for i in range (0,n):
        if Pmin<0:
            Pf1[i] = P[i] + Pmin
        else:
            Pf1[i] = P[i] - Pmin
        Pf2[i] = P[i] - Pmax
    Pf1 = np.array(Pf1)
    Pf2 = np.array(Pf2)

    while Pf2[min1]*Pf2[min2]>0:
        min2 = min2 + 10
    
    if max2<0:
        max2 = 0
    else:
        while Pf1[max1]*Pf1[max2]>0 and max2>=0:
            max2 = max2 - 1

    #Find initial guess for densities
    Pf1spln = splrep(rho,Pf1,k=3)
    Pf2spln = splrep(rho,Pf2,k=3)
    Pf1rr = splev(rho,Pf1spln)
    Pf2rr = splev(rho,Pf2spln)
    Pf1roots = InterpolatedUnivariateSpline(rho,Pf1rr).roots()
    Pf2roots = InterpolatedUnivariateSpline(rho,Pf2rr).roots()
    
    i = 0
    rho1=rho[0]-1
    while rho1<rho[0] and rho1<0:
        rho1 = Pf1roots[i]
        i = i+1
    if Pf1roots[0]>Pf2roots[0]:
        rho1 = 0.1
    rho1 = numerical.falsi_spline(rho,Pf1,rho[0],rhomax,1e-3)
    
    
    i = 0
    rho2 = rhomin-1
    while rho2<rhomin:
        rho2 = Pf2roots[i]
        i = i+1
    rho2 = numerical.falsi_spline(rho,Pf2,rhomin,rho[min2],1e-3)

    #Solve newton-raphson system
    tol = 1e-6
    drho1 = tol+1
    drho2 = tol+1
    du = tol+1
    dP = tol+1
    drho2old = tol+1
    drho1old = tol+1
    stop = 0.1
    counter = 0
    
    uspl = splrep(rho,u,k=3)
    Pspl = splrep(rho,P,k=3)
    fspl1 = splev(rho,fspl)
    uspl1 = splev(rho,uspl)
    dudrho = splev(rho,uspl,der=1)
    dPdrho = splev(rho,Pspl,der=1)
    
    while (abs(du)>tol or abs(dP)>tol) and (abs(Pmax-Pmin)>1e-3):
        f1 = InterpolatedUnivariateSpline(rho,fspl1,k=3)(rho1)
        f2 = InterpolatedUnivariateSpline(rho,fspl1,k=3)(rho2)
        u1 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho1)
        u2 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho2)
        P1 = -f1+rho1*u1
        P2 = -f2+rho2*u2

        du1 = InterpolatedUnivariateSpline(rho,dudrho,k=3)(rho1)
        du2 = InterpolatedUnivariateSpline(rho,dudrho,k=3)(rho2)
        dP1 = InterpolatedUnivariateSpline(rho,dPdrho,k=3)(rho1)
        dP2 = InterpolatedUnivariateSpline(rho,dPdrho,k=3)(rho2)
        detJ = -dP2*du1+dP1*du2

        drho2 = -du1/detJ*(P1-P2)+dP1/detJ*(u1-u2)
        drho1 = -du2/detJ*(P1-P2)+dP2/detJ*(u1-u2)
        
        rho1 = rho1 + stop*drho1
        rho2 = rho2 + stop*drho2
        
        if counter>0 and (drho1>drho1old and drho2>drho2old):
            rho1 = rho1 - stop*drho1
            rho2 = rho2 - stop*drho2
            stop = stop/1.05 #Break
            rho1 = rho1 + stop*drho1/2
            rho2 = rho2 + stop*drho2/2
            #print stop,counter
        
        du = abs(u1-u2)
        dP = abs(P1-P2)
        print rho1,rho2,du,dP
        
        counter = counter+1
        drho1old = drho1
        drho2old = drho2

    if abs(Pmax-Pmin)<1e-3:
        rho1 = (rhomax+rhomin)/2
        rho2 = rho1

    f1 = InterpolatedUnivariateSpline(rho,fspl1,k=3)(rho1)
    f2 = InterpolatedUnivariateSpline(rho,fspl1,k=3)(rho2)
    u1 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho1)
    u2 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho2)
    P1 = -f1+rho1*u1
    P2 = -f2+rho2*u2

    dens = []
    dens.append(rho1)
    dens.append(rho2)
    dens.append(P1)
    dens.append(P2)
    dens.append(u1)
    dens.append(u2)
    
    return dens
#======================================================================================

#Given initial T, using renormalization method, build P-rho envelope-------------------
def PV_envelope(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n):
    
    env = []
    Tv = []
    rhov = []
    rhol = []
    Pv = []
    
    print 'T:   rhov:   rhol:   P:'
    while T<=Tfinal:
        ren = renormalization.renorm(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n)
        dens = coexistence_dens(ren[2],ren[0])
        Tv.append(T)
        rhov.append(dens[0])
        rhol.append(dens[1])
        Pv.append(dens[2])
        print T,dens[0],dens[1],dens[2]
        T = T + stepT
        
    env.append(Tv)
    env.append(rhov)
    env.append(rhol)
    env.append(Pv)
    return env
#======================================================================================  

#Report PV envelope--------------------------------------------------------------------
def report_PV(data,options,title,print_options):
    n = len(data[0])
    T = np.array(menus.flatten(data[0]))
    rhov = np.array(menus.flatten(data[1]))
    rhol = np.array(menus.flatten(data[2]))
    P = np.array(menus.flatten(data[3]))
    header = 'T;rhov;rhol;P\n'
    savedir = str('../output/%s' %title)
    with open(savedir,'w') as file:
        file.write(' ')
        file.write('Defined Configuration:----------------------\n')
        file.write('Number of components: %i \n' %print_options[0])
        file.write('Equation of State:    %s \n' %print_options[2])
        file.write('Mixing Rule:          %s \n' %print_options[3])
        file.write('Components:           %s \n' %', '.join(print_options[1]))
        file.write('Feed Composition:     %s \n' %print_options[4])
        file.write('Association Rule:     %s \n' %', '.join(print_options[5]))
        file.write('Combining Rule:       %s \n' %print_options[6])
        file.write('============================================\n')
        file.write(header)
        file.write('\n')
        for i in range(0,n):
                lin1 = [str(round(T[i],9)),str(round(rhov[i],9)),str(round(rhol[i],9)),str(round(P[i],9))]
                lin = ('%s\n' % ';'.join(lin1))
                file.write(lin)
#======================================================================================

#Plot PV envelope----------------------------------------------------------------------
def plot_PV(title,xtitle,ytitle,figname,ydata,x1data,x2data):
    
    x1data = menus.flatten(x1data)
    x2data = menus.flatten(x2data)
    ydata = menus.flatten(ydata)
    
    x1data = np.array(x1data)
    x2data = np.array(x2data)
    ydata = np.array(ydata)
    size = x1data.shape[0]
    
    #Create figure
    fig, ax = plt.subplots(1,1)
    
    #Data
    plt.plot(x1data,ydata)
    plt.plot(x2data,ydata)
    
    #Labels
    plt.xlabel(xtitle)      #x axis label
    plt.ylabel(ytitle)      #y axis label
    plt.title(title)        #title
    
    #Axis size
    plt.ylim([ydata[0]-20,ydata[size-1]+20])
    
    #Save figure with 1000 dpi resolution
    figsavedir = str('../output/%s' %figname)
    fig.savefig(figsavedir, dpi=1000)
#======================================================================================

#Define and handle which envelope to calculate-----------------------------------------
def calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type):
    
    if env_type==1:
        #Calculate PT envelope with continuation method*********************************
        print '\nCalculating PT envelope'
        env_PT = PT_envelope(z,T,P,IDs,EoS,MR,kij,en_auto,beta_auto,CR,SM)
        print 'PT envelope calculated'

        print 'Creating PT report'
        reportname = str('PT_%s.csv' %('_'.join(print_options[1])))
        reportname = str('PT_%s.csv' %('_'.join(print_options[1])))
        report_PT(env_PT,print_options,reportname,print_options)
        print ('Report %s saved successfully' %reportname)

        print 'Starting to plot PT envelope'
        title = str('PT envelope\n%s' %(' + '.join(print_options[1])))
        figname = str('PT_%s.png' %('_'.join(print_options[1])))
        zlist = " ,".join(format(x, ".2f") for x in z)
        boxtext = str('z= [%s]' %zlist)
        #boxtext = 'box'
        plot_PT(title,'T(K)','P(MPa)',figname,boxtext,env_PT[2],env_PT[3])
        print ('Figure %s saved successfully' %figname)
        #*******************************************************************************

    if env_type==2:
        #Calculate pure PV envelope*****************************************************
        nd = 400
        nx = 200
        n = 6
        finalT = 520.0
        stepT = 2.5
        print '\nCalculating PV envelope'
        env_PV = PV_envelope(EoS,IDs,MR,T,finalT,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n)
        print 'PV envelope calculated'

        print 'Creating pure PV report'
        reportname = str('PV_%s.csv' %('_'.join(print_options[1])))
        report_PV(env_PV,print_options,reportname,print_options)
        print ('Report %s saved successfully' %reportname)

        print 'Starting to plot PV envelope'
        title = str('PV envelope\n%s' %(' + '.join(print_options[1])))
        figname = str('PV_%s.png' %('_'.join(print_options[1])))
        plot_PV(title,'Density (mol/m3)','T (K)',figname,env_PV[0],env_PV[1],env_PV[2])
        print ('Figure %s saved successfully' %figname)
        #*******************************************************************************

    if env_type==3:
        #Calculate Pxy envelope*********************************************************
        print '\nCalculating Pxy envelope'
        env_pxy = Pxy_envelope(T,IDs,EoS,MR,kij,nc,AR,CR,SM)
        print 'Pxy envelope calculated'

        print 'Creating Pxy report'
        reportname = str('Pxy_%s_%.2fK.csv' %('_'.join(print_options[1]),T))
        report_Pxy(env_pxy,print_options,reportname,print_options)
        print ('Report %s saved successfully' %reportname)

        print 'Starting to plot Pxy envelope'
        title = str('Pxy envelope\n%s' %(' + '.join(print_options[1])))
        figname = str('Pxy_%s_%.2fK.png' %('_'.join(print_options[1]),T))
        boxtext = str('T=%.2fK' %T)
        plot_Pxy(title,'x1,y1','P(MPa)',figname,boxtext,env_pxy[0],env_pxy[1],env_pxy[2])
        print ('Figure %s saved successfully' %figname)
        #*******************************************************************************
#======================================================================================

