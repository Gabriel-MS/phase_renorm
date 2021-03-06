#VLE algorithms to calculate and plot envelopes
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy
import data
import timeit
import menus
import eos
import math
import correlations
import association
import numerical
import renormalization
import critical
import derivativeprop
import PSO
import time
import estimation
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline, splrep, splev, interp1d
from scipy import optimize

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

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
    lnfugcoef_v = eos.lnfugcoef_func(IDs,EoS,MR,P,T,x,kij, 1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Vapor
    lnfugcoef_l = eos.lnfugcoef_func(IDs,EoS,MR,P,T,x,kij,-1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Liquid
    K = np.exp(lnfugcoef_v-lnfugcoef_l)
    lnK = np.log(K)
    
    #Continuation method specifications--------------------------------
    F = np.empty((nc+2))
    alfa = np.empty((nc+2))
    J = np.empty((nc+2,nc+2))
    X = []
    wlist = []
    rlist = []
    w2list = []
    r2list = []
    Tlist = []
    Plist = []
    tpd_LLV_list = []
    phase_list = []
    
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
    Plow  = 0.001 #Stop condition for low pressure
    Phigh = 100.0 #Stop condition for high pressure
    Tmax_step = 0.5 #Max step for temperature    #5.0
    lnKmax_step = 0.03 #Max step for lnK    #0.05
    
    #Critical point detector
    critK = 0.03 #Reference value for approaching critical point    #0.05
    
    h = 1e-5 #Value to calculate numerical derivatives    #1e-6
    
    tolN = 1e-6 #Tolerance    #1e-7
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
            lnfugcoef_r = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij,-phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Reference
            lnfugcoef_w = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Incipient
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

        #Check tpd for LL equilibria--------------------------------------------------------------
        nphases = 2
        tpd_LL = 0.0
        lnfugcoef_wL = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij,-1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Trial Liquid Phase
        lnfugcoef_rL = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij,-1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Trial Liquid Phase
        for jj in range(0,nphases-1): #Number of incipient phases
            for ii in range(0,nc): 
                tpd_LL += w[ii]*(np.log(w[ii])+lnfugcoef_w[ii]-(np.log(r[ii])+lnfugcoef_rL[ii]))

        #Check tpd for LLV equilibria--------------------------------------------------------------
        nphases = 2
        tpd_LLV = 0.0
        #Trial L composition
        k = np.ones((nc))
        k = k/nc*1.0
        lnfugcoef_k = eos.lnfugcoef_func(IDs,EoS,MR,P,T,k,kij,-1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Trial Liquid Phase
        lnfugcoef_w = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij,-1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Trial Liquid Phase
        lnfugcoef_r = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij, 1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Trial Vapor Phase
        for jj in range(0,nphases-1): #Number of incipient phases
            for ii in range(0,nc): 
                tpd_LLV += w[ii]*(np.log(w[ii])+lnfugcoef_w[ii]-(np.log(r[ii])+lnfugcoef_r[ii]))
                tpd_LLV += k[ii]*(np.log(k[ii])+lnfugcoef_k[ii]-(np.log(r[ii])+lnfugcoef_r[ii]))
        #if tpd_LLV<0:
        #   a = PT3_envelope(r,w,k,T,P,IDs,EoS,MR,kij,en_auto,beta_auto,CR,SM)

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
        
        print 'Phase:',phase,'lnK:',lnK[0],lnK[1],'T:',T,'P:',P,'w1:',w[0],'r1:',r[0],'tpd LL:',tpd_LL,'tpd LLV',tpd_LLV
        w = w/np.sum(w)
        wlist.append(w[0])
        rlist.append(r[0])
        w2list.append(w[1])
        r2list.append(r[1])
        Tlist.append(T)
        Plist.append(P)
        tpd_LLV_list.append(tpd_LLV)
        phase_list.append(phase)
        itN = 0
        
        pt = pt+1
                        
    out = []
    out.append(wlist)
    out.append(rlist)
    out.append(Tlist)
    out.append(Plist)
    out.append(tpd_LLV_list)
    out.append(phase_list)
    out.append(w2list)
    out.append(r2list)
    return out
#======================================================================================

#Given initial three-phase equilibria point, calculate PT envelope three-phase lines using continuation method--------------
def PT3_envelope(r,w,z,T,P,IDs,EoS,MR,kij,en_auto,beta_auto,CR,SM):
    print "ENTERED THREE-PHASE EQUILIBRIA"
    
    #Initial configuration in first step
    y = np.array(z)
    nc = r.shape[0]
    lnfugcoefr_v = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij, 1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Vapor
    lnfugcoefw_l = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij,-1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Liquid
    lnfugcoefz_l = eos.lnfugcoef_func(IDs,EoS,MR,P,T,z,kij,-1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Liquid
    K = np.exp(abs(lnfugcoefr_v)-abs(lnfugcoefw_l)-abs(lnfugcoefz_l))
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
    Plow  = 0.001 #Stop condition for low pressure
    Phigh = 100.0 #Stop condition for high pressure
    Tmax_step = 2.0 #Max step for temperature    #5.0
    lnKmax_step = 0.03 #Max step for lnK    #0.05
    
    #Critical point detector
    critK = 0.03 #Reference value for approaching critical point    #0.05
    
    h = 1e-5 #Value to calculate numerical derivatives    #1e-6
    
    tolN = 1e-6 #Tolerance    #1e-7
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
            lnfugcoef_r = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij,-phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Reference
            lnfugcoef_w = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Incipient
            lnfugcoef_z = eos.lnfugcoef_func(IDs,EoS,MR,P,T,z,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Incipient
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

        #Check tpd for LL equilibria--------------------------------------------------------------
        nphases = 2
        tpd_LL = 0.0
        lnfugcoef_wL = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij,-1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Trial Liquid Phase
        lnfugcoef_rL = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij,-1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Trial Liquid Phase
        for jj in range(0,nphases-1): #Number of incipient phases
            for ii in range(0,nc): 
                tpd_LL += w[ii]*(np.log(w[ii])+lnfugcoef_w[ii]-(np.log(r[ii])+lnfugcoef_rL[ii]))

        #Check tpd for LLV equilibria--------------------------------------------------------------
        nphases = 2
        tpd_LLV = 0.0
        #Trial L composition
        k = np.ones((nc))
        k = k/nc*1.0
        lnfugcoef_k = eos.lnfugcoef_func(IDs,EoS,MR,P,T,k,kij,-1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Trial Liquid Phase
        lnfugcoef_w = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij,-1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Trial Liquid Phase
        lnfugcoef_r = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij, 1,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0] #Trial Vapor Phase
        for jj in range(0,nphases-1): #Number of incipient phases
            for ii in range(0,nc): 
                tpd_LLV += w[ii]*(np.log(w[ii])+lnfugcoef_w[ii]-(np.log(r[ii])+lnfugcoef_r[ii]))
                tpd_LLV += k[ii]*(np.log(k[ii])+lnfugcoef_k[ii]-(np.log(r[ii])+lnfugcoef_r[ii]))

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
        
        print '333 ','Phase:',phase,'lnK:',lnK[0],lnK[1],'T:',T,'P:',P,'w1:',w[0],'r1:',r[0],'tpd LL:',tpd_LL,'tpd LLV',tpd_LLV
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
            lnphi_2 = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0]
            #print 'lnphi2',lnphi_2[i],lnphi_2[i],w[j],lnphi_2
            w[j] = w_orig - w_orig*h
            lnphi_1 = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0]
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
        lnphi_2w = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0]
        lnphi_2r = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij,-phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0]
        T = np.exp(T_orig*(1-h))
        lnphi_1w = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0]
        lnphi_1r = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij,-phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0]
            
        J[i][nc+1-1] = ((lnphi_2w[i]-lnphi_2r[i])-(lnphi_1w[i]-lnphi_1r[i]))/(2*T_orig*h)
        T = np.exp(X[nc+1-1])
            
    #NC equations derivatives with respect to lnP
    for i in range(0,nc):
        P_orig = X[nc+2-1]
        P = np.exp(P_orig*(1+h))
        lnphi_2w = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0]
        lnphi_2r = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij,-phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0]
        P = np.exp(P_orig*(1-h))
        lnphi_1w = eos.lnfugcoef_func(IDs,EoS,MR,P,T,w,kij, phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0]
        lnphi_1r = eos.lnfugcoef_func(IDs,EoS,MR,P,T,r,kij,-phase,0.0,en_auto,beta_auto,CR,SM,0,0,0)[0]
                
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
    TPD3 = np.array(menus.flatten(data[4]))
    phase = np.array(menus.flatten(data[5]))
    w2 = np.array(menus.flatten(data[6]))
    r2 = np.array(menus.flatten(data[7]))
    header = 'Phase 1;Phase 2;T(K);P(MPa);TPD-3-phase;Phase_ID;inc2;ref2;\n'
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
                lin1 = [str(round(p1[i],9)),str(round(p2[i],9)),str(round(T[i],9)),str(round(P[i],9)),str(round(TPD3[i],9)),str(round(phase[i],9)),str(round(w2[i],9)),str(round(r2[i],9))]
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
def Pxy_envelope(T,IDs,EoS,MR,kij,nc,AR,CR,SM,r_data):
    
    Pvec = []
    xvec = []
    yvec = []
    
    #Definitions--------------------------------------------------
    #Main iteration conditions
    x = np.array([0.005,0.995]) #x array
    xf = 0.995                    #Main stop condition
    stepx = 5e-4                #Main step
    it = 0                      #Iteration counter
    ity = 0                     #Iteration over y loop counter
    pt = 0                      #Point counter
    itmax = 100
    
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
    print Psat,x,P,y
    Vv = R*T/P #Not going to be used, just starting
    Vl = 0.99  #Not going to be used, just starting
    #P = P/100
    print 'P=',P
    #=============================================================
    
    #Main iteration start, range x1------------------------------------------
    while x[0]<xf:
    
        #Iteration Kx start--------------------------------------------------
        errK = tolK+1 #Force enter iteration
        it = 0
        while errK>tolK and it<itmax:
            func_l = eos.lnfugcoef_func(IDs,EoS,MR,P,T,x,kij,-1,Vl,en_auto,beta_auto,CR,SM,it,pt,r_data) #Liquid
            lnfugcoef_l = func_l[0]
            Vl = func_l[1]
            func_v = eos.lnfugcoef_func(IDs,EoS,MR,P,T,y,kij, 1,Vv,en_auto,beta_auto,CR,SM,it,pt,r_data) #Vapor
            lnfugcoef_v = func_v[0]
            Vv = func_v[1]
            K = np.exp(lnfugcoef_l-lnfugcoef_v)
            Kx = K*x
            sumKx = np.sum(Kx)
            
            #print 'phil = ',np.exp(lnfugcoef_l),Vl
            #print 'phiv = ',np.exp(lnfugcoef_v),Vv
            #print 'sumKx = ',sumKx,K,Kx
            
            #Iteration y start-----------------------------------------------
            erry = toly+1
            ity = 0
            while (erry>toly and ity<100) or ity<2:
                y = Kx/sumKx
                sumKxold = sumKx
                func_v = eos.lnfugcoef_func(IDs,EoS,MR,P,T,y,kij, 1,Vv,en_auto,beta_auto,CR,SM,it,pt,r_data) #Vapor
                lnfugcoef_v = func_v[0]
                Vv = func_v[1]
                K = np.exp(lnfugcoef_l-lnfugcoef_v)
                Kx = K*x
                sumKx = np.sum(Kx)
                erry = abs(sumKx-sumKxold)/sumKx
                ity = ity+1
                #print 'erry',erry,y,K
                #print 'y inside erry',y
            #Iteration y end=================================================
            
            errK = abs(sumKx-1)/sumKx
            y = Kx/sumKx
            P = P*sumKx
            it = it+1
            #print 'y',y
            #print 'P',P,sumKx
            #print 'it',it
            print 'errK---------------',errK,tolK,P,y,it
            if it>=itmax:
                    Psat = np.array(correlations.Psat_antoine(IDs,T))
                    P = np.dot(Psat,x)
                    y = Psat*x/P
                    Vv = R*T/P #Not going to be used, just starting
                    Vl = 0.99  #Not going to be used, just starting
            #raw_input('...')
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
        if stepx==5E-3 and pt==0:
            x[0] = 5E-3
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

#Given vector of T, iterate over T to calculate TV envelope--------------------------------------
def TV_envelope(Tvec,IDs,EoS,MR,kij,nc,AR,CR,SM,r_data):

    Pvec = []
    xvec = []
    yvec = []
    Vlvec = []
    Vvvec = []
    
    #Definitions--------------------------------------------------
    #Main iteration conditions
    x = np.array([0.0001,0.9999]) #x array
    xf = 1.0                    #Main stop condition
    stepx = 1e-3                #Main step
    it = 0                      #Iteration counter
    ity = 0                     #Iteration over y loop counter
    pt = 0                      #Point counter
    itmax = 100                 #Max number of iterations
    
    #tolerances
    tolK = 1e-2             #iteration over Kx
    toly = 1e-3             #iteration over y
    
    #CPA auto-association configurations
    auto = []
    auto = association.CPA_auto(AR,nc,IDs)
    en_auto = auto[0]
    beta_auto = auto[1]
    #=============================================================

    #Main iteration start, range T------------------------------------------
    for i in range(0,len(Tvec)):
        x = np.array([0.0001,0.9999]) #x array
        T = Tvec[i]

        #Initial guess------------------------------------------------
        #Saturation pressure
        Psat = np.array(correlations.Psat_antoine(IDs,T))

        #Initial guess
        if pt==0:
            P = np.dot(Psat,x)
            y = Psat*x/P
            Vv = R*T/P #Not going to be used, just starting
            Vl = 0.99  #Not going to be used, just starting
        #print Psat,x,P,y
        #=============================================================

        #Iteration Kx start--------------------------------------------------
        errK = tolK+1 #Force enter iteration
        while errK>tolK and it<itmax:
            func_l = eos.lnfugcoef_func(IDs,EoS,MR,P,T,x,kij,-1,Vl,en_auto,beta_auto,CR,SM,it,pt,r_data) #Liquid
            lnfugcoef_l = func_l[0]
            Vl = func_l[1]
            func_v = eos.lnfugcoef_func(IDs,EoS,MR,P,T,y,kij, 1,Vv,en_auto,beta_auto,CR,SM,it,pt,r_data) #Vapor
            lnfugcoef_v = func_v[0]
            Vv = func_v[1]
            K = np.exp(lnfugcoef_l-lnfugcoef_v)
            Kx = K*x
            sumKx = np.sum(Kx)
            
            #print 'phil = ',np.exp(lnfugcoef_l),Vl
            #print 'phiv = ',np.exp(lnfugcoef_v),Vv
            #print 'sumKx = ',sumKx,K,Kx
            """
            #Iteration y start-----------------------------------------------
            erry = toly+1
            ity = 0
            while (erry>toly and ity<200) or ity<2:
                y = Kx/sumKx
                sumKxold = sumKx
                func_v = eos.lnfugcoef_func(IDs,EoS,MR,P,T,y,kij, 1,Vv,en_auto,beta_auto,CR,SM,it,pt,r_data) #Vapor
                lnfugcoef_v = func_v[0]
                Vv = func_v[1]
                K = np.exp(lnfugcoef_l-lnfugcoef_v)
                Kx = K*x
                sumKx = np.sum(Kx)
                erry = abs(sumKx-sumKxold)/sumKx
                ity = ity+1
                print 'erry',erry,y,K
                #print 'y inside erry',y
            #Iteration y end=================================================
            """
            errK = abs(sumKx-1)/sumKx
            y = Kx/sumKx
            P = P*sumKx
            it = it+1
            #print 'y',y
            #print 'P',P,sumKx
            #print 'it',it
            #print 'errK---------------',errK,tolK,P,y
            #raw_input('...')
        #Iteration Kx end====================================================
        
        #Save P, x1 and y1
        Pvec.append(P)
        xvec.append(x[0])
        yvec.append(y[0])
        Vlvec.append(Vl)
        Vvvec.append(Vv)
        
        #Update main condition
        #print'======='
        #print 'x',x
        #print 'P',P
        #print 'it',it
        #print 'y',y
        #input('End a point--------------')
        #print x[0],y[0],P,it
        x[0] = x[0]+stepx
        if stepx==5E-3 and pt==0:
            x[0] = 5E-3
        x[1] = 1-x[0]
        pt = pt+1       #add counter pt
        it = 0          #zero counter it
        
        #print 'pt',pt,T,P,1/Vl
    #Main iteration end, range T============================================
    
    #input('...')
    out = []
    out.append(Pvec)
    out.append(xvec)
    out.append(yvec)
    out.append(Vlvec)
    out.append(Vvvec)
    
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

    #--------
    rhomax = rho[n-1]
    rhomaxmax = rhomax
    rho2 = np.empty((10000))
    #print n,rhomax
    for i in range(0,10000):
        rho2[i] = i*rhomax/10000
        #print rho2[i],i,float(i/10000)
    rho2[0] = 1e-5
    f2 = splev(rho2,fspl)        #Evaluate Cubic Spline First derivative
    f = f2
    rho = rho2
    fspl = splrep(rho,f,k=3)         #Cubic Spline Representation
    n = 10000
    #--------

    u = splev(rho,fspl,der=1)        #Evaluate Cubic Spline First derivative
    
    drho = rho[n/2]-rho[n/2-1]
    #for i in range(1,n-2):
    #    u[i] = (f[i+1]-f[i-1])/(2*drho)
    #u[n-1] = (f[n-1]-f[n-2])/drho
    #u[0] = (f[1]-f[0])/drho
    
    #V = 1/rho2
    #A = f2/rho2
    #print V,len(V)
    #print A,len(A)
    #dAdV = np.gradient(A,V,edge_order=2)
    #Pp = -dAdV
    
    #Calcule pressure
    n = 10000
    P = -f+rho*u
    a = (P[n/2-15]-P[n/2+15])/(rho[n/2-15]-rho[n/2+15])
    b = P[n/2+15]-a*rho[n/2+15]
    i = n/2-15
    while i<(n/2+15):
        P[i]=a*rho[i]+b
        i = i + 1

    #data.modify_isotherm(P,rho) #Change reported isotherm
    #raw_input('isotherm modified')

    #P = Pp

    #First derivative
    #drho = rho[1]-rho[0]
    Pspl = splrep(rho,P,k=3)         #Cubic Spline Representation
    P = splev(rho,Pspl,der=0)
    dPdrho = splev(rho,Pspl,der=1)        #Evaluate Cubic Spline First derivative

    #for i in range(0,n):
    #    print i,rho[i],P[i],dPdrho[i]
    #plt.plot(rho,P)
    #plt.ylim(0,5)
    #plt.show()

    #Find max and min pressure of isotherm inside binodal curve
    max1 = int(numerical.bin_max(dPdrho))
    if max1==n:
        dens = []
        dens.append(0)
        dens.append(0)
        dens.append(0)
        dens.append(0)
        dens.append(0)
        dens.append(0)
        return dens
    min1 = int(numerical.bin_min(dPdrho))
    rhomax = rho[max1]
    rhomin = rho[min1]
    Pmax = P[max1]
    Pmin = P[min1]
    min2 = min1+30
    max2 = max1-30

    #print Pmax,Pmin,rhomax,max2,max1,min1,min2

    #Pmin = Pmax+1 #method below seems always better, trying forcing it everytime
    if Pmin>Pmax:
        min1 = numerical.bin_min_seed(dPdrho,max1)
        if max1==n:
            dens = []
            dens.append(0)
            dens.append(0)
            dens.append(0)
            dens.append(0)
            dens.append(0)
            dens.append(0)
            return dens
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
        min2 = min2 + 2

    if min2>10000:
        min2 = 10000

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
    rho1 = np.amin(Pf1roots)
    rho2 = np.amax(Pf2roots)

    #print 'roots',Pf1roots,Pf2roots
    #input('...')
    
    #i = 0
    #rho1=rho[0]-1
    #while rho1<rho[0] and rho1<0:
    #    rho1 = Pf1roots[i]
    #    i = i+1
    #if Pf1roots[0]>Pf2roots[0]:
    #    rho1 = 0.1
    #rho1 = numerical.falsi_spline(rho,Pf1,rho[0],rhomax,1e-7)
     
    
    #i = 0
    #rho2 = rhomin-1
    #while rho2<rhomin:
    #    rho2 = Pf2roots[i]
    #    i = i+1
    #rho2 = numerical.falsi_spline(rho,Pf2,rhomin,rho[min2],1e-7)

    #print 'rho1,rho2',rho1,rho2
    #print 'max',rho[0],rhomax
    #print 'min',rhomin,rho[min2]
    #print 'Pmin/max',Pmin,Pmax

    #plt.plot(rho,P,'r-')
    #plt.plot(rho,Pp,'g-.')
    #plt.ylim(0,8)
    #plt.show()

    #Solve newton-raphson system
    tol = 1e-10
    drho1 = tol+1
    drho2 = tol+1
    du = tol+1
    dP = tol+1
    drho2old = tol+1
    drho1old = tol+1
    stop = 1.0
    counter = 0

    uspl = splrep(rho,u,k=3)
    Pspl = splrep(rho,P,k=3)
    fspl1 = splev(rho,fspl)
    uspl1 = splev(rho,uspl)
    Pspl1 = splev(rho,Pspl)
    dudrho = splev(rho,uspl,der=1)
    dPdrho = splev(rho,Pspl,der=1)
    Nitmax = 1000
    Nit = 0

    #print '................'
    #print rho1,rho2

    #plt.plot(rho,P)
    #plt.ylim(8,10)
    #plt.show()

    #plt.plot(rho,dPdrho)
    #plt.ylim(-0.001,0.001)
    #plt.show()
    
    """
    start1 = timeit.default_timer()
    rho1_orig = rho1
    rho2_orig = rho1
    
    def uP_calc(rhov,rho,uspl1,Pspl1):
        rho1 = rhov[0]
        rho2 = rhov[1]
        u1 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho1)
        u2 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho2)
        P1 = InterpolatedUnivariateSpline(rho,Pspl1,k=3)(rho1)
        P2 = InterpolatedUnivariateSpline(rho,Pspl1,k=3)(rho2)
        du = abs(u1-u2)
        dP = abs(P1-P2)
        uP = du+dP
        return uP
    rhos = []
    rhos.append(rho1)
    rhos.append(rho2)
    wow = optimize.minimize(uP_calc,rhos,args=(rho,uspl1,Pspl1),method='Nelder-Mead',options={'maxiter':10000})
    rho1 = wow.x[0]
    rho2 = wow.x[1]
    #print 'Nelder:',rho1,rho2
    #rho1 = rho1_orig
    #rho2 = rho2_orig
    """
    
    while (abs(du)>tol or abs(dP)>tol) and (abs(Pmax-Pmin)>1e-4) and (Nit<Nitmax):
    #while (abs(drho1)>tol or abs(drho2)>tol) and (abs(Pmax-Pmin)>1e-3) and (Nit<Nitmax):
        Nit = Nit+1
        rho1 = rho1 + stop*drho1
        rho2 = rho2 + stop*drho2

        #f1 = InterpolatedUnivariateSpline(rho,fspl1,k=3)(rho1)
        #f2 = InterpolatedUnivariateSpline(rho,fspl1,k=3)(rho2)
        u1 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho1)
        u2 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho2)
        #P1 = -f1+rho1*u1
        #P2 = -f2+rho2*u2
        P1 = InterpolatedUnivariateSpline(rho,Pspl1,k=3)(rho1)
        P2 = InterpolatedUnivariateSpline(rho,Pspl1,k=3)(rho2)

        du1 = InterpolatedUnivariateSpline(rho,dudrho,k=3)(rho1)
        du2 = InterpolatedUnivariateSpline(rho,dudrho,k=3)(rho2)
        dP1 = InterpolatedUnivariateSpline(rho,dPdrho,k=3)(rho1)
        dP2 = InterpolatedUnivariateSpline(rho,dPdrho,k=3)(rho2)
        detJ = -dP2*du1+dP1*du2

        drho2 = -du1/detJ*(P1-P2)+dP1/detJ*(u1-u2)
        drho1 = -du2/detJ*(P1-P2)+dP2/detJ*(u1-u2)

        #if drho1>500:
        #    drho1=500
        #if drho1<-500:
        #    drho1=-500
        #if drho2>500:
        #    drho2=500
        #if drho2<-500:
        #    drho2=-500
        if rho1>rhomax:
            #if P2>Pmax:
            #    Pf1 = P - Pmax
            #else:
            #    Pf1 = P - P2
            rho1 = numerical.falsi_spline(rho,Pf1,rho[0],rhomax,1e-5)
        if rho2<rhomin:
            #if P1<Pmin:
            #    Pf2 = P - Pmin
            #else:
            #    Pf2 = P - P1
            rho2 = numerical.falsi_spline(rho,Pf2,rhomin,rho[min2],1e-5)
        
        du = abs(u1-u2)
        dP = abs(P1-P2)
        #print P1,P2,rho1,rho2,du,dP,drho1,drho2,stop,Nit

        if counter>0 and (drho1>drho1old and drho2>drho2old):
            rho1 = rho1 - stop*drho1
            rho2 = rho2 - stop*drho2
            stop = stop/1.005 #Break
            rho1 = rho1 + stop*drho1/2
            rho2 = rho2 + stop*drho2/2
            #print stop,counter
        
        counter = counter+1
        drho1old = drho1
        drho2old = drho2
        duold = du
        dPold = dP

    if abs(Pmax-Pmin)<1e-4:
        rho1 = (rhomax+rhomin)/2
        rho2 = rho1
    
    #f1 = InterpolatedUnivariateSpline(rho,fspl1,k=3)(rho1)
    #f2 = InterpolatedUnivariateSpline(rho,fspl1,k=3)(rho2)
    u1 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho1)
    u2 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho2)
    #P1 = -f1+rho1*u1
    #P2 = -f2+rho2*u2
    P1 = InterpolatedUnivariateSpline(rho,Pspl1,k=3)(rho1)
    P2 = InterpolatedUnivariateSpline(rho,Pspl1,k=3)(rho2)

    dens = []
    dens.append(rho1)
    dens.append(rho2)
    dens.append(P1)
    dens.append(P2)
    dens.append(u1)
    dens.append(u2)
    #input('...')
    
    return dens
#======================================================================================

#Given pure component isotherm, calculates phase coexistence densities-----------------
def coexistence_dens2(rho1,f1,rhol_guess,rhov_guess):

    rho1 = np.array(rho1)
    f1 = np.array(f1)
    n = f1.shape[1]
    rho = rho1.flatten()
    f = f1.flatten()

    #Spline to get chemical potential
    fspl = splrep(rho,f,k=3)         #Cubic Spline Representation

    #--------
    rhomax = rho[n-1]
    rhomaxmax = rhomax
    rho2 = np.empty((10000))
    #print n,rhomax
    for i in range(0,10000):
        rho2[i] = i*rhomax/10000
        #print rho2[i],i,float(i/10000)
    rho2[0] = 1e-5
    f2 = splev(rho2,fspl)        #Evaluate Cubic Spline First derivative
    f = f2
    rho = rho2
    fspl = splrep(rho,f,k=3)         #Cubic Spline Representation
    n = 10000
    #--------

    u = splev(rho,fspl,der=1)        #Evaluate Cubic Spline First derivative
    
    drho = rho[n/2]-rho[n/2-1]
    
    #Calcule pressure
    n = 10000
    P = -f+rho*u
    a = (P[n/2-15]-P[n/2+15])/(rho[n/2-15]-rho[n/2+15])
    b = P[n/2+15]-a*rho[n/2+15]
    i = n/2-15
    while i<(n/2+15):
        P[i]=a*rho[i]+b
        i = i + 1

    #First derivative
    #drho = rho[1]-rho[0]
    Pspl = splrep(rho,P,k=3)         #Cubic Spline Representation
    P = splev(rho,Pspl,der=0)
    dPdrho = splev(rho,Pspl,der=1)        #Evaluate Cubic Spline First derivative

    #Solve newton-raphson system
    tol = 1e-10
    drho1 = tol+1
    drho2 = tol+1
    du = tol+1
    dP = tol+1
    drho2old = tol+1
    drho1old = tol+1
    stop = 1.0
    counter = 0

    uspl = splrep(rho,u,k=3)
    Pspl = splrep(rho,P,k=3)
    fspl1 = splev(rho,fspl)
    uspl1 = splev(rho,uspl)
    Pspl1 = splev(rho,Pspl)
    dudrho = splev(rho,uspl,der=1)
    dPdrho = splev(rho,Pspl,der=1)
    Nitmax = 1000
    Nit = 0
    
    drho1 = 0.0
    drho2 = 0.0
    rho1 = rhov_guess
    rho2 = rhol_guess
    
    while (abs(du)>tol or abs(dP)>tol) and (Nit<Nitmax):
    #while (abs(du)>tol or abs(dP)>tol) and (abs(Pmax-Pmin)>1e-5) and (Nit<Nitmax):
        Nit = Nit+1
        rho1 = rho1 + stop*drho1
        rho2 = rho2 + stop*drho2

        u1 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho1)
        u2 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho2)
        P1 = InterpolatedUnivariateSpline(rho,Pspl1,k=3)(rho1)
        P2 = InterpolatedUnivariateSpline(rho,Pspl1,k=3)(rho2)

        du1 = InterpolatedUnivariateSpline(rho,dudrho,k=3)(rho1)
        du2 = InterpolatedUnivariateSpline(rho,dudrho,k=3)(rho2)
        dP1 = InterpolatedUnivariateSpline(rho,dPdrho,k=3)(rho1)
        dP2 = InterpolatedUnivariateSpline(rho,dPdrho,k=3)(rho2)
        detJ = -dP2*du1+dP1*du2

        drho2 = -du1/detJ*(P1-P2)+dP1/detJ*(u1-u2)
        drho1 = -du2/detJ*(P1-P2)+dP2/detJ*(u1-u2)
        
        du = abs(u1-u2)
        dP = abs(P1-P2)

        if counter>0 and (drho1>drho1old and drho2>drho2old):
            rho1 = rho1 - stop*drho1
            rho2 = rho2 - stop*drho2
            stop = stop/1.005 #Break
            rho1 = rho1 + stop*drho1/2
            rho2 = rho2 + stop*drho2/2
        
        counter = counter+1
        drho1old = drho1
        drho2old = drho2
        duold = du
        dPold = dP
    
    u1 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho1)
    u2 = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho2)
    P1 = InterpolatedUnivariateSpline(rho,Pspl1,k=3)(rho1)
    P2 = InterpolatedUnivariateSpline(rho,Pspl1,k=3)(rho2)

    dens = []
    dens.append(rho1)
    dens.append(rho2)
    dens.append(P1)
    dens.append(P2)
    dens.append(u1)
    dens.append(u2)
    #input('...')
    
    return dens
#======================================================================================

#Given pure component isotherm, calculates phase coexistence densities-----------------
def inflex_search(rho1,f1,tolI):

    rho1 = np.array(rho1)
    f1 = np.array(f1)
    n = f1.shape[1]
    rho = rho1.flatten()
    f = f1.flatten()

    #Spline to get chemical potential
    fspl = splrep(rho,f,k=3)         #Cubic Spline Representation

    #--------
    rhomax = rho[n-1]
    rhomaxmax = rhomax
    rho2 = np.empty((10000))
    #print n,rhomax
    for i in range(0,10000):
        rho2[i] = i*rhomax/10000
    rho2[0] = 1e-5
    f2 = splev(rho2,fspl)        #Evaluate Cubic Spline First derivative
    f = f2
    rho = rho2
    fspl = splrep(rho,f,k=3)         #Cubic Spline Representation
    n = 10000
    #--------

    u = splev(rho,fspl,der=1)        #Evaluate Cubic Spline First derivative
    
    #drho = rho[n/2]-rho[n/2-1]
    
    #Calcule pressure
    n = 10000
    P = -f+rho*u
    #a = (P[n/2-15]-P[n/2+15])/(rho[n/2-15]-rho[n/2+15])
    #b = P[n/2+15]-a*rho[n/2+15]
    #i = n/2-15
    #while i<(n/2+15):
    #    P[i]=a*rho[i]+b
    #    i = i + 1

    Pspl = splrep(rho,P,k=3)         #Cubic Spline Representation
    P = splev(rho,Pspl,der=0)
    dPdrho = splev(rho,Pspl,der=1)        #Evaluate Cubic Spline First derivative
    d2Pdrho2 = splev(rho,Pspl,der=2)
    
    #data.modify_isotherm2(rho,P,dPdrho,d2Pdrho2)
    #raw_input('......')
    
    """
    #Create figure
    fig, ax = plt.subplots(1,1)
    
    #Data
    plt.plot(rho,P)
    
    #Axis size
    #plt.ylim([ydata[0]-20,ydata[size-1]+20])
    
    #Save figure with 1000 dpi resolution
    figsavedir = str('../output/test.png')
    fig.savefig(figsavedir, dpi=1000)
    """
    
    #for i in range(0,len(P)):
    #    if i%50==0:
    #        print 'P',i,rho[i],P[i],dPdrho[i],d2Pdrho2[i]
    
    """
    inflex = False

    dens = []
    minP = np.argmin(d2Pdrho2)
    dist = d2Pdrho2[minP]
    if d2Pdrho2[minP]>0: #There is no inflexion
        rhoc = rho2(minP)
        Pc = P(minP)
        inflex = True
        dens.append(rhoc)
        dens.append(rhoc)
        dens.append(Pc)
        dens.append(Pc)
        dens.append(0)
        dens.append(dist)
        dens.append(inflex)
    else:
        dens.append(0)
        dens.append(0)
        dens.append(0)
        dens.append(0)
        dens.append(0)
        dens.append(dist)
        dens.append(inflex)
    """
    
    dens = []
    tol = 1e-10

    #print 'search minimum'
    min_d2P = numerical.bisect_spline(rho,d2Pdrho2,rho[n/10],rho[len(rho)-n/10],tol,200)
    dP_min = InterpolatedUnivariateSpline(rho,dPdrho,k=3)(min_d2P)
    inflex = False
    inflex2 = False
    
    #print 'min=',min_d2P,dP_min
    
    if abs(dP_min)<=tolI:  #T=Tc
        rhoc = min_d2P
        Pc = InterpolatedUnivariateSpline(rho,P,k=3)(min_d2P)
        inflex = True
        inflex2 = False
        dens.append(rhoc)
        dens.append(rhoc)
        dens.append(Pc)
        dens.append(Pc)
        dens.append(0)
        dens.append(dP_min)
        dens.append(inflex)
        dens.append(inflex2)
    if abs(dP_min)>tolI:                   #T<Tc or T>Tc
        inflex = False
        inflex2 = False
        dens.append(0)
        dens.append(0)
        dens.append(0)
        dens.append(0)
        dens.append(0)
        dens.append(dP_min)
        dens.append(inflex)
        dens.append(inflex2)
    
    return dens
#======================================================================================

#Given initial T, using renormalization method, build P-rho envelope-------------------
def PV_envelope(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n):
    
    env = []
    Tv = []
    rhov = []
    rhol = []
    Pv = []
    
    print 'T:   dens_vap:   dens_liq:   P:'
    while T<=Tfinal:
        ren = renormalization.renorm(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
        dens = coexistence_dens(ren[2],ren[0])
        Tv.append(T)
        rhov.append(dens[0])
        rhol.append(dens[1])
        Pv.append(dens[2])
        print T,dens[0],dens[1],dens[2]
        if abs(dens[0]-dens[1])<1.0:
	        break
        T = T + stepT
        
    env.append(Tv)
    env.append(rhov)
    env.append(rhol)
    env.append(Pv)
    return env
#======================================================================================  

#Given initial T, using renormalization method, build P-rho envelope-------------------
def PV_vec_envelope(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n):
    
    env = []
    Tv = []
    rhov = []
    rhol = []
    Pv = []
    size = len(T)

    #print 'T:   dens_vap:   dens_liq:   P:'
    for i in range(0,size):
        ren = renormalization.renorm(EoS,IDs,MR,T[i],nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
        dens = coexistence_dens(ren[2],ren[0])
        Tv.append(T)
        rhov.append(dens[0])
        rhol.append(dens[1])
        Pv.append(dens[2])
        #print i,T[i],dens[1],dens[2]

    env.append(Tv)
    env.append(rhov)
    env.append(rhol)
    env.append(Pv)
    return env
#======================================================================================  

#Given initial T, using renormalization method, build P-rho envelope-------------------
def PV_findTc_envelope(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,crit_bool,L__est,phi__est):
    
    env = []
    Tv = []
    rhov = []
    rhol = []
    Pv = []
    Fobj = 2.0
    i = 0
    h = 0
    step = 1.0
    flag0 = False
    TTc = Tfinal
    print 'param:',L__est,phi__est
    
    print 'T:   dens_vap:   dens_liq:   Fobj:   Fobjder:   step:'
    while Fobj>0.1:
        ren = renormalization.renorm(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,L__est,phi__est)
        dens = coexistence_dens(ren[2],ren[0])
        if (dens[2]!=0):
            Tv.append(T)
            rhov.append(dens[0])
            rhol.append(dens[1])
            Pv.append(dens[2])
        Fobj = abs(dens[0]-dens[1])

        #ren_p = renormalization.renorm(EoS,IDs,MR,T*(1+h),nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n)
        #dens_p = coexistence_dens(ren_p[2],ren_p[0])
        #Fobj_plus = abs(dens_p[0]-dens_p[1])

        #ren_m = renormalization.renorm(EoS,IDs,MR,T*(1-h),nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n)
        #dens_m = coexistence_dens(ren_m[2],ren_m[0])
        #Fobj_minus = abs(dens_m[0]-dens_m[1])
        
        #Fobj_derivative = (Fobj_plus-Fobj_minus)/(2*T*h)
        #if Fobj_derivative==0:
        #    Fobj_derivative=-100
        #step = -Fobj/Fobj_derivative

        #print T,dens[0],dens[1],Fobj,Fobj_plus,Fobj_minus,Fobj_derivative,step
        print T,dens[0],dens[1],dens[2]

        if dens[2]!=0:
            flag0 = False
        if flag0 == True:
            T = T - 1.0
        if dens[2]==0:
            TTc = T
            T = T - step
            #print 'TTc',TTc
            Fobj = 50
            if flag0 == True:
                T = T - 1.0
            else:
                flag0 = True

        if Fobj>3000:
            step = 2.0#1.0
        if Fobj>1000 and Fobj<=3000:
            step = 1.0#0.5
        if Fobj>100 and Fobj<=1000:
            step = 0.5#0.05
        if Fobj>10 and Fobj<=100:
            step = 0.1#0.01
        if Fobj<10:
            step = 0.01#0.001
        step = Fobj/2e3
        T = T + step
        i = i+1
    
    if crit_bool==True:
        crit_exp = renormalization.critical_exponents(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,L__est,phi__est,T-step,dens[0])
        report_crit(T-step,dens[2],dens[0],IDs,L__est,phi__est,crit_exp[0])
    else:
        report_crit(T-step,dens[2],dens[0],IDs,0,0,0)
    env.append(Tv)
    env.append(rhov)
    env.append(rhol)
    env.append(Pv)
    if crit_bool==True:
        env.append(crit_exp[0])
    else:
        env.append(0.000)
    return env
#======================================================================================  

#Given initial T, using renormalization method, build P-rho envelope-------------------
def PV_findTc2_envelope(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,L__est,phi__est):
    
    env = []
    Tv = []
    rhov = []
    rhol = []
    Pv = []
    Fobj = 2.0
    i = 0
    h = 0
    step = 1.0
    flag0 = False
    flagI = False
    tol = 1e-8
    
    print 'T:   dens_vap:   dens_liq:   Fobj:   Fobjder:   step:'
    while flagI==False:
        ren = renormalization.renorm(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,L__est,phi__est)
        ans = inflex_search(ren[2],ren[0],1e-10)
        Fobj = ans[5]
        flagI = ans[6]

        print T,ans[0],ans[1],ans[2],ans[5],ans[6],i

        if Fobj>1e-5:
            step = -1e0
        if Fobj>1e-7 and Fobj<=1e-5:
            step = -1e-1
        if Fobj>1e-10 and Fobj<=1e-7:
            step = -1e-2

        if Fobj<-1e-5:
            step = 1e0
        if Fobj<-1e-7 and Fobj>=-1e-5:
            step = 1e-1
        if Fobj<-1e-10 and Fobj>=-1e-7:
            step = 1e-2
        step = -4e4*ans[5]
        #if abs(ans[5])<1e-10:
        #    step = -3e6*ans[5]
        T = T + step
        i = i+1
    
    if estimate_bool==True:
        crit_exp = renormalization.critical_exponents(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,L__est,phi__est,T-step,ans[0])
        report_crit(T-step,ans[2],ans[0],IDs,L__est,phi__est,crit_exp[0])
    else:
        report_crit(T-step,ans[2],ans[0],IDs,0,0,0)
    env.append(T-step)
    env.append(ans[0])
    env.append(ans[0])
    env.append(ans[2])
    if estimate_bool==True:
        env.append(crit_exp[0])
    else:
        env.append(0.000)
    return env
#======================================================================================  

#Given initial T, using renormalization method, build P-rho envelope-------------------
def PV_findTc3_envelope(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,crit_bool,L__est,phi__est):
    
    env = []
    Tv = []
    rhov = []
    rhol = []
    Pv = []
    Fobj = 2.0
    i = 0
    h = 0
    step = 1.0
    flag0 = False
    flagmax = False
    TT = T
    Tmax = 0
    Fobjold = 1./3.
    Tmin = 0
    #print 'param:',L__est,phi__est
    rholl = 1.0
    rhovv = 1.0
    
    #print 'T:   dens_vap:   dens_liq:   Fobj:   Fobjder:   step:'
    while Fobj>0.1:
        ren = renormalization.renorm(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,L__est,phi__est)
        #if i==0 or rholl==0:
        dens = coexistence_dens(ren[2],ren[0])
        #if (i>0 and rholl!=0):
        #    dens = coexistence_dens2(ren[2],ren[0],rholl,rhovv)
        if (dens[2]!=0):
            Tv.append(T)
            rhov.append(dens[0])
            rhol.append(dens[1])
            Pv.append(dens[2])
            rhovv = dens[0]
            rholl = dens[1]
        Fobj = abs(dens[0]-dens[1])
        print T,dens[0],dens[1],dens[2],Tmax,Tmin
        #raw_input('...')

        if dens[2]!=0:
            flag0 = False
            Tmin = T
        if dens[2]==0:
            flag0 = True

        if flag0==True:
            Tmax = T
            T = T-step/2
            Tmin = T
            flagmax = True
            Fobj = 100

        if flag0==False and flagmax==False:
            step = Fobj/2e3
            T = T + step

        if flag0==False and flagmax==True:
            step = (Tmax-Tmin)
            T = T + step/2
            T = Tmin+step/2
 
        if Fobjold==Fobj:
            T = T+0.1
        Fobjold = Fobj

        if abs(dens[1])>1e10:
            T = TT+1.0

        i = i+1
    
    if crit_bool==True:
        crit_exp = renormalization.critical_exponents(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,L__est,phi__est,T-step,dens[0])
        report_crit(T-step,dens[2],dens[0],IDs,L__est,phi__est,crit_exp[0])
    else:
        report_crit(T-step,dens[2],dens[0],IDs,0,0,0)
    env.append(Tv)
    env.append(rhov)
    env.append(rhol)
    env.append(Pv)
    if crit_bool==True:
        env.append(crit_exp[0])
    else:
        env.append(0.000)
    return env
#======================================================================================  

#Report found critical data------------------------------------------------------------
def report_crit(Tc,Pc,rhoc,IDs,Lc,phic,beta):
    title = 'critical.log'
    if Lc==0 and phic==0:
        Lc = float(data.L(IDs)[0])
        phic = float(data.phi(IDs)[0])
    savedir = str('../output/%s' %title)
    with open(savedir,'a') as file:
        lin1 = [str(round(Tc,9)),str(round(Pc,9)),str(round(rhoc,9)),str('%.3E' %Lc),str(round(phic,9)),str(round(beta,9)),str(round(IDs[0],9))]
        lin = ('%s\n' % ';'.join(lin1))
        file.write(lin)
#======================================================================================

#Report found critical data------------------------------------------------------------
def report_crit_map(Lc,phic,Tc,Pc,rhoc,IDs):
    title = 'critical_map.log'
    savedir = str('../output/%s' %title)
    for i in range(0,len(Lc)):
        for j in range(0,len(phic)):
            with open(savedir,'a') as file:
                lin1 = [str(round(Tc[i][j],9)),str(round(Pc[i][j],9)),str(round(rhoc[i][j],9)),str('%.2E' %Lc[i]),str(round(phic[j],9))]
                lin = ('%s\n' % ';'.join(lin1))
                file.write(lin)
#======================================================================================

#Given initial T, using renormalization method, estimate L and phi parameters----------
def PV_estimate_Tc_envelope(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n):
    
    env = []
    Tv = []
    rhov = []
    rhol = []
    Pv = []

    diffs = []
    Ls = []
    phis = []
    Fobj_minus_old = 1e-5

    tol = 1e-7
    Tinit = T
    Tfinal = 1000.0
    phi = 1.0
    L = 1.0e-10
    Lsize = 60
    phisize = 50
    Ls = np.linspace(3.00e-10,9.00e-10,Lsize)
    phis = np.linspace(0.01,10.0,phisize)
    Tcmap = np.empty((Lsize,phisize))
    Pcmap = np.empty((Lsize,phisize))
    rhocmap = np.empty((Lsize,phisize))
    estimate_bool = True
    crit_bool = True
    i = 0
    print 'L  |  phi  |  Tc  |  Pc  |  rhoc'
    for j in range(0,phisize):
        T = Tinit
        for i in range(0,Lsize):
            env = PV_findTc_envelope(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,estimate_bool,crit_bool,Ls[i],phis[j])
            size = len(env[0])
            Tc = env[0][size-1]
            rhoc = env[1][size-1]
            Pc = env[3][size-1]
            #report_crit(Tc,Pc,rhoc,IDs,Ls[i],phis[j])
            Tcmap[i][j] = Tc
            Pcmap[i][j] = Pc
            rhocmap[i][j] = rhoc
            T = Tc-1
            print Ls[i],phis[j],Tcmap[i][j],Pcmap[i][j],rhocmap[i][j],env[4]

    report_crit_map(Ls,phis,Tcmap,Pcmap,rhocmap,IDs)
    crit_map = []
    crit_map.append(Tcmap)
    crit_map.append(Pcmap)
    crit_map.append(rhocmap)
    #env.append(Tv)
    #env.append(rhov)
    #env.append(rhol)
    #env.append(Pv)
    return crit_map
#====================================================================================== 

#Given initial T, using renormalization method, calculate derivative properties--------
def PV_deriv_calc_envelope(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,P,dtype):
    
    fres0 = []
    fres1 = []
    fres2 = []
    
    f0 = []
    f1 = []
    f2 = []
    
    P0 = []
    P1 = []
    P2 = []
    
    rho0 = []
    rho1 = []
    rho2 = []
    
    T0 = []
    T1 = []
    T2 = []
    
    h = 5e-2

    step = 1.0
    finalT = T
    
    if EoS==2 or EoS==4 or EoS==6:
        #print 'T:   dens_vap:   dens_liq:   P:'
        while T<=finalT:

            T0.append(T-T*h)
            T1.append(T)
            T2.append(T+T*h)
        
            ren = renormalization.renorm(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
            fres1.append(ren[7])
            rho1.append(ren[2])
            P1.append(ren[8])
            f1.append(ren[0])
            #print 'central',T,dens[0],dens[1],dens[2],Fobj

            ren = renormalization.renorm(EoS,IDs,MR,T+T*h,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
            fres2.append(ren[7])
            rho2.append(ren[2])
            P2.append(ren[8])
            f2.append(ren[0])
            #print 'plus',T,dens[0],dens[1],dens[2],Fobj_plus

            ren = renormalization.renorm(EoS,IDs,MR,T-T*h,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
            fres0.append(ren[7])
            rho0.append(ren[2])
            P0.append(ren[8])
            f0.append(ren[0])
            #print 'minus',T,dens[0],dens[1],dens[2],Fobj_minus

            T = T+step
    
        T_list = []
        fres_list = []
        P_list = []
        rho_list = []
        f_list = []
    
        T_list.append(T0)
        T_list.append(T1)
        T_list.append(T2)
    
        fres_list.append(fres0)
        fres_list.append(fres1)
        fres_list.append(fres2)
        
        f_list.append(f0)
        f_list.append(f1)
        f_list.append(f2)
    
        P_list.append(P0)
        P_list.append(P1)
        P_list.append(P2)
    
        rho_list.append(rho0)
        rho_list.append(rho1)
        rho_list.append(rho2)
    
        if dtype==1: #Isothermal
            der_prop = derivativeprop.calc_isothermal_dev_prop_pure(T_list,f_list,P_list,rho_list,T*h,IDs)
        if dtype==2: #isobaric
            der_prop = derivativeprop.calc_isobaric_dev_prop_pure(T_list,f_list,P_list,rho_list,T*h,IDs,P)
        
    else:
        n = 500
        b = eos.b_calc(IDs,EoS)
        x = np.array([0.999,0.001])
        bmix = eos.bmix_calc(MR,b,x)
        a = eos.a_calc(IDs,EoS,T)
        amix = eos.amix_calc(MR,a,x,kij)
        V = np.empty((n))
        rho = np.empty((n))
        for i in range(0,n):
            rho[i] = np.array(float(i)/n/bmix)
        rho[0] = 1e-8
        V = 1/rho
        f = rho*R*T*np.log(rho/(1-rho*bmix))-rho*amix/bmix/np.sqrt(8)*np.log((1+rho*bmix*(1+np.sqrt(2)))/(1+rho*bmix*(1-np.sqrt(2))))
        A = f*V
        
        print T,R,bmix,amix
        
        #rho = np.array(menus.flatten(rho))
        #f = np.array(menus.flatten(f))
        
        fspl = splrep(rho,f,k=3)
        u = splev(rho,fspl,der=1)
        
        P = -f+rho*u
        
        der_prop = derivativeprop.calc_isothermal_dev_prop_pure_analytical(T,A,P,V,IDs,EoS,MR,kij)
    
    return der_prop
#====================================================================================== 

#Given initial T, using renormalization method, calculate entropy curve--------
def PV_entropy(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,P,dtype):
    
    fres0 = []
    fres1 = []
    fres2 = []
    
    f0 = []
    f1 = []
    f2 = []
    
    A0 = []
    A1 = []
    A2 = []
    
    P0 = []
    P1 = []
    P2 = []
    
    rho0 = []
    rho1 = []
    rho2 = []
    
    T0 = []
    T1 = []
    T2 = []
    
    h = 5e-2

    step = 1.0
    finalT = T
    
    if EoS==2 or EoS==4 or EoS==6:
        #print 'T:   dens_vap:   dens_liq:   P:'
        while T<=finalT:

            T0.append(T-T*h)
            T1.append(T)
            T2.append(T+T*h)
        
            ren = renormalization.renorm(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
            fres1.append(ren[7])
            rho1.append(ren[2])
            P1.append(ren[8])
            f1.append(ren[0])
            A1 = f1[0][0]/rho1[0][0]
            #print 'central',T,dens[0],dens[1],dens[2],Fobj

            ren = renormalization.renorm(EoS,IDs,MR,T+T*h,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
            fres2.append(ren[7])
            rho2.append(ren[2])
            P2.append(ren[8])
            f2.append(ren[0])
            A2 = f2[0][0]/rho2[0][0]
            #print 'plus',T,dens[0],dens[1],dens[2],Fobj_plus

            ren = renormalization.renorm(EoS,IDs,MR,T-T*h,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
            fres0.append(ren[7])
            rho0.append(ren[2])
            P0.append(ren[8])
            f0.append(ren[0])
            A0 = f0[0][0]/rho0[0][0]
            #print 'minus',T,dens[0],dens[1],dens[2],Fobj_minus

            T = T+step
    
        T_list = []
        fres_list = []
        P_list = []
        rho_list = []
        f_list = []
        A_list = []
    
        T_list.append(T0)
        T_list.append(T1)
        T_list.append(T2)
    
        fres_list.append(fres0)
        fres_list.append(fres1)
        fres_list.append(fres2)
        
        f_list.append(f0)
        f_list.append(f1)
        f_list.append(f2)
        
        A_list.append(A0)
        A_list.append(A1)
        A_list.append(A2)
    
        P_list.append(P0)
        P_list.append(P1)
        P_list.append(P2)
    
        rho_list.append(rho0)
        rho_list.append(rho1)
        rho_list.append(rho2)
    
        if dtype==1: #Isothermal
            s_prop = renormalization.calc_entropy(T_list,f_list,P_list,rho_list,T*h,IDs)
        if dtype==2: #isobaric
            s_prop = renormalization.calc_entropy(T_list,f_list,P_list,rho_list,T*h,IDs,P)
    
    return s_prop
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

#Report estimation history-------------------------------------------------------------
def report_param(data,filename):
    #df = pd.DataFrame(data)
    #df.to_csv(filename,index=False,header=False,index_label=False,sep=';')
    df_a = pd.DataFrame(data)
    df = pd.read_csv(filename,sep=';')
    with open(filename, 'a') as f:
        df_a.to_csv(f,index=False,header=False,index_label=False,sep=';')
#======================================================================================

#Report estimation history-------------------------------------------------------------
def report_df(data,filename):
    df_a = pd.DataFrame(data)
    df = pd.read_csv(filename,sep=';')
    with open(filename, 'a') as f:
        df_a.to_csv(f,index=False,header=False,index_label=False,sep=';')
#======================================================================================

#Define and handle which envelope to calculate-----------------------------------------
def calc_env(user_options,print_options,nc,IDs,EoS,MR,z,AR,CR,P,T,kij,auto,en_auto,beta_auto,SM,env_type):
    
    #INDEX START------------------------------
    #1 - PT Envelope for mixture
    #2 - PV Envelope for pure
    #3 - PV Envelope estimate
    #4 - Critical Point for pure
    #5 - Derivative properties
    #6 - Pxy Envelope for binary mixture
    #7 - Critical point for mixture
    #8 - Renormalization Parameters estimation
    #9 - Association Parameters estimation
    #10 - CPA Parameters estimation
    #INDEX END--------------------------------

    if env_type==1:
        #Calculate PT envelope with continuation method*********************************
        print '\nCalculating PT envelope'
        env_PT = PT_envelope(z,T,P,IDs,EoS,MR,kij,en_auto,beta_auto,CR,SM)
        print 'PT envelope calculated'

        print 'Creating PT report'
        reportname = str('PT_%s.csv' %('_'.join(print_options[1])))
        aa=['{:.4f}'.format(x) for x in print_options[4]]
        reportname = str('PT_%s.csv' %('_'.join(print_options[1]+[print_options[2]]+aa)))
        report_PT(env_PT,print_options,reportname,print_options)
        print ('Report %s saved successfully' %reportname)

        print 'Starting to plot PT envelope'
        title = str('PT envelope\n%s' %(' + '.join(print_options[1]+[print_options[2]]+aa)))
        figname = str('PT_%s.png' %('_'.join(print_options[1]+[print_options[2]]+aa)))
        zlist = " ,".join(format(x, ".2f") for x in z)
        boxtext = str('z= [%s]' %zlist)
        #boxtext = 'box'
        plot_PT(title,'T(K)','P(MPa)',figname,boxtext,env_PT[2],env_PT[3])
        print ('Figure %s saved successfully' %figname)
        #*******************************************************************************

    if env_type==5:
        #Calculate pure PV derivative properties*****************************************************
        nd = 1000
        nx = 200
        n = 8
        stepT = 0.5
        finalT = T
        print '\nCalculating pure isothermal derivative properties'
        dp_dat = PV_deriv_calc_envelope(EoS,IDs,MR,T,finalT,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,P,1) #Isothermal
        #dp_dat = PV_deriv_calc_envelope(EoS,IDs,MR,T,finalT,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,P,2) #Isobaric
        print 'pure isothermal derivative properties calculated'

        print 'Creating pure isothermal derivative properties report'
        reportname = str('Deriv_Prop_%s.csv' %('_'.join(print_options[1])))
        derivativeprop.report_isothermal_dev_prop_pure(reportname,dp_dat[0],dp_dat[1],dp_dat[2],dp_dat[3],dp_dat[4],dp_dat[5],dp_dat[6],
                                                     dp_dat[7],dp_dat[8],dp_dat[9],dp_dat[10],dp_dat[11],dp_dat[12],
                                                     dp_dat[13],dp_dat[14],dp_dat[15],print_options)
        print ('Report %s saved successfully' %reportname)

        print 'Starting to plot pure isothermal derivative properties'
        title = str('Derivative Properties\n%s' %(' + '.join(print_options[1])))
        figname = str('Deriv_Prop_%s.png' %('_'.join(print_options[1])))
        derivativeprop.plot_isothermal_dev_prop_pure(dp_dat[0],dp_dat[1],dp_dat[2],dp_dat[3],dp_dat[4],dp_dat[5],dp_dat[6],
                                                     dp_dat[7],dp_dat[8],dp_dat[9],dp_dat[10],dp_dat[11],dp_dat[12],
                                                     dp_dat[13],print_options,figname)
        print ('Figure %s saved successfully' %figname)
        #*******************************************************************************
        
    if env_type==2 or env_type==3 or env_type==4:
        #Calculate pure PV envelope*****************************************************
        nd = 500
        nx = 200
        n = 5
        finalT = 540.0
        stepT = 2.5
        print '\nCalculating PV envelope'
        if env_type==2:
            env_PV = PV_envelope(EoS,IDs,MR,T,finalT,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n)
        if env_type==3:
            env_PV = PV_estimate_Tc_envelope(EoS,IDs,MR,T,finalT,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n)
        if env_type==4:
            env_PV = PV_findTc3_envelope(EoS,IDs,MR,T,finalT,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,False,0,0)
            #env_PV = PV_findTc2_envelope(EoS,IDs,MR,T,finalT,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
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

    if env_type==6:
        #Calculate Pxy envelope*********************************************************

        if EoS==2 or EoS==4 or EoS==6:
            print '\nCalculating renormalized helmholtz energy surface'
            nd = 500
            nx = 200
            n = 5
            r_data = renormalization.renorm(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
            print '\nHelmholtz Energy Surface calculated and reported'
        else:
            r_data = []

        print '\nCalculating Pxy envelope'
        env_pxy = Pxy_envelope(T,IDs,EoS,MR,kij,nc,AR,CR,SM,r_data)
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

    if env_type==7:
        #Calculate critical point for mixture*******************************************************
        if EoS==2 or EoS==4 or EoS==6:
            print '\nCalculating renormalized helmholtz energy surface'
            nd = 40
            nx = 40
            n = 5
            #r_data = renormalization.renorm(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
            print '\nHelmholtz Energy Surface calculated and reported'
        else:
            nd = 40
            nx = 40
            n = 5
            r_data = []
        
        x = np.array([0.1,0.9])
        renormalization.crit_mix(False,0,0,EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,x)
        #critical.sadus_hicks_young(r_data[4],r_data[3],r_data[1],r_data[9])
        #*******************************************************************************

    if env_type==8 or env_type==9 or env_type==10:
        #Estimate Renormalization Parameters********************************************
        if env_type==8:
            nd = 400
            nx = 200
            n = 5
            finalT = 530.0
            stepT = 0.5
            expname = []
            expname.append('sat_exp_crit.data')

            print '\nEstimating renormalization Parameters'
            param = estimation.Estimate_Parameters_crit(EoS,IDs,MR,T,finalT,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,expname,True,False,AR)
            print 'Renormalization Parameters Estimated'

        #Estimate Association Parameters************************************************
        if env_type==9:
            nd = 200
            nx = 200
            n = 5
            finalT = 530.0
            stepT = 0.5
            expname = []
            expname.append('sat_exp.data')

            print '\nEstimating CPA Association Parameters'
            param = estimation.Estimate_Parameters_assoc(EoS,IDs,MR,T,finalT,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,expname,True,False,AR)
            print 'CPA Association Parameters Estimated'

        #Estimate CPA Parameters************************************************
        if env_type==10:
            nd = 200
            nx = 200
            n = 5
            finalT = 530.0
            stepT = 0.5
            expname = []
            expname.append('sat_exp.data')

            print '\nEstimating CPA Parameters'
            param = estimation.Estimate_Parameters_CPA(EoS,IDs,MR,T,finalT,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,expname,True,False,AR)
            print 'CPA Parameters Estimated'

            print 'Creating estimation report'
            #reportname = str('../output/OBJ_CPA_%s.csv' %('_'.join(print_options[1])))
            #report_param(param[2],reportname,print_options)
            #reportname = str('../output/param_CPA_%s.csv' %('_'.join(print_options[1])))
            #report_param(param[1],reportname,print_options)
            print ('Estimation reports saved successfully')
        
    #Analyze Entropy*****************************************************
    if env_type==11:
        nd = 500
        nx = 200
        n = 5
        stepT = 0.5
        finalT = T
        S_renorm = []
        print '\nCalculating entropy curve'
        S_renorm = PV_entropy(EoS,IDs,MR,T,finalT,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,P,1) #Isothermal
        #dp_dat = PV_deriv_calc_envelope(EoS,IDs,MR,T,finalT,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,P,2) #Isobaric
        print 'Entropy curve calculated'

        print 'Creating entropy report'
        reportname = str('Entropy_%s.csv' %('_'.join(print_options[1])))
        renormalization.report_entropy(reportname,S_renorm[0],S_renorm[1],S_renorm[2],print_options)
        print ('Report %s saved successfully' %reportname)

        print 'Starting to plot entropy curve'
        title = str('Entropy\n%s' %(' + '.join(print_options[1])))
        figname = str('Entropy_%s.png' %('_'.join(print_options[1])))
        renormalization.plot_entropy(S_renorm[0],S_renorm[1],S_renorm[2],print_options,figname)
        print ('Figure %s saved successfully' %figname)
        #*******************************************************************************
#======================================================================================
