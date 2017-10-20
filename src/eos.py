#File dedicated to functions for equation of state related functions
import data
import numpy as np
import math
import association
import menus
import numerical
import renormalization

#Ideal gas constant
R = 8.314462175e-6 #m3.MPa/K/mol

#EoS parameters-----------------------------------------------
def parameters(EoS):
    param = []
    
    sig = {
        1: 1, #SRK
        2: 1, #SRK+RG
        3: 1+2**(0.5), #PR
        4: 1+2**(0.5), #PR+RG
        5: 1, #CPA
        6: 1  #CPA+RG
    }.get(EoS,'NULL')
    
    eps = {
        1: 0, #SRK
        2: 0, #SRK+RG
        3: 1-2**(0.5), #PR
        4: 1-2**(0.5), #PR+RG
        5: 0, #CPA
        6: 0, #CPA+RG
    }.get(EoS,'NULL')
    
    OMEGAa = {
        1: 0.42748, #SRK
        2: 0.42748, #SRK+RG
        3: 0.45724, #PR
        4: 0.45724, #PR+RG
        5: 0.42748, #CPA
        6: 0.42748  #CPA+RG
    }.get(EoS,'NULL')
    
    OMEGAb = {
        1: 0.08664, #SRK
        2: 0.08664, #SRK+RG
        3: 0.07780, #PR
        4: 0.07780, #PR+RG
        5: 0.08664, #CPA
        6: 0.08664  #CPA+RG
    }.get(EoS,'NULL')
    
    param = [sig,eps,OMEGAa,OMEGAb]
    
    return param
#=============================================================

#Cubic, calculates ac parameter-------------------------------
def ac_calc(IDs,EoS,T):
    
    Tc = np.array(data.Tc(IDs))
    Pc = np.array(data.Pc(IDs))
    alfa = np.array(alfa_calc(IDs,EoS,T))
    param = parameters(EoS)
    OMEGAa = param[2]
    
    #Calculate ac
    ac = {
        1: OMEGAa*R*R*Tc*Tc/Pc, #SRK
        2: OMEGAa*R*R*Tc*Tc/Pc, #SRK+RG
        3: OMEGAa*R*R*Tc*Tc/Pc, #PR
        4: OMEGAa*R*R*Tc*Tc/Pc, #PR+RG
        5: np.array(data.a0(IDs)), #CPA
        6: np.array(data.a0(IDs))  #CPA+RG
    }.get(EoS,'NULL')
    
    return ac
#=============================================================

#Cubic, calculates a parameter--------------------------------
def a_calc(IDs,EoS,T):
    Tc = np.array(data.Tc(IDs))
    ac = np.array(ac_calc(IDs,EoS,T))
    alfa = np.array(alfa_calc(IDs,EoS,T))
    
    #Calculate a
    a = {
        1: ac*alfa, #SRK
        2: ac*alfa, #SRK+RG
        3: ac*alfa, #PR
        4: ac*alfa, #PR+RG
        5: ac*alfa, #CPA
        6: ac*alfa #CPA+RG
    }.get(EoS,'NULL')
    
    return a
#=============================================================

#Cubic, calculates a prime parameter--------------------------
def a_pr_calc(IDs,EoS,T):
    Tc = np.array(data.Tc(IDs))
    ac = np.array(ac_calc(IDs,EoS,T))
    omega = np.array(data.omega(IDs))
    
    #Calculate kapa
    kapa = {
        1: 0.48508+1.55171*omega-0.15613*omega*omega,#SRK
        2: 0.48508+1.55171*omega-0.15613*omega*omega, #SRK+RG
        3: 0.37464+1.54226*omega-0.26992*omega*omega, #PR
        4: 0.37464+1.54226*omega-0.26992*omega*omega #PR+RG
    }.get(EoS,'NULL')
    
    #Calculate a
    a = {
        1: -ac*kapa*((1+kapa)/np.sqrt(T/Tc)-kapa)/Tc,#SRK
        2: -ac*kapa*((1+kapa)/np.sqrt(T/Tc)-kapa)/Tc, #SRK+RG
        3: -ac*kapa*((1+kapa)/np.sqrt(T/Tc)-kapa)/Tc, #PR
        4: -ac*kapa*((1+kapa)/np.sqrt(T/Tc)-kapa)/Tc, #PR+RG
    }.get(EoS,'NULL')
    
    return a
#=============================================================

#Cubic, calculates b parameter--------------------------------
def b_calc(IDs,EoS):
    Tc = np.array(data.Tc(IDs))
    Pc = np.array(data.Pc(IDs))
    param = parameters(EoS)
    OMEGAb = param[3]
    
    #Calculate b
    b = {
        1: OMEGAb*R*Tc/Pc,#SRK
        2: OMEGAb*R*Tc/Pc, #SRK+RG
        3: OMEGAb*R*Tc/Pc, #PR
        4: OMEGAb*R*Tc/Pc, #PR+RG
        5: np.array(data.bCPA(IDs)), #CPA
        6: np.array(data.bCPA(IDs))  #CPA+RG
    }.get(EoS,'NULL')
    
    return b
#=============================================================

#Cubic, calculates alfa parameter-----------------------------
def alfa_calc(IDs,EoS,T):
    omega = np.array(data.omega(IDs))
    Tc = np.array(data.Tc(IDs))
    Tr = T/Tc
    
    #Calculate kapa
    if omega[0]<0.2:
        kapa = {
            1: 0.48508+1.55171*omega-0.15613*omega*omega, #SRK
            2: 0.48508+1.55171*omega-0.15613*omega*omega, #SRK+RG
            3: 0.37640+1.54230*omega-0.26990*omega*omega, #PR
            4: 0.37640+1.54230*omega-0.26990*omega*omega, #PR+RG
            5: np.array(data.c1(IDs)), #CPA
            6: np.array(data.c1(IDs))  #CPA+RG
        }.get(EoS,'NULL')
    
    else:
        kapa = {
            1: 0.48508+1.55171*omega-0.15613*omega*omega, #SRK
            2: 0.48508+1.55171*omega-0.15613*omega*omega, #SRK+RG
            3: 0.3796+(1.4850+(-0.1644+0.01667*omega)*omega)*omega, #PR
            4: 0.3796+(1.4850+(-0.1644+0.01667*omega)*omega)*omega, #PR+RG
            5: np.array(data.c1(IDs)), #CPA
            6: np.array(data.c1(IDs))  #CPA+RG
        }.get(EoS,'NULL')
        
    alfa = (1.0+kapa*(1.0-(Tr**0.5)))**2
    
    return alfa
#=============================================================

#Cubic, calculates alfa prime parameter-----------------------
def alfa_pr_calc(IDs,EoS,T):
    omega = np.array(data.omega(IDs))
    Tc = np.array(data.Tc(IDs))
    Tr = T/Tc
    
    #Calculate kapa
    kapa = {
        1: 0.48508+1.55171*omega-0.15613*omega*omega, #SRK
        2: 0.48508+1.55171*omega-0.15613*omega*omega, #SRK+RG
        3: 0.37464+1.54226*omega-0.26992*omega*omega, #PR
        4: 0.37464+1.54226*omega-0.26992*omega*omega  #PR+RG
    }.get(EoS,'NULL')
    
    alfa_pr = kapa*(kapa-(1+kapa)*(Tr**(-0.5)))
    return alfa_pr
#=============================================================

#Cubic, calculates 'a' mixture parameter----------------------
def amix_calc(MR,a,x,kij):
    a = np.array(a)
    x = np.array(x)
    nc = x.shape[0]
    
    xi = np.zeros((nc,nc))
    xj = np.zeros((nc,nc))
    ai = np.zeros((nc,nc))
    aj = np.zeros((nc,nc))
    aij = np.zeros((nc,nc))
    
    for i in range(0,nc):
        ai[:,i] = a
        xi[:,i] = x
    
    aj = np.transpose(ai)
    xj = np.transpose(xi)
    aij = (np.sqrt(ai*aj))*(1-kij)
    
    #Calculate bmix
    amix = {
        1: np.sum(xi*xj*aij) #VdW1f
    }.get(MR,'NULL')
    
    return amix
#=============================================================

#Cubic, calculates 'b' mixture parameter----------------------
def bmix_calc(MR,b,x):
    b = np.array(b)
    x = np.array(x)
    
    #Calculate bmix
    bmix = {
        1: np.sum(x*b) #VdW1f
    }.get(MR,'NULL')
    return bmix
#=============================================================

#Cubic, molar volume roots calculation------------------------
def V_calc(EoS,P,T,amix,bmix,phase):

    #EoS parameters
    param = parameters(EoS)
    sig = param[0]
    eps = param[1]
    
    #Cubic coefficients
    alfa = (sig+eps-1)*bmix-R*T/P
    beta = sig*eps*bmix*bmix-(R*T/P+bmix)*(sig+eps)*bmix+amix/P
    gama = -(R*T/P+bmix)*sig*eps*bmix*bmix-bmix*amix/P
    
    coef = [alfa,beta,gama]
    
    #Find volume roots
    V = numerical.cubic_v_roots(coef,R,T,P)
    
    #phase 1 indicates vapor phase, -1 indicates liquid phase
    if phase==1:
        V = np.amax(V)
    else:
        V = np.amin(V)
        
    return V
#=============================================================

#CPA, objective function to calculate molar volume------------
def V_objf(nc,V,R,P,amix,bmix,T,x,X,Bcpa,iota,i,a):
    V = bmix/iota;
    one4 = np.ones((4))
    x4 = np.kron(x,one4)
    one4nc = np.ones((4*nc))

    pre_F_v = np.sum(x4*(one4nc-X))
    F = ((R*T/(V-bmix) - amix/(V*(V+bmix)) - 0.5*R*T/V * (1+0.475*Bcpa/(V-0.475*Bcpa)) * pre_F_v) - P)
    F = (1-iota)*(F)
    return F
#=============================================================

#CPA, molar volume and X calculation--------------------------
def V_CPA(EoS,P,T,amix,bmix,b,phase,Vinit,CR,en_auto,beta_auto,x,a,SM):

    x = np.array(x)
    nc = x.shape[0]
    
    if nc==1:
        nc=2
        
    one4nc = np.ones((4*nc))
    X = np.ones((4*nc))
    k = 0
    deltaV = 0
    Bcpa = np.sum(x*b);
    
    #Initial guess for V comes from SRK
    #Initial guess for iota
    if phase==1: #Vapor phase
        iota = bmix/Vinit #If first iteration, iota = 0.99
    else:
        iota = bmix/Vinit #If first iteration, iota = bmix/(bmix+(R*T/P))
        
    V = bmix/iota
    
    #Definitions
    iota_min = 0
    iota_max = 1
    i=0
    iteri = 0
    tolV = 1e-6 #original 1e-6
    cond_iota = tolV+1
    deltaV = 0
    
    while cond_iota>tolV:
        if i==0:
            Xf = association.frac_nbs(nc,V,CR,en_auto,beta_auto,b,bmix,X,i,x,deltaV,T,SM)
            X = Xf
        
        F_obj = V_objf(nc,V,R,P,amix,bmix,T,x,X,Bcpa,iota,i,a);
        F_obj_plus = V_objf(nc,V,R,P,amix,bmix,T,x,X,Bcpa,iota+iota*1e-6,i,a);
        F_obj_minus = V_objf(nc,V,R,P,amix,bmix,T,x,X,Bcpa,iota-iota*1e-6,i,a);
        F_obj_derivative =(F_obj_plus-F_obj_minus)/(2*iota*1e-6);
        
        if i==0:
            F_obj_derivative = 1e6
        
        iota_old = iota
        F_obj_old = F_obj
        
        if F_obj>0:
            iota_max = iota
        else:
            iota_min = iota

        div = F_obj/F_obj_derivative
        iota_new = iota - div

        if iota_min<iota_new and iota_new<iota_max:
            iota2 = iota_new
        else:
            iota2 = (iota_min+iota_max)/2

        cond_iota = abs(iota2 - iota)/iota
        cond_iota = abs(F_obj)

        if i>0:
            cond_iota = abs(F_obj/F_obj_derivative)

        i = i+1
        iota = iota2
        deltaV = bmix/iota - V
        V = bmix/iota
        V_obj = F_obj

        if i!=0:
            Xf = association.frac_nbs(nc,V,CR,en_auto,beta_auto,b,bmix,X,i,x,deltaV,T,SM)
            X = Xf
            
        if i>500:
            cond_iota = tolV-1
            V = Vinit
            print 'VOLUME MAX ITER REACHED \n'
            i = 0
        
        if math.isnan(V) or math.isinf(V):
            V = Vinit+Vinit*0.0001*k;
            iota = bmix/V;
            Xf = association.frac_nbs(nc,V,CR,en_auto,beta_auto,b,bmix,X,i,x,deltaV,T,SM)
            X = Xf
            cond_iota = tolV+1 #original is tolV+1
            print 'VOLUME = NAN or INF',k
            Vinit = V
            cond_iota = tolV-1
            k = k+1
    
    out = []
    out.append(V)
    out.append(X)
    return out
#=============================================================

#Molar volume calculation-------------------------------------
def V_func(EoS,P,T,amix,bmix,b,phase,Vinit,CR,en_auto,beta_auto,x,a,SM,r_data):

    if EoS==6:
        V = {
            1: V_calc(EoS,P,T,amix,bmix,phase),
            2: V_calc(EoS,P,T,amix,bmix,phase),
            3: V_calc(EoS,P,T,amix,bmix,phase),
            4: V_calc(EoS,P,T,amix,bmix,phase),
            5: V_CPA(EoS,P,T,amix,bmix,b,phase,Vinit,CR,en_auto,beta_auto,x,a,SM),
            6: renormalization.volume_renorm(phase,x[0],P,bmix,R,T,r_data)
        }.get(EoS,'NULL')
    else:
        V = {
            1: V_calc(EoS,P,T,amix,bmix,phase),
            2: V_calc(EoS,P,T,amix,bmix,phase),
            3: V_calc(EoS,P,T,amix,bmix,phase),
            4: V_calc(EoS,P,T,amix,bmix,phase),
            5: V_CPA(EoS,P,T,amix,bmix,b,phase,Vinit,CR,en_auto,beta_auto,x,a,SM)
        }.get(EoS,'NULL')

    return V
#=============================================================

#Cubic, residual Enthalpy calculation-------------------------
def resH_calc(EoS,P,T,V,a,a_pr,b):
    
    #EoS parameters
    param = parameters(EoS)
    sig = param[0]
    eps = param[1]
    
    PHI = 1/(b*(sig-eps))*np.log((V+sig*b)/(V+eps*b))
    
    Hres = P*V-R*T+PHI*(T*a_pr-a)

    return Hres
#=============================================================

#Cubic, residual Entropy calculation--------------------------
def resS_calc(EoS,P,T,V,a,a_pr,b):
    
    #EoS parameters
    param = parameters(EoS)
    sig = param[0]
    eps = param[1]
    
    PHI = 1/(b*(sig-eps))*np.log((V+sig*b)/(V+eps*b))
    
    Sres = R*np.log(P*(V-b)/(R*T))+a_pr*PHI

    return Sres
#=============================================================

#Cubic, residual Gibbs energy calculation---------------------
def resG_calc(EoS,P,T,V,a,a_pr,b):
    
    Hres = resH_calc(EoS,P,T,V,a,a_pr,b)
    Sres = resS_calc(EoS,P,T,V,a,a_pr,b)
    Gres = Hres-T*Sres

    return Gres
#=============================================================

#Cubic, residual internal energy calculation------------------
def resU_calc(EoS,P,T,V,a,a_pr,b):
    
    Hres = resH_calc(EoS,P,T,V,a,a_pr,b)
    Ures = Hres-P*V

    return Ures
#=============================================================

#FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Cubic, residual Gibbs energy calculation---------------------
def resA_calc(EoS,IDs,P,T,V,a,a_pr,b,x,kij):
    
    #Ures = resU_calc(EoS,P,T,V,a,a_pr,b)
    #Sres = resS_calc(EoS,P,T,V,a,a_pr,b)
    #Ares = Ures-T*Sres
    
    a1 = a_calc(IDs,EoS,T)
    b1 = b_calc(IDs,EoS)
    
    x = np.array(x)
    nc = x.shape[0]
    ai = np.zeros((nc,nc))
    aj = np.zeros((nc,nc))
    aij = np.zeros((nc,nc))
    bi = np.zeros((nc,nc))
    bj = np.zeros((nc,nc))
    bij = np.zeros((nc,nc))
    
    for i in range(0,nc):
        ai[:,i] = a1
        bi[:,i] = b1
    
    aj = np.transpose(ai)
    bj = np.transpose(bi)
    bij = (bi+bj)/2
    aij = (np.sqrt(ai*aj))*(1-kij)
    
    n = np.sum(x)
    Vt = V*n
    Biv1 = 2*np.dot(bij,x)
    Biv = (Biv1-b)/n
    Div = 2*np.dot(aij,x)
    Fn = -(np.log(1-b/V));
    f = (np.log(1+b/Vt))/(R*b)
    fV = -(1/(R*Vt*(Vt+b)))
    D_T = a
    FB = n/(Vt-b)+D_T*(f+Vt*fV)/(T*b)
    FD = -f/T
    Ares = Fn + FB*Biv + FD*Div

    return Ares
#=============================================================

#Cubic, ln fugacity of component in a mixture calculation-----
def lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,phase):
    
    #EoS parameters
    param = parameters(EoS)
    sig = param[0]
    eps = param[1]
    
    x = x/np.sum(x)
    
    a = a_calc(IDs,EoS,T)
    b = b_calc(IDs,EoS)
    
    amix = amix_calc(MR,a,x,kij)
    bmix = bmix_calc(MR,b,x)
    
    V = V_calc(EoS,P,T,amix,bmix,phase)
    Z = P*V/(R*T)
    
    C_ = {
        1: -0.69314, #SRK
        2: -0.69314, #SRK+RG
        3: -0.62323, #PR
        4: -0.62323  #PR+RG
    }.get(EoS,'NULL')
    
    B = bmix*P/(R*T);
    qe = (amix/(bmix*R*T));
    
    q_ = {
        1: qe*(b/bmix-2*np.sqrt(a/amix)), #VdW1f
        #2: -((a/b)/(R*T)+lngamma/C_) #HV
    }.get(MR,'NULL')
    
    I = (1/(sig-eps))*(np.log((Z+B*sig)/(Z+B*eps)))
    logZB = np.log(Z-B)
    q_I = q_*I
    Z1 = Z-1
    bZ1 = b/bmix*Z1
    lnfugcoef = bZ1-logZB+q_I
    
    return lnfugcoef
#=============================================================

#CPA, ln fugacity of component in a mixture calculation-------
def lnfugcoef_CPA(IDs,EoS,MR,P,T,x,kij,phase,V,X):
    
    x = np.array(x)
    nc = np.shape(x)[0]
    
    x = x/np.sum(x) #EXTRA
    
    if nc==1:
        nc=2
    
    one4 = np.ones((4))
    one4c = np.ones((4*nc,nc))
    one4nc = np.ones((4*nc))
    onenc = np.ones((nc))
    
    for i in range(0,4*nc):
        for j in range(0,nc):
            if math.trunc(i/4) == j:
                one4c[i][j] = 1
            else:
                one4c[i][j] = 0
    
    a = a_calc(IDs,EoS,T)
    b = b_calc(IDs,EoS)
    
    amix = amix_calc(MR,a,x,kij)
    bmix = bmix_calc(MR,b,x)
    
    bi = np.zeros((nc,nc))
    bij = np.zeros((nc,nc))
    ai = np.zeros((nc,nc))
    aij = np.zeros((nc,nc))
    for i in range(0,nc):
        bi[:,i] = b
        ai[:,i] = a
    bij = (bi+bi.T)/2
    aij = (np.sqrt(ai*(ai.T)))*(1-kij)
    
    n = np.sum(x)
    
    Bcpa = np.dot(x,b)
    Biv = (2*(np.dot(bij,x))-Bcpa)/n
    Div = 2*(np.dot(aij,x))
    
    Vt = V*n
    
    dlngdn = b*(0.475/(V-0.475*Bcpa))
    p_dlng_dp = 0.475*bmix/(V-0.475*bmix)
    eta = bmix/(4*V)
    g = 1/(1-1.9*eta)
    x4 = np.kron(x,one4)
    h = np.sum(x4*(one4nc-X))

    #Associative compressibility factor
    Z_assoc = -0.5*(1+p_dlng_dp)*h;

    #Physical(SRK) compressibility factor
    Z_SRK = V/(V-bmix)-amix/(R*T*(V+bmix))

    #CPA compressibility factor
    Z_CPA = Z_SRK + Z_assoc;

    #Association residual chemical potential
    lnX = np.log(X)
    u_assoc_1 = np.dot(one4nc,np.dot(np.diag(lnX),one4c))
    u_assoc_2 = h*dlngdn
    u_assoc = u_assoc_1-0.5*u_assoc_2
    
    #Physical(SRK) residual chemical potential
    Fn = -(np.log(1-Bcpa/Vt))
    f = (np.log(1+Bcpa/Vt))/(R*Bcpa)
    fV = -(1/(R*Vt*(Vt+Bcpa)))
    D_T = amix
    FB = n/(Vt-Bcpa)+D_T*(f+Vt*fV)/(T*Bcpa)
    FD = -f/T
    u_SRK = Fn + FB*Biv + FD*Div
    
    #CPA residual chemical potential
    u_CPA = u_SRK + u_assoc
    
    lnfugcoef = u_CPA - np.log(Z_CPA)
    #print 'ures=',u_CPA*R*T,Z_CPA,R,T
    
    return lnfugcoef
#=============================================================

#ln fugacity of component in a mixture calculation------------
def lnfugcoef_func(IDs,EoS,MR,P,T,x,kij,phase,V,en_auto,beta_auto,CR,SM,it,pt,r_data):
    
    X = np.empty(4)
    if EoS==5 or EoS==6:
        a = a_calc(IDs,EoS,T)
        b = b_calc(IDs,EoS)
        amix = amix_calc(MR,a,x,kij)
        bmix = bmix_calc(MR,b,x)
        if it<1 or it>=1:
            if phase==1:
                V = bmix+(R*T/P)
            else:
                V = bmix/0.99
        cpa = V_func(EoS,P,T,amix,bmix,b,phase,V,CR,en_auto,beta_auto,x,a,SM,r_data)
        V = cpa[0]
        X = cpa[1]
    
    out = []
    
    if EoS==1:
        lnfugcoef = lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,phase)
        out.append(lnfugcoef)
    elif EoS==2:
        lnfugcoef = lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,phase)
        out.append(lnfugcoef)
    elif EoS==3:
        lnfugcoef = lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,phase)
        out.append(lnfugcoef)
    elif EoS==4:
        lnfugcoef = lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,phase)
        out.append(lnfugcoef)
    elif EoS==5:
        lnfugcoef = lnfugcoef_CPA(IDs,EoS,MR,P,T,x,kij,phase,V,X)
        out.append(lnfugcoef)
        out.append(V)
    elif EoS==6:
        lnfugcoef = renormalization.lnfugcoef_renorm(P,T,x,phase,V,r_data,bmix)
        out.append(lnfugcoef)
        out.append(V)
    else:
        out.append('NULL')
    
    #print x,lnfugcoef,(1/V)*bmix,V,phase
    #raw_input('...')
    return out
#=============================================================

#ln fugacity of component in a mixture calculation------------
def lnfugcoef_funcFAKE(IDs,EoS,MR,P,T,x,kij,phase,V,en_auto,beta_auto,CR,SM,it,pt):
    
    X = np.empty(4)
    if EoS==5:
        a = a_calc(IDs,EoS,T)
        b = b_calc(IDs,EoS)
        amix = amix_calc(MR,a,x,kij)
        bmix = bmix_calc(MR,b,x)
        if it<1 or it>=1:
            if phase==1:
                V = bmix+(R*T/P)
            else:
                V = bmix/0.99
        cpa = V_func(EoS,P,T,amix,bmix,b,phase,V,CR,en_auto,beta_auto,x,a,SM)
        V = cpa[0]
        X = cpa[1]
    
    out = []
    
    if EoS<5:
        lnfugcoef = {
            1: lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,phase),
            2: lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,phase),
            3: lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,phase),
            4: lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,phase)
        }.get(EoS,'NULL')
        out.append(lnfugcoef)
    
    else:
        lnfugcoef = {
            1: lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,phase),
            2: lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,phase),
            3: lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,phase),
            4: lnfugcoef_calc(IDs,EoS,MR,P,T,x,kij,phase),
            5: lnfugcoef_CPA(IDs,EoS,MR,P,T,x,kij,phase,V,X)
        }.get(EoS,'NULL')
        out.append(lnfugcoef)
        out.append(V)
    
    return out
#=============================================================
