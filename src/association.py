import numpy as np
import scipy
import data
import menus

R = 8.314462175e-6 #m3.MPa/K/mol

#Build auto association matrix based on AR vector------------------------------------------
def CPA_auto(AR,nc,IDs):

    if nc==1:
        nc =2
    
    enval = data.en(IDs)
    betaval = data.beta(IDs)
    
    en = np.zeros((4*nc,4*nc))
    beta = np.zeros((4*nc,4*nc))
    enij = np.zeros((4,4))
    betaij = np.zeros((4,4))
    
    for i in range(0,nc):
        enij = auto_assoc(enval[i],AR[i])
        betaij = auto_assoc(betaval[i],AR[i])
        if i==0:
            en = enij
            beta = betaij
        else:
            en = np.hstack((en,enij))
            beta = np.hstack((beta,betaij))
            
    en1 = en
    beta1 = beta
    for i in range(0,nc-1):
        en = np.vstack((en,en1))
        beta = np.vstack((beta,beta1))
    
    auto = []
    auto.append(en)
    auto.append(beta)
    return auto
#==========================================================================================

#Gives back auto association matrix with defined association rule--------------------------
def auto_assoc(beta,AR):
    A1 =    np.array([[beta,    0,    0,    0],
                      [   0,    0,    0,    0],
                      [   0,    0,    0,    0],
                      [   0,    0,    0,    0]])

    A2 =    np.array([[beta, beta,    0,    0],
                      [beta, beta,    0,    0],
                      [   0,    0,    0,    0],
                      [   0,    0,    0,    0]])
                  
    B2 =    np.array([[   0, beta,    0,    0],
                      [beta,    0,    0,    0],
                      [   0,    0,    0,    0],
                      [   0,    0,    0,    0]])
                  
    A3 =    np.array([[beta, beta, beta,    0],
                      [beta, beta, beta,    0],
                      [beta, beta, beta,    0],
                      [   0,    0,    0,    0]])
                  
    B3 =    np.array([[   0,    0, beta,    0],
                      [   0,    0, beta,    0],
                      [beta, beta,    0,    0],
                      [   0,    0,    0,    0]])

    A4 =    np.array([[beta, beta, beta, beta],
                      [beta, beta, beta, beta],
                      [beta, beta, beta, beta],
                      [beta, beta, beta, beta]])

    B4 =    np.array([[   0,    0,    0, beta],
                      [   0,    0,    0, beta],
                      [   0,    0,    0, beta],
                      [beta, beta, beta,    0]])

    C4 =    np.array([[   0, beta, beta,    0],
                      [beta,    0,    0, beta],
                      [beta,    0,    0, beta],
                      [   0, beta, beta,    0]])

    auto = {
        1: A1,
        2: A2,
        3: B2,
        4: A3,
        5: B3,
        6: A4,
        7: B4,
        8: C4
    }.get(AR,'NULL')
    
    return auto
#==========================================================================================

#Calculates cross-association Delta--------------------------------------------------------
def delta_calc(nc,V,CR,en_auto,beta_auto,b,bmix,T,SM):
    
    if nc==1:
        nc =2
    
    b = np.array(b)
    
    #Reduced density
    eta = bmix/(4*V)
    
    #Radial distribution function (sCPA)
    g = 1/(1-1.9*eta)
    
    #bij
    bi = np.zeros((nc,nc))
    bij = np.zeros((nc,nc))
    one4 = np.ones((4,4))
    for i in range(0,nc):
        bi[:,i] = b
    bij = (bi+bi.T)/2
    Bij = np.kron(bij,one4)
    
    #cross-association arrays
    en_cross = {
        1: (en_auto+en_auto.T)/2, #CR-1
        2: (en_auto+en_auto.T)/2, #CR-2
        3: np.sqrt(en_auto*en_auto.T), #CR-3
        4: np.sqrt(en_auto*en_auto.T), #CR-4
    }.get(CR,'NULL')
    
    beta_cross = {
        1: np.sqrt(beta_auto*beta_auto.T), #CR-1
        2: (beta_auto+beta_auto.T)/2     , #CR-2
        3: np.sqrt(beta_auto*beta_auto.T), #CR-3
        4: (beta_auto+beta_auto.T)/2     , #CR-4
    }.get(CR,'NULL')
    
    #Apply solvation and cross-association considerations
    one4 = np.ones((4,4))
    np.fill_diagonal(SM,1) #Make sure auto-associations are preserved
    SM = np.kron(SM,one4)
    en_cross = en_cross*SM
    beta_cross = beta_cross*SM
    
    if CR!=5:
        delta = g*(np.exp(en_cross/(R*T))-1)*Bij*beta_cross
    else:
        delta_auto = g*(np.exp(en_auto/(R*T))-1)*Bij*beta_auto
        delta = np.sqrt(delta_auto*delta_auto.T)
    
    return delta
#==========================================================================================

#Calculates Fraction of non-bonded sites---------------------------------------------------
def frac_nbs(nc,V,CR,en_auto,beta_auto,b,bmix,X,Viter,x,deltaV,T,SM):
    
    if nc==1:
        nc=2
    
    #Calculate delta
    delta = delta_calc(nc,V,CR,en_auto,beta_auto,b,bmix,T,SM)
    
    x = np.array(x)
    X = np.array(X)
    
    #auxiliary arrays
    one4 = np.ones((4))
    one44 = np.ones((4,4))
    one4nc = np.ones((4*nc))
    x4 = np.kron(x.T,one4)
    Xx4 = X*x4
    Xx4m = np.zeros((nc*4,nc*4))
    x4m = np.zeros((nc*4,nc*4))
    for i in range(0,nc*4):
        Xx4m[:,i] = Xx4
        x4m[:,i]  = x4

    #Q derivatives and X adjustment with volume
    xXD = Xx4m*delta
    pre_Xnew = np.dot(one4nc.T,xXD)/V
    QXV = x4*pre_Xnew*(1/V-1)
    X2 = X**2
    K = x4m*x4m.T*delta/V
    QXX1 = np.diag(x4/(X**2))
    QXX = -QXX1-K
    dXdV = np.dot((np.linalg.inv(QXX)),(-QXV))
    X = X + dXdV*deltaV

    #If first iteration on volume, X have all values set to 0.2
    if Viter==0:
        for i in range(0,nc*4):
            X[i] = 0.2
            
    #Main iteration definitions
    k = 1                   #Iteration counter
    tolX = 1e-7             #Tolerace
    Xmax = tolX+1           #Force condition to enter
    I = np.identity(4*nc)   #Identity matrix
    gmax = tolX+1
    
    while (Xmax>tolX or gmax>tolX):
        #First steps using sucessive substitution
        if k<=4: 
            delta = delta_calc(nc,V,CR,en_auto,beta_auto,b,bmix,T,SM)
            Xx4 = X*x4
            for i in range(0,nc*4):
                Xx4m[:,i] = Xx4
            xXD = Xx4m*delta
            pre_Xnew = np.dot(one4nc,xXD)/V
            Xnew = 0.8*(1/(pre_Xnew+one4nc))+0.2*X
            Xcond = np.absolute(Xnew-X)
            Xmax = np.amax(Xcond)
            Xmax = Xmax/np.amax(X)
            X = Xnew
        
        #Other steps using Hessian Matrix
        else:
            delta = delta_calc(nc,V,CR,en_auto,beta_auto,b,bmix,T,SM)
            Xx4 = X*x4
            for i in range(0,nc*4):
                Xx4m[:,i] = Xx4
            K = x4m*x4m.T*delta/V
            H_1 = np.diag((x4+np.dot(X,K))/X)
            H = -H_1-K
            g = np.dot((np.diag(1/X)-I),x4.T).T-np.dot(X,K)
            gmax = np.amax(g)

            deltaX = np.dot(np.linalg.inv(H),(-g))
            
            X1 = np.log(X)-X+1
            
            Qold = np.dot(X,(np.dot(X,K).T))
            Qold = np.sum(x4*X1)-0.5*Qold
            Q = Qold-1 #Force enter loop

            while Q<Qold:
                Xnew = X+deltaX
                for l in range (0,4*nc):
                    Xnew[l] = np.maximum(Xnew[l],0.2*X[l])
                X1 = np.log(Xnew)-Xnew+1
                Q = np.dot(Xnew,(np.dot(Xnew,K).T))
                Q = np.sum(x4*X1)-0.5*Q
                if Q<Qold:
                    deltaX = deltaX/2
            
            Xcond = np.absolute(Xnew-X)
            Xmax = np.amax(Xcond)
            X = Xnew
            
        k = k + 1
        
    return X
#==========================================================================================