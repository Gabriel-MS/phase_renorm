#############HOLDS OLDER OR NON-WORKING VERSIONS OF SOME FUNCTIONS###############

#Given pure component isotherm, calculates phase coexistence densities-----------------
def coexistence_dens_2(rho1,f1):
    
    rho1 = np.array(rho1)
    f1 = np.array(f1)
    n = f1.shape[1]
    rho = rho1.flatten()
    f = f1.flatten()
    #Spline to get chemical potential
    fspl = splrep(rho,f,k=3)         #Cubic Spline Representation
    u = splev(rho,fspl,der=1)        #Evaluate Cubic Spline First derivative
    
    drho = rho[n/2]-rho[n/2-1]
    for i in range(1,n-2):
        u[i] = (f[i+1]-f[i-1])/(2*drho)
    u[n-1] = (f[n-1]-f[n-2])/drho
    u[0] = (f[1]-f[0])/drho
    
    #Calcule pressure
    P = -f+rho*u
    a = (P[n/2-15]-P[n/2+15])/(rho[n/2-15]-rho[n/2+15])
    b = P[n/2+15]-a*rho[n/2+15]
    i = n/2-15
    while i<(n/2+15):
        P[i]=a*rho[i]+b
        i = i + 1
        
    #plt.plot(rho,P)
    #plt.ylim(0,20)
    #plt.show()

    #First derivative
    #drho = rho[1]-rho[0]
    Pspl = splrep(rho,P,k=3,s=2)         #Cubic Spline Representation
    dPdrho = splev(rho,Pspl,der=1)        #Evaluate Cubic Spline First derivative
    
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
    
    #print '/---'
    #print rho[max2],rho[max1],rho[min1],rho[min2]
    #print Pmax,Pmin
    #print '-----/'
    #input ('...')

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

    #print 'rhoss',rho[0],rhomax,rhomin,rho[min2]
    #input('...')
    
    i = 0
    rho1=rho[0]-1
    while rho1<rho[0] and rho1<0:
        rho1 = Pf1roots[i]
        i = i+1
    if Pf1roots[0]>Pf2roots[0]:
        rho1 = 0.1
    rho1 = numerical.falsi_spline(rho,Pf1,rho[0],rhomax,1e-5)
    
    
    i = 0
    rho2 = rhomin-1
    while rho2<rhomin:
        rho2 = Pf2roots[i]
        i = i+1
    rho2 = numerical.falsi_spline(rho,Pf2,rhomin,rho[min2],1e-5)

    #print 'rho1,rho2',rho1,rho2
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
    #Pspl = splrep(rho,P,k=3)
    fspl1 = splev(rho,fspl)
    uspl1 = splev(rho,uspl)
    dudrho = splev(rho,uspl,der=1)
    dPdrho = splev(rho,Pspl,der=1)
    Nitmax = 10000
    Nit = 0    

    while (abs(du)>tol or abs(dP)>tol) and (abs(Pmax-Pmin)>1e-3) and (Nit<Nitmax):
    #while (abs(drho1)>tol or abs(drho2)>tol) and (abs(Pmax-Pmin)>1e-3) and (Nit<Nitmax):
        Nit = Nit+1
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
        
        if drho1>200:
            drho1=200.0
        if drho1<-200:
            drho1=-200.0
        if drho2>200:
            drho2=200.0
        if drho2<-200:
            drho2=-200.0
        
        rho1 = rho1 + stop*drho1
        rho2 = rho2 + stop*drho2
        
        #if counter>0 and (drho1>drho1old and drho2>drho2old):
        #    rho1 = rho1 - stop*drho1
        #    rho2 = rho2 - stop*drho2
        #    stop = stop/1.05 #Break
        #    rho1 = rho1 + stop*drho1/2
        #    rho2 = rho2 + stop*drho2/2
        #    print 'stop',stop,counter
        
        du = abs(u1-u2)
        dP = abs(P1-P2)
        #print rho1,rho2,du,dP,drho1,drho2,stop
        
        counter = counter+1
        drho1old = drho1
        drho2old = drho2
        duold = du
        dPold = dP

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
    #input('...')
    
    return dens
#======================================================================================