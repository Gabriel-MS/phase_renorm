#############HOLDS OLDER OR NON-WORKING VERSIONS OR PIECES OF SOME FUNCTIONS###############

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
    T = 512.6
    phi = 1.0
    L = 4.0e-10
    step = tol+1
    i = 0
    
    print 'T:   dens_vap:   dens_liq:   P:'
    while math.fabs(step)>tol:
        ren = renormalization.renorm(EoS,IDs,MR,T,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
        dens = coexistence_dens(ren[2],ren[0])
        Tv.append(T)
        rhov.append(dens[0])
        rhol.append(dens[1])
        Pv.append(dens[2])
        Fobj = math.fabs(dens[0]-dens[1])
        #print 'central',T,dens[0],dens[1],dens[2],Fobj

        ren = renormalization.renorm_est(EoS,IDs,MR,T+1,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
        dens = coexistence_dens(ren[2],ren[0])
        Tv.append(T)
        rhov.append(dens[0])
        rhol.append(dens[1])
        Pv.append(dens[2])
        Fobj_plus = math.fabs(dens[0]-dens[1])
        #print 'plus',T,dens[0],dens[1],dens[2],Fobj_plus

        ren = renormalization.renorm_est(EoS,IDs,MR,T-1,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
        dens = coexistence_dens(ren[2],ren[0])
        Tv.append(T)
        rhov.append(dens[0])
        rhol.append(dens[1])
        Pv.append(dens[2])
        Fobj_minus = math.fabs(dens[0]-dens[1])
        #print 'minus',T,dens[0],dens[1],dens[2],Fobj_minus

        Fobj_derivative = (Fobj_plus-Fobj_minus)/(2*phi*1e-5)
        if i==0:
            Fobj_derivative=1e6
        #step = -Fobj/Fobj_derivative
        #phi = phi + step
        if phi<0:
            phi = phi - step
            phi = phi + step/1e3
        #print 'given step',phi,step,Fobj,Fobj_derivative
        i = i+1
        print L,phi,Fobj,Fobj_minus,Fobj_plus
        
        L = L+0.1e-10
        if Fobj_minus<Fobj_minus_old:
            phi = phi+0.01
            L = L-0.1e-10
        
        diffs.append(Fobj)
        Ls.append(L)
        phis.append(phi)

        Fobj_old = Fobj
        Fobj_plus_old = Fobj_minus
        Fobj_minus_old = Fobj_plus

        #if L>9e-10:
        #    phi = phi + 0.1
        #    L = 1e-10
        
        #if phi>8.0:
        #    phi = 1.0
        #    L = L + 0.1e-10
        
        #raw_input('...')
        
        
    env.append(Tv)
    env.append(rhov)
    env.append(rhol)
    env.append(Pv)
    return env
#====================================================================================== 


############TRIED TO IMPLEMENT PSO TO USE MAXWELL CONSTRUCTION====================================================

    nmap = 50
    P11 = np.empty((nmap,nmap))
    P22 = np.empty((nmap,nmap))
    u11 = np.empty((nmap,nmap))
    u22 = np.empty((nmap,nmap))
    dmap = np.empty((nmap,nmap))
    dens1 = np.empty((nmap))
    dens2 = np.empty((nmap))
    dens1 = np.linspace(rho1,rhomax,nmap)
    dens2 = np.linspace(rhomin,rho2,nmap)
    #for i in range(0,nmap):
    #    print i
    #    for j in range(0,nmap):
    #        rho1 = dens1[i]
    #        rho2 = dens2[j]
    #        P11[i][j] = InterpolatedUnivariateSpline(rho,Pspl1,k=3)(rho1)
    #        P22[i][j] = InterpolatedUnivariateSpline(rho,Pspl1,k=3)(rho2)
    #        u11[i][j] = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho1)
    #        u22[i][j] = InterpolatedUnivariateSpline(rho,uspl1,k=3)(rho2)
    #        dP = abs(P11[i][j]-P22[i][j])*1e5
    #        du = abs(u11[i][j]-u22[i][j])*1e5
    #        dmap[i][j] = dP+du

    #TRY MAXWELL'S CONSTRUCTION
    #Plin = np.linspace(Pmin,Pmax,10)
    #area = np.empty((10))
    #for i in range(0,10):
    #    print 'i',i
    #    Pint = P-Plin[i]
    #    rho1 = numerical.falsi_spline(rho,Pint,rho[0],rhomax,1e-7)
    #    print 'rho1',rho1
    #    rho2 = numerical.falsi_spline(rho,Pint,rhomax,rhomin,1e-7)
    #    print 'rho2',rho2
    #    rho3 = numerical.falsi_spline(rho,Pint,rhomin,rho[min2],1e-7)
    #    print 'rho3',rho3

    #    ind1 = False
    #    ind2 = False
    #    ind3 = False
    #    for j in range(0,10000):
    #        if rho1-rho[j]<0 and (ind1==False):
    #            rho1ind = j-1
    #            ind1 = True
    #        if rho2-rho[j]<0 and (ind2==False):
    #            rho2ind = j-1
    #            ind2 = True
    #        if rho3-rho[j]<0 and (ind3==False):
    #            rho3ind = j-1
    #            ind3 = True
    #    print 'inds',rho1ind,rho2ind,rho3ind

    #    if abs(rho1-rho2)<1e-3:
    #        area1 = 0
    #    else:
    #        area1 = numerical.trapezoidal(rho,P,rho1ind,rho2ind)
    #    if abs(rho2-rho3)<1e-3:
    #        area2 = 0
    #    else:
    #        area2 = numerical.trapezoidal(rho,P,rho2ind,rho3ind)
    #    area[i] = abs(area1-area2)
    #    print area1,area2,area[i]
    #print area

    nparticles = 2 #number of particles
    nparameter = 1 #number of parameters to adjust
    ndata = 1
    Plin = np.linspace(Pmin,Pmax,nparticles)
    argvec = []
    argvec.append(rho)
    argvec.append(P)
    argvec.append(rho[0])
    argvec.append(rhomax)
    argvec.append(rhomin)
    argvec.append(rho[min2])
    Fmin = np.empty((1))
    Fmin[0] = 1e-10
    #Optimize parameters with PSO
    P_list = []
    P_list.append(Plin)
    #P_opt = PSO.PSO_1(nparticles,nparameter,ndata,P_list,Fmin,numerical.maxwell_area,argvec,Pmax,Pmin)
    #print P_opt
#############################################################=============================================


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
