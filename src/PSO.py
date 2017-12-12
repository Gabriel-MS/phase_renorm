import numpy as np
import math
import numerical

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

#Standard PSO routine-------------------------------------------------------------------------------
def PSO(nparam,ndata,nswarm,objFunc,args,p,bmin,bmax):
    #nparam  - number of parameters
    #ndata   - number of data points
    #nswarm  - number of particles
    #objFunc - Objective Function
    #args    - Objective Function arguments

    #Organize parameters and arguments------------------------------------
    xdata = args
    #=====================================================================

    #Initialize PSO parameters--------------------------------------------
    k = 1           #iteration counter
    kmax = 100    #max iterations allowed
    c1 = 0.5        #weight of particle's best position distance
    c2 = 1.0        #weight of swarm's best position distance
    w = 0.5        #weight of particle's last velocity
    tol = 1e-4     #tolerance
    #=====================================================================

    #Initialize solutions-------------------------------------------------
    best_swarm_pos    = np.array(p[0])             #Best Swarm position, choosing "randomly the first particle"
    best_swarm_obj    = np.array(objFunc(p[0],xdata))    #Objective Function of Best Swarm position
    best_particle_obj = np.empty((nswarm))          #Initialize Best Particle objective function
    best_particle_pos = np.empty((nswarm,nparam))          #Initialize Best Particle position
    Fobj              = np.empty((nswarm))          #Initialize objective function
    v                 = np.empty((nswarm,nparam))   #Initial velocities
    scale = np.ones((nparam))                       #Parameters scale
    for j in range(0,nparam):
        scale[j] = numerical.magnitude(p[0][j])

    for i in range(0,nswarm):
        for j in range(0,nparam):
            v[i][j]          = np.random.uniform(0,1)*(1*(10**(scale[j])))/100
    #print '---best---'
    #for i in range(0,nswarm):
    #    best_particle_obj[i] = np.array(objFunc(p[i],xdata))[0] #Objective Function of Best Particle position
    #    best_particle_pos[i] = np.array(p[i])          #Objective Function of Best Particle position
    #    print i,best_particle_obj[i]
    #print '----------'
    Fobj_old = best_particle_obj
    #=====================================================================

    #MAIN LOOP------------------------------------------------------------
    #Calculate Objective function for all particles
    while (k<kmax) and (best_swarm_obj>tol):
        
        #Calculate Objective Function for all particles
        for i in range(0,nswarm):
            #print 'part i',p[i]
            Fobj[i] = objFunc(p[i],xdata)[0]
            if k==1:
                best_particle_obj[i] = Fobj[i] #Objective Function of Best Particle position
                best_particle_pos[i] = np.array(p[i])          #Objective Function of Best Particle position
            #print 'i',i,Fobj[i]
            #Update each particle best position
            if Fobj[i]<=Fobj_old[i]:
                best_particle_obj[i] = Fobj[i]
                for j in range(0,nparam):
                    best_particle_pos[i][j] = p[i][j]
                
            #Update swarm best position
            if Fobj[i]<=best_swarm_obj:
                best_swarm_obj = Fobj[i]
                for j in range(0,nparam):
                    best_swarm_pos[j] = p[i][j]
                
        Fobj_old = Fobj
            
        #Update positions
        for i in range(0,nswarm):
            for j in range(0,nparam):
                v[i][j] = w*v[i][j] + c1*np.random.rand()*(best_particle_pos[i][j]-p[i][j]) + c2*np.random.rand()*(best_swarm_pos[j]-p[i][j])
                p[i][j] = p[i][j] + v[i][j]
        
        #Check boundaries
        for i in range(0,nswarm):
            for j in range(0,nparam):
                if p[i][j]>bmax[j]:
                    p[i][j]=bmax[j]
                if p[i][j]<bmin[j]:
                    p[i][j]=bmin[j]
        
        #Update iteration counter
        k = k+1
        print 'k',k-1,best_swarm_pos,best_swarm_obj
    #=====================================================================
    
    for i in range(0,nswarm):
        print p[i],Fobj[i]
    return best_swarm_pos  
#===================================================================================================

#Lennard-Jones PSO routine-------------------------------------------------------------------------------
def LJPSO(nparam,ndata,nswarm,objFunc,args,p):
    #nparam  - number of parameters
    #ndata   - number of data points
    #nswarm  - number of particles
    #objFunc - Objective Function
    #args    - Objective Function arguments
    
    #Organize parameters and arguments------------------------------------
    xdata = args
    #=====================================================================   
    
    #Initialize PSO parameters--------------------------------------------
    k = 1           #iteration counter
    kmax = 10000     #max iterations allowed
    c1 = 1.5        #weight of particle's best position distance
    c2 = 2.5        #weight of swarm's best position distance
    w = 0.5         #weight of particle's last velocity
    tol = 1e-5     #tolerance
    eps = 1.0       #epsilon LJ parameter
    sig = 1.0       #sigma LJ parameter
    rc  = 1.0       #cutoff radius
    m   = np.linspace(1.0,1.0,nswarm) #mass of particles
    flagcount = 0   #Number of stationary particles
    scale = np.ones((nparam)) #Parameters scale
    for i in range(0,nparam):
        scale[i] = numerical.magnitude(p[0][i])
    #====================================================================
    
    #Initialize solutions-------------------------------------------------
    best_swarm_pos    = np.array(p[0])             #Best Swarm position, choosing "randomly the first particle"
    best_swarm_obj    = np.array(objFunc(p[0],xdata))    #Objective Function of Best Swarm position
    best_particle_obj = np.empty((nswarm))          #Initialize Best Particle objective function
    best_particle_pos = np.empty((nswarm,nparam))          #Initialize Best Particle position
    Fobj              = np.empty((nswarm))          #Initialize objective function
    v                 = np.empty((nswarm,nparam))   #Initial velocities
    a                 = np.zeros((nswarm,nparam))   #Initial accelerations
    flag              = np.empty((nswarm))
    for i in range(0,nswarm):
        for j in range(0,nparam):
            v[i][j]          = np.random.uniform(0,1)*(1*(10**(scale[j]))) 
        flag[i]       = False
    for i in range(0,nswarm):
        best_particle_obj[i] = abs(np.array(objFunc(p[i],xdata))) #Objective Function of Best Particle position
        best_particle_pos[i] = np.array(p[i])          #Objective Function of Best Particle position
    Fobj_old = best_particle_obj
    #=====================================================================
    
    #MAIN LOOP------------------------------------------------------------
    #Calculate Objective function for all particles
    #while (k<kmax) and (np.amax(best_particle_obj)>tol):
    while (k<kmax) and (flagcount<=nswarm/2):
        
        #Calculate Objective Function for all particles
        for i in range(0,nswarm):
            if flag[i]==False:
                Fobj[i] = objFunc(p[i],xdata)[0]
                Fobj[i] = abs(Fobj[i])
            
                #Update each particle best position
                if Fobj[i]<=Fobj_old[i]:
                    best_particle_obj[i] = Fobj[i]
                    for j in range(0,nparam):
                        best_particle_pos[i][j] = p[i][j]
                    if best_particle_obj[i]<tol:
                        flag[i] = True
                        p[i] = best_particle_pos[i]
                        flagcount = flagcount+1
                
                #Update swarm best position
                if Fobj[i]<=best_swarm_obj:
                    best_swarm_obj = Fobj[i]
                    for j in range(0,nparam):
                        best_swarm_pos[j] = p[i][j]
                
                #New position is the best in the swarm, get swarm near to this newly converged particle
                if flag[i]==True:
                    for j in range(0,nparam):
                        best_swarm_pos[j] = p[i][j]
                
            Fobj_old = Fobj
        
        #First particle gets straight to global best
        #for j in range(0,nparam):
        #    p[0][j] = best_swarm_pos[j]


        #Update positions
        for i in range(0,nswarm):
            if flag[i]==False:
                for j in range(0,nparam):
                    v[i][j] = w*v[i][j] + c1*np.random.rand()*(best_particle_pos[i][j]-p[i][j]) + c2*np.random.rand()*(best_swarm_pos[j]-p[i][j])
                    p[i][j] = p[i][j] + v[i][j]

        #Forces and acceleration calculation if distance < cut-off
        F = np.zeros((nswarm,nparam))
        for i in range(0,nswarm):
            for l in range(i+1,nswarm):
                if flag[i]==True and flag[l]==True:
                    flag[i] = True
                else:
                    d = 0
                    for j in range(0,nparam):
                        d = d + (p[i][j]-p[l][j])**2
                    rij = d**0.5
                    #print rij,i,l
                    if rij<rc: #distance between particles is less than cut-off radius
                        for j in range(0,nparam):
                            Fij = np.random.rand()*3.14/rc*math.cos(3.14*rij/rc)*(p[i][j]-p[l][j])
                            F[i][j] = F[i][j] + Fij
                            F[l][j] = F[l][j] - Fij

        #Update velocities according to forces
        #print '=========='
        for i in range(0,nswarm):
            if flag[i]==False:
                for j in range(0,nparam):
                    a[i][j] = F[i][j]/m[i]
                #    p[i][j] = p[i][j] + a[i][j]*(1*(10**(scale[j])))   #Should we add it or subtract it???
                #print a[i][j],i,j

        #Plot actual status
        h = k
        if h%1==0:
            xa = np.empty((nswarm))
            ya = np.empty((nswarm))
            xc = np.empty((flagcount))
            yc = np.empty((flagcount))
            qq = 0
            ww = 0
            for ww in range(0,nswarm):
                if flag[ww]==True:
                    xc[qq] = best_particle_pos[ww][0]
                    yc[qq] = best_particle_pos[ww][1]
                    qq = qq+1
                else:
                    xa[ww] = best_particle_pos[ww][0]
                    ya[ww] = best_particle_pos[ww][1]
            fig = plt.figure()
            plt.plot(xa,ya,'r.')
            plt.plot(xc,yc,'g^')
            plt.xlim(1e-10,9e-10)
            plt.ylim(0,10)
            plt.savefig('xyk.png')
        
        #Update iteration counter
        k = k+1       
        print 'k',k,best_swarm_pos,best_swarm_obj,np.amax(best_particle_pos),flagcount
        #raw_input('...')
    #=====================================================================
    
    
    for i in range(0,nswarm):
        print best_particle_pos[i],best_particle_obj[i]
    
    #plot solution
    x = np.empty((flagcount))
    y = np.empty((flagcount))
    k = 0
    for i in range(0,nswarm):
        if flag[i]==True:
            x[k] = best_particle_pos[i][0]
            y[k] = best_particle_pos[i][1]
            k = k+1
    fig = plt.figure()
    plt.plot(x,y,'r.')
    plt.savefig('../output/xy.png')

    pos = []
    pos.append(x)
    pos.append(y)

    return pos
#===================================================================================================
