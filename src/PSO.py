import numpy as np
import math
import numerical

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def PSO_1(n,npar,ndata,param,ydata,func,args,max_bound,min_bound):

    #n = swarm size
    #npar = number of parameters
    #ndata = size of data to fit

    #Initialize parameters
    k = 1
    kmax = 10000
    c1 = 0.5
    c2 = 1.0
    w = 0.05
    Fobj = 0
    tol = 1e-20
    
    obj_v = np.empty((n))
    
    px = []
    for j in range(0,len(param)):
        px.append(param[j])
    
    best_swarm_pos = np.empty((len(param)))
    best_swarm_obj = 5000.0
    best_part_pos = np.empty((npar,n))
    best_obj_v = np.empty((n))
    
    v = np.empty((npar,n))
    for i in range(0,n):
        for j in range(0,npar):
            v[j][i] = np.random.rand(1,n)[0][i]
    
    #Initialize best solution
    for j in range(0,npar):
        best_part_pos[j] = px[j][0]
        
    for j in range(0,len(param)):
        best_swarm_pos[j] = px[j][0] #Select random values for global best to begin
    
    for i in range(0,n):
        best_obj_v[i] = 5000.0
    
    #MAIN LOOP==================================================================
    while (k<kmax) and min(best_obj_v)>tol:
        #Evaluate positions and objective function
        for i in range(0,n):
            Fobj = 0
            for j in range(0,ndata):
               ycalc = func(args,param[j][i])
               Fobj = Fobj + abs(ydata[j]-ycalc)
            obj_v[i] = Fobj
            print 'Fobj',obj_v[i],i
    
        #Update particles own best position
        for i in range(0,n):
            if obj_v[i]<=best_obj_v[i]:
                for j in range(0,npar):
                    best_part_pos[j][i] = px[j][i] #best position parameter 1
                best_obj_v[i] = obj_v[i]
    
        #Update best solution among all particles
        for i in range(0,n):
            if obj_v[i]<=best_swarm_obj:
                best_swarm_obj = obj_v[i]
                for j in range(0,npar):
                    best_swarm_pos[j] = px[j][np.argmin(best_obj_v)] #best position parameter 1
    
        #update positions
        for i in range(0,n):
            for j in range(0,npar):
                v[j][i] = w*v[j][i] + c1*np.random.rand()*(best_part_pos[j][i]-px[j][i]) + c2*np.random.rand()*(best_swarm_pos[j]-px[j][i])
                px[j][i] = px[j][i] + v[j][i]

        #check boundaries
        for i in range(0,n):
            for j in range(0,npar):
                if px[j][i]>max_bound or px[j][i]<max_bound:
                    px[j][i] = np.random.rand()*(max_bound-min_bound)+min_bound
            
        k = k+1
        print k,best_swarm_pos,best_obj_v[i],px
        
    #print best_swarm_pos,min(best_obj_v),k
    print 'iter',k
    return best_swarm_pos
    #===========================================================================
