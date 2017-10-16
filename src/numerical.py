import numpy as np
import menus
import math
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline, splrep, splev, UnivariateSpline

#Calculate roots of given coefficients for volume------------------------------------------
def cubic_v_roots(coef,R,T,P):
    
    V = []
    
    coef = menus.flatten(coef)
    coef = np.array(coef)
    
    alfa = coef[0]
    beta = coef[1]
    gama = coef[2]
    
    #Newton Raphson to find first root
    #Initial guess
    V1 = R*T/P
    tol = 1e-10
    err = tol+1
    i = 1
    stop = 1.0
    while err>tol or i<3:
        V1old = V1
        F = V1**3+alfa*V1**2+beta*V1+gama
        dF = 3*V1**2+2*alfa*V1+beta
        V1 = V1 - stop*F/dF
        err = abs((V1-V1old)/V1)
        i = i+1
        if i>50:
            stop = stop/1.1
            i = 1
    
    #Solve quadratic equation to find other roots
    delta = (alfa+V1)**2-4*1*(V1**2+alfa*V1+beta)
    
    if delta>0:
        V2 = (-(alfa+V1)+(delta)**0.5)/2
        V3 = (-(alfa+V1)-(delta)**0.5)/2
    else:
        V2 = V1
        V3 = V1
    
    V.append(V1)
    V.append(V2)
    V.append(V3)
    
    return V
#==========================================================================================

#Calculate trapezoidal integral in given interval------------------------------------------
def trapezoidal(x,y,a,b):
    i = 0
    area = 0;
    if b==a:
        h = 1e-15
    else:
        h = (x[b]-x[a])/(b-a)
    area = area + y[a]
    i = a+1
    while i<b:
        area = area + 2*y[i]
        i = i+1
    area = area + y[b]
    area = area * h/2
    return area;
#==========================================================================================

#Finds maximum integer position of dPdV inside binodal curve-------------------------------
def bin_max(vec):
    vec = np.array(vec)
    size = vec.shape[0]
    i=int(0.02*size)
    while vec[i]<0:
        i = i-1
    pmax_cond=1
    while pmax_cond>0:
        pmax_cond = vec[i]
        i=i+1
        if i==size:
            i=size+1
            break
    return i-1
#==========================================================================================

#Finds minimum integer position of dPdV inside binodal curve-------------------------------
def bin_min(vec):
    vec = np.array(vec)
    size = vec.shape[0]
    i=int(0.90*size)
    pmin_cond=1
    while pmin_cond>0:
        pmin_cond = vec[i]
        i=i-1
    return i+1
#==========================================================================================

#Finds minimum integer position of dPdV inside binodal curve, given seed--------------------
def bin_min_seed(vec,seed):
    vec = np.array(vec)
    i=seed+1
    pmin_cond=-1
    while pmin_cond<0:
        pmin_cond = vec[i]
        i=i+1
    return i-1
#==========================================================================================

#Uses Regula Falsi method with cubic spline interpolation to find root in given interval---
def falsi_spline(x,y,a,b,tol):
    yc = tol+1
    it = 0
    a0 = a
    b0 = b
    #spl = splrep(x,y,k=3)
    #spl = UnivariateSpline(x,y)
    #spl.set_smoothing_factor(0.5)
    while abs(yc)>tol:
        ya = InterpolatedUnivariateSpline(x,y,k=3)(a)
        yb = InterpolatedUnivariateSpline(x,y,k=3)(b)
        c = b - yb*(a-b)/(ya-yb)
        yc = InterpolatedUnivariateSpline(x,y,k=3)(c)
        
        #ya = splev(a,spl)
        #yb = splev(b,spl)
        #c = b - yb*(a-b)/(ya-yb)
        #yc = splev(c,spl)
        
        #ya = spl(a)
        #yb = spl(b)
        #c = b - yb*(a-b)/(ya-yb)
        #yc = spl(c)
        #print c,yc,a,b,ya,yb
    
        if ya*yc<0:
            b = c
        else:
            a = c
            
        if c<a or c>b:
            c = bisect_spline(x,y,a0,b0,tol)
            yc = InterpolatedUnivariateSpline(x,y,k=3)(c)
            break
        it = it+1
    #input('...')
    return c
#==========================================================================================

#Uses Bisection method with cubic spline interpolation to find root in given interval------
def bisect_spline(x,y,a,b,tol):
    yc = tol+1
    it = 0
    a0 = a
    b0 = b
    while abs(yc)>tol:
        ya = InterpolatedUnivariateSpline(x,y,k=3)(a)
        yb = InterpolatedUnivariateSpline(x,y,k=3)(b)
        c = (a+b)/2
        
        if it>5:
            yc = InterpolatedUnivariateSpline(x,y,k=3)(c)
            
        if ya*yc<0:
            b = c
        else:
            a = c
            
        it = it+1
    return c
#==========================================================================================

#Uses Secant method with cubic spline interpolation to find root in given interval---------
def secant_spline(x,y,a,b,tol):
    xi0 = a
    xi = b
    fxi1 = tol+1
    fxi0 = InterpolatedUnivariateSpline(x,y,k=3)(xi0)
    fxi  = InterpolatedUnivariateSpline(x,y,k=3)(xi)
    while abs(fxi1)>tol:
        xi1 = xi - fxi*(xi0-xi)/(fxi0-fxi)
        fxi1 = InterpolatedUnivariateSpline(x,y,k=3)(xi1)
        xi0 = xi
        xi  = xi1
        fxi0 = fxi
        fxi = fxi1
        print fxi1,xi1
    return xi1
#==========================================================================================
