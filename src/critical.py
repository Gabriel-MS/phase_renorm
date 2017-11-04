#Algorithms to deal with crtical point calculation
import numpy as np
import menus


import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

from scipy.interpolate import RectBivariateSpline

#Hicks-Young Algorithm to calculate mixture critical point, modified by Sadus-----------
def sadus_hicks_young(f,rho,x,bm):

    Vb = 1/rho[0]

    nx = len(x)
    nv = len(Vb)
    A      = np.empty((nx,nv))
    fmat   = np.empty((nx,nv))
    d2AdV2 = np.empty((nx,nv))
    dAdV   = np.empty((nx,nv))
    d2Adx2 = np.empty((nx,nv))
    dAdx   = np.empty((nx,nv))
    d2AdVx = np.empty((nx,nv))
    d2AdxV = np.empty((nx,nv))
    dWdx   = np.empty((nx,nv))
    dWdV   = np.empty((nx,nv))
    W      = np.empty((nx,nv))
    Q      = np.empty((nx,nv))
    V      = np.empty((nx,nv))

    for i in range(0,nx):
        for j in range(0,nv):
            fmat[i][j] = f[i][j]
            V[i][j] = Vb[j]*bm[i]

    for i in range(0,nx):
        for j in range(0,nv):
            A[i][j] = fmat[i][j]*V[i][j]

    dAdx   = np.gradient(A,x,edge_order=2,axis=0)
    d2Adx2 = np.gradient(dAdx,x,edge_order=2,axis=0)
    dAdV   = np.gradient(A,Vb,edge_order=2,axis=1)
    for i in range(0,nx):
        for j in range(0,nv):
            dAdV[i][j] = dAdV[i][j]/bm[i]
    d2AdV2 = np.gradient(dAdV,Vb,edge_order=2,axis=1)
    for i in range(0,nx):
        for j in range(0,nv):
            d2AdV2[i][j] = d2AdV2[i][j]/bm[i]
    d2AdVx = np.gradient(dAdV,x,edge_order=2,axis=0)
    d2AdxV = np.gradient(dAdx,Vb,edge_order=2,axis=1)
    for i in range(0,nx):
        for j in range(0,nv):
            d2AdxV[i][j] = d2AdxV[i][j]/bm[i]
    W      = (d2Adx2*d2AdV2)-(d2AdVx**2)
    dWdx   = np.gradient(W,x,edge_order=2,axis=0)
    dWdV   = np.gradient(W,Vb,edge_order=2,axis=1)
    for i in range(0,nx):
        for j in range(0,nv):
            dWdV[i][j] = dWdV[i][j]/bm[i]
    A3x    = np.gradient(d2Adx2,x,edge_order=2,axis=0)
    A3V    = np.gradient(d2AdV2,Vb,edge_order=2,axis=1)
    for i in range(0,nx):
        for j in range(0,nv):
            A3V[i][j] = A3V[i][j]/bm[i]
    AV2x   = np.gradient(d2AdVx,x,edge_order=2,axis=0)
    A2Vx   = np.gradient(d2AdV2,x,edge_order=2,axis=0)
    
    for i in range(0,nx):
        for j in range(0,nv):
            Q[i][j] = A3x[i][j] - 3*AV2x[i][j]*(d2AdVx[i][j]/d2AdV2[i][j]) + 3*A2Vx[i][j]*((d2AdVx[i][j]/d2AdV2[i][j])**2) - A3V[i][j]*((d2AdVx[i][j]/d2AdV2[i][j])**3)
    Q = A3x - 3*AV2x*(d2AdVx/d2AdV2) + 3*A2Vx*((d2AdVx/d2AdV2)**2) - A3V*((d2AdVx/d2AdV2)**3)
    #Q = dWdV*d2Adx2-dWdx*d2AdxV
    dQdx   = np.gradient(Q,x,edge_order=2,axis=0)
    dQdV   = np.gradient(Q,Vb,edge_order=2,axis=1)
    for i in range(0,nx):
        for j in range(0,nv):
            dQdV[i][j] = dQdV[i][j]/bm[i]

    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #xplot, Vbplot = np.meshgrid(x, Vb)
    #surf = ax.plot_surface(xplot, Vbplot, W, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #ax.set_zlim(-0.5, 0.5)x
    #ax.zaxis.set_major_locator(LinearLocator(10))
    #ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    #plt.show()

    #Find zeros of W and Q
    #invert volume
    V1 = np.empty((nv))
    W1 = np.empty((nx,nv))
    Q1 = np.empty((nx,nv))
    dWdx1 = np.empty((nx,nv))
    dWdV1 = np.empty((nx,nv))
    dQdx1 = np.empty((nx,nv))
    dQdV1 = np.empty((nx,nv))

    for j in range(0,nv):
        V1[j] = Vb[j]
    for j in range(0,nv):
        Vb[j] = V1[nv-j-1]

    for i in range(0,nx):
        for j in range(0,nv):
            W1[i][j] = W[i][j]
            Q1[i][j] = Q[i][j]
            dWdx1[i][j] = dWdx[i][j]
            dWdV1[i][j] = dWdV[i][j]
            dQdx1[i][j] = dQdx[i][j]
            dQdV1[i][j] = dQdV[i][j]
    for i in range(0,nx):
        for j in range(0,nv):
            W[i][j] = W1[i][nv-j-1]
            Q[i][j] = Q1[i][nv-j-1]
            dWdx[i][j] = dWdx1[i][nv-j-1]
            dWdV[i][j] = dWdV1[i][nv-j-1]
            dQdx[i][j] = dQdx1[i][nv-j-1]
            dQdV[i][j] = dQdV1[i][nv-j-1]

    W2 = RectBivariateSpline(x,Vb,W)
    Q2 = RectBivariateSpline(x,Vb,Q)
    Wx = RectBivariateSpline(x,Vb,dWdx)
    Qx = RectBivariateSpline(x,Vb,dQdx)
    Wv = RectBivariateSpline(x,Vb,dWdV)
    Qv = RectBivariateSpline(x,Vb,dQdV)

    tolx = 1e-2
    tolV = 1e-2
    xint = x[nx/2]
    xint = 0.875
    Vint = Vb[nv-10]
    xintold = xint
    Vintold = Vint
    dx = tolx+1
    dV = tolV+1
    Wint = tolx+1
    Qint = tolV+1
    k = 0
    while (abs(Wint)>tolx) or (abs(Qint)>tolV):
        Wint = W2(xint,Vint)[0][0]
        Qint = Q2(xint,Vint)[0][0]
        Wxint = Wx(xint,Vint)[0][0]
        Wvint = Wv(xint,Vint)[0][0]
        Qxint = Qx(xint,Vint)[0][0]
        Qvint = Qv(xint,Vint)[0][0]

        J = np.matrix([[Wxint,Wvint],
                      [Qxint,Qvint]])

        F = np.matrix([[-Wint],
                      [-Qint]])

        dx = (J.I*F)[0]
        dV = (J.I*F)[1]

        if dx>0.1:
            dx = 0.01
        if dx<-0.1:
            dx = -0.01
        if dx<1e-8 and dx>=0:
            dx = dx*1000

        xint = xintold + dx
        Vint = Vintold + dV

        xintold = xint
        Vintold = Vint

        #if Vint>Vb[nv-1]:
        #    Vint = Vb[5]
        #if Vint<0:
        #    Vint = Vb[nv/2]
 
        if k%2000==0:
            print dx,dV,xint,Vint,Wint,Qint
        k = k+1
        #raw_input('---')


    print Wint,Qint

    #print W
    #print '----------------------------------------------------------------------------'
    #print Q
#======================================================================================




