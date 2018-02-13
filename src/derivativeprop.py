import numpy as np
import menus
import correlations
import data
import eos

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

R = 8.314462175e-6 #m3.MPa/K/mol

#Function to calculate derivative properties of pure compounds after renormalization--------------- 
def calc_isothermal_dev_prop_pure(T,f,P,rho,h,IDs):
    
    f0 = np.array(menus.flatten(f[0]))
    f1 = np.array(menus.flatten(f[1]))
    f2 = np.array(menus.flatten(f[2]))
    
    T0 = np.array(menus.flatten(T[0]))
    T1 = np.array(menus.flatten(T[1]))
    T2 = np.array(menus.flatten(T[2]))
    
    rho0 = np.array(menus.flatten(rho[0]))
    rho1 = np.array(menus.flatten(rho[1]))
    rho2 = np.array(menus.flatten(rho[2]))
    
    P0 = np.array(menus.flatten(P[0]))
    P1 = np.array(menus.flatten(P[1]))
    P2 = np.array(menus.flatten(P[2]))
    
    V0 = 1/rho0
    V1 = 1/rho1
    V2 = 1/rho2
    
    A0 = f0/rho0
    A1 = f1/rho1
    A2 = f2/rho2

    F0 = A0/(R*T0)
    F1 = A1/(R*T1)
    F2 = A2/(R*T2)
    
    n = len(P1[0])
    
    #Derivatives
    d2AdT2 = (A2 - 2*A1 + A0)/(h**2)
    d2FdT2 = (F2 - 2*F1 + F0)/(h**2)
    
    dAdT = (A2 - A0)/(2*h)
    dFdT = (F2 - F0)/(2*h)

    Pv = P1[0]
    rhov = rho1[0]
    Vv = V1[0]
    Fv = F1[0]

    dPdrho = np.empty((n))
    for i in range(1,n-1):
        dPdrho[i] = (Pv[i+1]-Pv[i-1])/(rhov[i+1]-rhov[i-1])
    dPdrho[0] = (Pv[1]-Pv[0])/(rhov[1]-rhov[0])
    dPdrho[n-1] = (Pv[n-1]-Pv[n-2])/(rhov[n-1]-rhov[n-2])

    dPdV = np.empty((n))
    for i in range(1,n-1):
        dPdV[i] = (Pv[i+1]-Pv[i-1])/(Vv[i+1]-Vv[i-1])
    dPdV[0] = (Pv[1]-Pv[0])/(Vv[1]-Vv[0])
    dPdV[n-1] = (Pv[n-1]-Pv[n-2])/(Vv[n-1]-Vv[n-2])

    d2PdV2 = np.empty((n))
    for i in range(1,n-1):
        d2PdV2[i] = (Pv[i+1]+2*Pv[i]-Pv[i-1])/(Vv[i+1]-Vv[i-1])
    d2PdV2[0] = (Pv[1]-Pv[0])/((Vv[1]-Vv[0])**2)
    d2PdV2[n-1] = (Pv[n-1]-Pv[n-2])/((Vv[n-1]-Vv[n-2])**2)

    d2FdV2 = np.empty((n))
    for i in range(1,n-1):
        d2FdV2[i] = (Fv[i+1]+2*Fv[i]-Fv[i-1])/(Vv[i+1]-Vv[i-1])
    d2FdV2[0] = (Fv[1]-Fv[0])/((Vv[1]-Vv[0])**2)
    d2FdV2[n-1] = (Fv[n-1]-Fv[n-2])/((Vv[n-1]-Vv[n-2])**2)

    #d2FdTV = np.empty((n))
    #for i in range(1,n-1):
    #    d2FdTV[i] = (dFdT[i+1]-dFdT[i-1])/(Vv[i+1]-Vv[i-1])
    #d2FdTV[0] = (dFdT[1]-dFdT[0])/(Vv[1]-Vv[0])
    #d2FdTV[n-1] = (dFdT[n-1]-dFdT[n-2])/(Vv[n-1]-Vv[n-2])

    d2PdT2 = (P2 - 2*P1 + P0)/(h**2)
    dPdT = (P2 - P0)/(2*h)

    #dPdV = -R*T*d2FdV2-R*T/(Vv**2)
    #dPdT = -R*T*d2FdTV+Pv/T1[0]

    #Isochoric Heat Capacity
    Cv = -T1*d2AdT2
    #Cv = -R*T1*T1*d2FdT2-2*R*T1*dFdT
    Cv = Cv[0]
    print T1
    
    #Isothermal compression coefficient
    inv_kT = rhov*dPdrho
    kT = 1/inv_kT
    lnkT = np.log(kT)
    
    #Thermal Expansion coefficient
    dPdT = (P2 - P0)/(2*h)
    alfa = kT*dPdT
    alfa = alfa[0]
    
    #Joule-Thomson coefficient
    uJT = T1[0]*dPdT-rhov*dPdrho
    uJT = uJT[0]
    
    #Isobaric Heat Capacity
    Cp = Cv + T1[0]*(alfa**2)/kT/rhov
    #Cp = Cv - T1[0]*(dPdT[0]**2)/dPdV-R

    #Speed of Sound
    Mw = data.mass(IDs)[0]
    Cp_ig = correlations.ideal_cp(IDs,T1[0])
    Cp1 = Cp + Cp_ig
    Cv_ig = Cp_ig - R
    Cv1 = Cv + Cv_ig
    w = (Cp1/Cv1*dPdrho/Mw*1e6)**(0.5)
    #w = (-(Vv**2)*Cp/Cv*dPdV/Mw)**(0.5)
    
    #for i in range(0,len(f1[0])):
    #    print rho1[0][i],f1[0][i],'wow'
    
    deriv_data = []
    deriv_data.append(Pv)
    deriv_data.append(rhov)
    deriv_data.append(d2AdT2[0])
    deriv_data.append(dAdT[0])
    deriv_data.append(Cv)
    deriv_data.append(dPdrho)
    deriv_data.append(inv_kT) #inv_kT
    deriv_data.append(kT) #kT
    deriv_data.append(lnkT) #lnkT
    deriv_data.append(dPdT[0])
    deriv_data.append(alfa)
    deriv_data.append(uJT)
    deriv_data.append(Cp)
    deriv_data.append(w) #w
    deriv_data.append(Cp1) #Cp1
    deriv_data.append(Cv1) #Cv1
    
    return deriv_data
#==================================================================================================

#Function to calculate derivative properties of pure compounds after renormalization--------------- 
def calc_isothermal_dev_prop_pure_analytical(T,A,P,V,IDs,EoS,MR,kij):
    
    b = eos.b_calc(IDs,EoS)
    x = np.array([0.999,0.001])
    bmix = eos.bmix_calc(MR,b,x)
    a = eos.a_calc(IDs,EoS,T)
    amix = eos.amix_calc(MR,a,x,kij)
    ac = eos.ac_calc(IDs,EoS,T)[0]
    m = eos.kapa_calc(IDs,EoS)[0]
    Tc = np.array(data.Tc(IDs))[0]
    rho = 1/V
    
    d2adT2 = ac*m*(1+m)*np.sqrt(Tc/T)/(2*T*Tc)
    B = bmix*P/(R*T)
    Z = P*V/(R*T)
    
    f = A/V
    
    noll = np.ones((len(V)))
    
    Cv_res = (T*d2adT2/(bmix*np.sqrt(8)))*np.log((Z+B*(1+np.sqrt(2)))/(Z+B*(1-np.sqrt(2))))
    
    d2AdT2 = -d2adT2/bmix/np.sqrt(8)*np.log((1+rho*bmix*(1+np.sqrt(2)))/(1+rho*bmix*(1-np.sqrt(2))))
    
    h = 5e-2
    T_p = T+T*h
    T_m = T-T*h
    
    a_p = eos.a_calc(IDs,EoS,T_p)
    amix_p = eos.amix_calc(MR,a_p,x,kij)

    a_m = eos.a_calc(IDs,EoS,T_m)
    amix_m = eos.amix_calc(MR,a_m,x,kij)
    
    A = R*T*np.log(rho/(1-rho*bmix))-amix/bmix/np.sqrt(8)*np.log((1+rho*bmix*(1+np.sqrt(2)))/(1+rho*bmix*(1-np.sqrt(2))))
    A_p = R*T_p*np.log(rho/(1-rho*bmix))-amix_p/bmix/np.sqrt(8)*np.log((1+rho*bmix*(1+np.sqrt(2)))/(1+rho*bmix*(1-np.sqrt(2))))
    A_m = R*T_m*np.log(rho/(1-rho*bmix))-amix_m/bmix/np.sqrt(8)*np.log((1+rho*bmix*(1+np.sqrt(2)))/(1+rho*bmix*(1-np.sqrt(2))))
        
    d2AdT2_2 = (A_p-2*A+A_m)/((T*h)**2)
    
    A = A*1e8
    A_p = A_p*1e8
    A_m = A_m*1e8
    d2AdT2 = d2AdT2*1e6
    d2AdT2_2 = d2AdT2_2*1e6
    
    deriv_data = []
    deriv_data.append(P)
    deriv_data.append(rho)
    deriv_data.append(A)
    deriv_data.append(A_p)
    deriv_data.append(Cv_res)
    deriv_data.append(A_m)
    deriv_data.append(d2AdT2_2) 
    deriv_data.append(d2AdT2)
    deriv_data.append(noll)
    deriv_data.append(noll)
    deriv_data.append(noll)
    deriv_data.append(noll)
    deriv_data.append(noll)
    deriv_data.append(noll)
    deriv_data.append(noll)
    deriv_data.append(noll)
    
    return deriv_data
#==================================================================================================

#Function to plot calculated derivative properties of pure compounds after renormalization--------- 
def plot_isothermal_dev_prop_pure(P,rho,d2AdT2,dAdT,Cv,dPdrho,inv_kT,kT,lnkT,dPdT,alfa,uJT,Cp,w,print_options,figname):
    
    n = len(P)
    
    #plot selection
    alfa_plot = np.empty((n-10))
    P_plot = np.empty((n-10))
    rho_plot = np.empty((n-10))
    Cv_plot = np.empty((n-10))
    Cp_plot = np.empty((n-10))
    lnkT_plot = np.empty((n-10))
    w_plot = np.empty((n-10))
    uJT_plot = np.empty((n-10))
    
    for i in range(5,n-5):
        j = i-5
        alfa_plot[j] = alfa[i]
        P_plot[j] = P[i]
        rho_plot[j] = rho[i]
        Cv_plot[j] = Cv[i]
        Cp_plot[j] = Cp[i]
        lnkT_plot[j] = lnkT[i]
        w_plot[j] = w[i]
        uJT_plot[j] = uJT[i]
    
    #plot variables
    fig, ax = plt.subplots(3,2,figsize=(10,10))
    
    ax[0,0].plot(P_plot,Cv_plot)
    ax[0,0].set_xlim(xmin=0, xmax=10.107)
    ax[0,0].set_ylabel('Cv')
    ax[0,0].set_title('Cv')
    
    ax[0,1].plot(P_plot,lnkT_plot)
    ax[0,1].set_xlim(xmin=0, xmax=10.107)
    ax[0,1].set_ylabel('lnkT')
    ax[0,1].set_title('lnkT')
    
    ax[1,0].plot(P_plot,alfa_plot)
    ax[1,0].set_xlim(xmin=0, xmax=10.107)
    ax[1,0].set_ylabel('alfa')
    ax[1,0].set_title('alfa')
    
    ax[1,1].plot(P_plot,w_plot)
    ax[1,1].set_xlim(xmin=0, xmax=10.107)
    ax[1,1].set_ylim(ymin=70, ymax=590)
    ax[1,1].set_ylabel('w')
    ax[1,1].set_title('w')
    
    ax[2,0].plot(P_plot,uJT_plot)
    ax[2,0].set_xlim(xmin=0, xmax=10.107)
    ax[2,0].set_ylim(ymin=0, ymax=10)
    ax[2,0].set_ylabel('uJT')
    ax[2,0].set_title('uJT')
    
    ax[2,1].plot(P_plot,Cp_plot)
    ax[2,1].set_xlim(xmin=0, xmax=10.107)
    ax[2,1].set_ylabel('Cp')
    ax[2,1].set_title('Cp')
    
    fig.savefig('../output/%s' %figname)
#==================================================================================================

#Function to report calculated derivative properties of pure compounds after renormalization------- 
def report_isothermal_dev_prop_pure(title,P,rho,d2AdT2,dAdT,Cv,dPdrho,inv_kT,kT,lnkT,dPdT,alfa,uJT,Cp,w,Cp1,Cv1,print_options):
    n = len(P)
    header = 'P(MPa);rho(mol/m3);Cv_res;Cp_res;alfa;uJT;w;lnkT;inv_kT;kT;d2AdT2;dAdT;dPdrho;dPdT;Cv;Cp\n'
    savedir = str('../output/%s' %title)
    with open(savedir,'w') as file:
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
        for i in range(0,n):
                lin1 = [str(round(P[i],9)),str(round(rho[i],9)),str(round(Cv[i],9)),str(round(Cp[i],9)),
                str(round(alfa[i],9)),str(round(uJT[i],9)),str(round(w[i],9)),str(round(lnkT[i],9)),str(round(inv_kT[i],9)),
                str(round(kT[i],9)),str(round(d2AdT2[i],9)),str(round(dAdT[i],9)),str(round(dPdrho[i],9)),str(round(dPdT[i],9)),
                str(round(Cv1[i],9)),str(round(Cp1[i],9))]
                lin = ('%s\n' % ';'.join(lin1))
                file.write(lin)
#==================================================================================================
