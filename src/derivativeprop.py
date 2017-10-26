import numpy as np
import menus

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = 8.314462175e-6 #m3.MPa/K/mol

#Function to calculate derivative properties of pure compounds after renormalization--------------- 
def calc_isothermal_dev_prop_pure(T,f,P,rho,h):
    
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
    
    n = len(P1[0])
    
    #Derivatives
    d2AdT2 = (A2 - 2*A1 + A0)/(h**2)
    
    dAdT = (A2 - A0)/(2*h)

    Pv = P1[0]
    rhov = rho1[0]
    Vv = V1[0]

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

    d2PdT2 = (P2 - 2*P1 + P0)/(h**2)

    #Isochoric Heat Capacity
    Cv = -T1*d2AdT2
    #Cv = -R*T1*T1*d2AdT2-2*R*T1*dAdT
    Cv = Cv[0]
    
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
    #Cp = Cv - T1[0]*d2PdT2[0]/dPdV-R

    #Speed of Sound
    Mw = 0.05812
    w = (Cp/Cv*dPdrho/Mw*1e6)**(0.5)
    #w = (-(Vv**2)*Cp/Cv*dPdV/Mw)**(0.5)
    
    deriv_data = []
    deriv_data.append(Pv)
    deriv_data.append(rhov)
    deriv_data.append(d2AdT2[0])
    deriv_data.append(dAdT[0])
    deriv_data.append(Cv)
    deriv_data.append(dPdrho)
    deriv_data.append(inv_kT)
    deriv_data.append(kT)
    deriv_data.append(lnkT)
    deriv_data.append(dPdT[0])
    deriv_data.append(alfa)
    deriv_data.append(uJT)
    deriv_data.append(Cp)
    deriv_data.append(w)
    
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
    ax[0,0].set_xlim(xmin=0, xmax=25)
    ax[0,0].set_ylabel('Cv')
    ax[0,0].set_title('Cv')
    
    ax[0,1].plot(P_plot,lnkT_plot)
    ax[0,1].set_xlim(xmin=0, xmax=25)
    ax[0,1].set_ylabel('lnkT')
    ax[0,1].set_title('lnkT')
    
    ax[1,0].plot(P_plot,alfa_plot)
    ax[1,0].set_xlim(xmin=0, xmax=25)
    ax[1,0].set_ylabel('alfa')
    ax[1,0].set_title('alfa')
    
    ax[1,1].plot(P_plot,w_plot)
    ax[1,1].set_xlim(xmin=0, xmax=25)
    ax[1,1].set_ylim(ymin=0, ymax=500)
    ax[1,1].set_ylabel('w')
    ax[1,1].set_title('w')
    
    ax[2,0].plot(P_plot,uJT_plot)
    ax[2,0].set_xlim(xmin=0, xmax=25)
    ax[2,0].set_ylim(ymin=0, ymax=10)
    ax[2,0].set_ylabel('uJT')
    ax[2,0].set_title('uJT')
    
    ax[2,1].plot(P_plot,Cp_plot)
    ax[2,1].set_xlim(xmin=0, xmax=25)
    ax[2,1].set_ylabel('Cp')
    ax[2,1].set_title('Cp')
    
    fig.savefig('../output/%s' %figname)
#==================================================================================================

#Function to report calculated derivative properties of pure compounds after renormalization------- 
def report_isothermal_dev_prop_pure(title,P,rho,d2AdT2,dAdT,Cv,dPdrho,inv_kT,kT,lnkT,dPdT,alfa,uJT,Cp,w,print_options):
    n = len(P)
    header = 'P(MPa);rho(mol/m3);Cv;Cp;alfa;uJT;w;lnkT;inv_kT;kT;d2AdT2;dAdT;dPdrho;dPdT;\n'
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
                str(round(kT[i],9)),str(round(d2AdT2[i],9)),str(round(dAdT[i],9)),str(round(dPdrho[i],9)),str(round(dPdT[i],9))]
                lin = ('%s\n' % ';'.join(lin1))
                file.write(lin)
#==================================================================================================
