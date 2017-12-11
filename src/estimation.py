import association
import data
import envelope
import numpy as np

#Calculate Objective function based on Critical Temperature and Critical Pressure---------
def objFunc_dens_Psat(par,argss):

    #Parameters
    EoS     = argss[0]
    IDs     = argss[1]
    MR     = argss[2]
    T      = argss[3]
    Tfinal  = argss[4]
    stepT   = argss[5]
    nd      = argss[6]
    nx      = argss[7]
    kij     = argss[8]
    nc      = argss[9]
    CR      = argss[10]
    en_auto = argss[11]
    beta_auto = argss[12]
    SM      = argss[13]
    n       = argss[14]
    estimate_bool = argss[15]
    crit_bool = argss[16]
    expfile = argss[17]
    AR      = argss[18]

    #Modify parameters in properties data bank
    data.modify_assoc(IDs,par[0],par[1])

    #Calculate new auto-association parameters
    auto = []
    auto = association.CPA_auto(AR,nc,IDs)
    en_auto = auto[0]
    beta_auto = auto[1]

    #Recover Experimental data for liquid saturated density
    T_exp    = data.loadexp2(expfile[0])[0] #Temperatures
    dens_exp = data.loadexp2(expfile[0])[1] #Liquid Saturated Density
    Psat_exp = data.loadexp2(expfile[1])[2] #Saturated Vapor Pressure

    #Calculate Saturated Vapor Pressure and Liquid Saturated Density for given parameters
    env = envelope.PV_vec_envelope(EoS,IDs,MR,T_exp,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,0)
    dens_calc = env[2]
    Psat_calc = env[3]

    #Calculate Objective Function
    Fobj_dens = np.sum(((dens_calc-dens_exp)/dens_exp)**2)
    Fobj_P    = np.sum(((Psat_calc-Psat_exp)/Psat_exp)**2)
    Fobj      = Fobj_dens + Fobj_P

    """
    print '--------------------------------'
    print 'Parameters:',L__est,phi__est
    print 'Critical Temperature:',Tc_calc,Tc_exp
    print 'Critical Pressure:',Pc_calc,Pc_exp
    print 'Objective Function:',Fobj,(Tc_calc-Tc_exp)**2/var_Tc_exp,(Pc_calc-Pc_exp)**2/var_Pc_exp
    print '--------------------------------\n'
    """

    out = []
    out.append(Fobj)
    out.append(Tc_calc)
    out.append(Pc_calc)
    out.append(rhoc_calc)
    out.append(Fobj_T)
    out.append(Fobj_P)
    return out
#=========================================================================================

#Given initial T, using renormalization method, estimate association parameters----------
def Estimate_Parameters_assoc(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,expfile,estimate_bool,crit_bool,AR):

    #Parameters for PSO
    nswarm = 20
    nparameter = 2
    ndata = 1

    #Create Particles
    p = np.empty((nswarm,nparameter))
    for i in range(0,nswarm):
        p[i][0] = np.random.uniform(1e-10,9e-10)
        p[i][1] = np.random.uniform(0.01,10)
    print 'particles'
    print p

    #Organize Parameters
    argss = []
    argss.append(EoS)
    argss.append(IDs)
    argss.append(MR)
    argss.append(T)
    argss.append(Tfinal)
    argss.append(stepT)
    argss.append(nd)
    argss.append(nx)
    argss.append(kij)
    argss.append(nc)
    argss.append(CR)
    argss.append(en_auto)
    argss.append(beta_auto)
    argss.append(SM)
    argss.append(n)
    argss.append(estimate_bool)
    argss.append(crit_bool)
    argss.append(expfile)
    argss.append(AR)

    #Initialize PSO method
    PSO.PSO(nparameter,ndata,nswarm,objFunc_dens_Psat,argss,p)

    par = 1
    return par
#====================================================================================== 
