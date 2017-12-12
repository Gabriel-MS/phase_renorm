import association
import data
import envelope
import numpy as np
import PSO

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
    Psat_exp = data.loadexp2(expfile[0])[1] #Liquid Saturated Density
    Psat_var = data.loadexp2(expfile[0])[2] #Liquid Saturated Density
    dens_exp = data.loadexp2(expfile[0])[3] #Saturated Vapor Pressure
    dens_var = data.loadexp2(expfile[0])[4] #Saturated Vapor Pressure

    #Calculate Saturated Vapor Pressure and Liquid Saturated Density for given parameters
    r_data = []
    env = envelope.TV_envelope(T_exp,IDs,EoS,MR,kij,nc,AR,CR,SM,r_data)
    dens_calc = 1/np.array((env[3]))
    Psat_calc = np.array((env[0]))
    #env = envelope.PV_vec_envelope(EoS,IDs,MR,T_exp,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,0)
    #dens_calc = env[2]
    #Psat_calc = env[3]

    #Calculate Objective Function
    dens_calc = np.array(dens_calc)
    Psat_calc = np.array(Psat_calc)
    dens_exp = np.array(dens_exp)
    Psat_exp = np.array(Psat_exp)

    Fobj_dens = np.sum((dens_calc-dens_exp)**2/dens_var)
    Fobj_P    = np.sum((Psat_calc-Psat_exp)**2/Psat_var)
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
    #out.append(Tc_calc)
    #out.append(Pc_calc)
    #out.append(rhoc_calc)
    #out.append(Fobj_T)
    #out.append(Fobj_P)
    return out
#=========================================================================================

#Calculate Objective function based on Critical Temperature and Critical Pressure---------
def objFunc_dens_Psat_full(par,argss):

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
    data.modify_CPA(IDs,par)

    #Calculate new auto-association parameters
    auto = []
    auto = association.CPA_auto(AR,nc,IDs)
    en_auto = auto[0]
    beta_auto = auto[1]

    #Recover Experimental data for liquid saturated density
    T_exp    = data.loadexp2(expfile[0])[0] #Temperatures
    Psat_exp = data.loadexp2(expfile[0])[1] #Liquid Saturated Density
    Psat_var = data.loadexp2(expfile[0])[2] #Liquid Saturated Density
    dens_exp = data.loadexp2(expfile[0])[3] #Saturated Vapor Pressure
    dens_var = data.loadexp2(expfile[0])[4] #Saturated Vapor Pressure
    

    #Calculate Saturated Vapor Pressure and Liquid Saturated Density for given parameters
    r_data = []
    env = envelope.TV_envelope(T_exp,IDs,EoS,MR,kij,nc,AR,CR,SM,r_data)
    dens_calc = 1/np.array((env[3]))
    Psat_calc = np.array((env[0]))
    #env = envelope.PV_vec_envelope(EoS,IDs,MR,T_exp,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,0)
    #dens_calc = env[2]
    #Psat_calc = env[3]

    #Calculate Objective Function
    dens_calc = np.array(dens_calc)
    Psat_calc = np.array(Psat_calc)
    dens_exp = np.array(dens_exp)
    Psat_exp = np.array(Psat_exp)

    Fobj_dens = np.sum((dens_calc-dens_exp)**2/dens_var)
    Fobj_P    = np.sum((Psat_calc-Psat_exp)**2/Psat_var)
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
    #out.append(Tc_calc)
    #out.append(Pc_calc)
    #out.append(rhoc_calc)
    #out.append(Fobj_T)
    #out.append(Fobj_P)
    return out
#=========================================================================================

#Given initial T, using renormalization method, estimate association parameters----------
def Estimate_Parameters_assoc(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,expfile,estimate_bool,crit_bool,AR):

    #Parameters for PSO
    nswarm = 10
    nparameter = 2
    ndata = 1

    #Create boundaries
    bmax = np.empty((2))
    bmin = np.empty((2))
    bmin[0] = 0.02 #0.02
    bmax[0] = 0.03 #0.03
    bmin[1] = 0.01 #0.01
    bmax[1] = 0.02 #0.02

    #Create Particles
    p = np.empty((nswarm,nparameter))
    for i in range(0,nswarm):
        for j in range(0,nparameter):
            p[i][j] = np.random.uniform(bmin[j],bmax[j])
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
    PSO.PSO(nparameter,ndata,nswarm,objFunc_dens_Psat,argss,p,bmin,bmax)

    par = 1
    return par
#====================================================================================== 

#Given initial T, using renormalization method, estimate association parameters----------
def Estimate_Parameters_CPA(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,expfile,estimate_bool,crit_bool,AR):

    #Parameters for PSO
    nswarm = 15
    nparameter = 5
    ndata = 1

    #Create boundaries
    """
    bmax = np.empty((5))
    bmin = np.empty((5))
    bmin[0] = 3.9e-7
    bmax[0] = 4.1e-7
    bmin[1] = 3.0e-5
    bmax[1] = 3.2e-5
    bmin[2] = 0.4
    bmax[2] = 0.45
    bmin[3] = 0.02
    bmax[3] = 0.03
    bmin[4] = 0.01
    bmax[4] = 0.02
    """
    bmax = np.empty((5))
    bmin = np.empty((5))
    bmin[0] = 3.85e-7
    bmax[0] = 4.15e-7
    bmin[1] = 2.85e-5
    bmax[1] = 3.35e-5
    bmin[2] = 0.3
    bmax[2] = 0.5
    bmin[3] = 0.015
    bmax[3] = 0.035
    bmin[4] = 0.01
    bmax[4] = 0.03

    #Create Particles
    p = np.empty((nswarm,nparameter))
    for i in range(0,nswarm):
        for j in range(0,nparameter):
            p[i][j] = np.random.uniform(bmin[j],bmax[j])
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
    best = PSO.PSO(nparameter,ndata,nswarm,objFunc_dens_Psat_full,argss,p,bmin,bmax)
    best_param = best[0]
    param_list = best[1]
    param_Fobj = best[2]

    #Modify parameters in properties data bank, using best found solution
    data.modify_CPA(IDs,best[0])
    
    print param_list
    raw_input('...paramlist...')
    print param_Fobj
    raw_input('...Fobj...')

    return best
#====================================================================================== 
