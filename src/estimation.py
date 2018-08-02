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
    out.append(Fobj_dens)
    out.append(Fobj_P)
    return out
#=========================================================================================

#Calculate Objective function based on saturated liquid density and saturation pressure---
def objFunc_dens_Psat_crit(par,argss):

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
    
    #Particle parameters
    L__est = par[0]
    phi__est = par[1]
    
    #Modify renormalization parameters
    data.modify_renorm(IDs,par)

    #Recover Experimental data for liquid saturated density
    T_exp    = data.loadexp3(expfile[0])[0] #Temperatures
    Psat_exp = data.loadexp3(expfile[0])[1] #Saturated Vapor Pressure
    #Psat_var = data.loadexp2(expfile[0])[2] #Liquid Saturated Density
    dens_liq_exp = data.loadexp3(expfile[0])[2] #Liquid Saturated Density
    dens_vap_exp = data.loadexp3(expfile[0])[3] #Vapor Saturated Density
    #dens_var = data.loadexp2(expfile[0])[4] #Saturated Vapor Pressure
    

    #Calculate Saturated Vapor Pressure and Liquid Saturated Density for given parameters
    r_data = []
    #env = envelope.TV_envelope(T_exp,IDs,EoS,MR,kij,nc,AR,CR,SM,r_data)
    #dens_calc = 1/np.array((env[3]))
    #Psat_calc = np.array((env[0]))
    env = envelope.PV_vec_envelope(EoS,IDs,MR,T_exp,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n)
    dens_liq_calc = env[2]
    dens_vap_calc = env[1]
    Psat_calc = env[3]

    #Calculate Objective Function
    dens_liq_calc = np.array(dens_liq_calc)
    dens_vap_calc = np.array(dens_vap_calc)
    Psat_calc = np.array(Psat_calc)
    dens_liq_exp = np.array(dens_liq_exp)
    dens_vap_exp = np.array(dens_vap_exp)
    Psat_exp = np.array(Psat_exp)

    Fobj_dens_liq = np.sum(((dens_liq_calc-dens_liq_exp)/dens_liq_exp)**2)
    Fobj_dens_vap = np.sum(((dens_vap_calc-dens_vap_exp)/dens_vap_exp)**2)
    Fobj_P    = np.sum(((Psat_calc-Psat_exp)/Psat_exp)**2)
    #for i in range(0,len(dens_liq_calc)):
        #if dens_liq_calc[i]<0.1:
        #    Fobj_dens_liq = 100.0
        #if dens_vap_calc[i]<0.1:
        #    Fobj_dens_vap = 100.0
    Fobj      = Fobj_dens_liq + Fobj_dens_vap + 5*Fobj_P
    Fobj      = Fobj_dens_liq + Fobj_P
    
    #Calculate critical point
    finalT = 800.0
    crit_pt = envelope.PV_findTc3_envelope(EoS,IDs,MR,T,finalT,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,True,0,0)
    Tc_calc = crit_pt[0][len(crit_pt[0])-1]
    Pc_calc = crit_pt[3][len(crit_pt[3])-1]
    rhoc_calc = crit_pt[1][len(crit_pt[1])-1]
    
    data.surface_out(L__est,phi__est,Fobj_dens_liq,Fobj_dens_vap,Fobj_P,Tc_calc,Pc_calc,rhoc_calc)

    print '--------------------------------'
    print 'Parameters:',L__est,phi__est
    print 'Critical Point:',Tc_calc,Pc_calc,rhoc_calc
    print 'Objective Function:',Fobj,Fobj_dens_liq,Fobj_dens_vap,Fobj_P
    print '--------------------------------\n'
    
    out = []
    out.append(Fobj)
    #out.append(Tc_calc)
    #out.append(Pc_calc)
    #out.append(rhoc_calc)
    out.append(Fobj_dens_liq)
    out.append(Fobj_dens_vap)
    out.append(Fobj_P)
    return out
#=========================================================================================

#Calculate Objective function based on critical point---
def objFunc_dens_Psat_Tc(par,argss):

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

    #Particle parameters
    L__est = par[0]
    phi__est = par[1]

    #Modify renormalization parameters
    data.modify_renorm(IDs,par)
    
    #Recover Experimental data for critical point
    T_exp    = data.loadexp3(expfile[0])[0] #Temperatures
    Psat_exp = data.loadexp3(expfile[0])[1] #Saturated Vapor Pressure
    dens_liq_exp = data.loadexp3(expfile[0])[2] #Liquid Saturated Density
    
    Tc_exp = T_exp[len(T_exp)-1]
    Pc_exp = Psat_exp[len(Psat_exp)-1]
    rhoc_exp = dens_liq_exp[len(dens_liq_exp)-1]
    
    #Calculates critical point
    Crit = envelope.PV_findTc2_envelope(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
    Tc_calc = Crit[0]
    Pc_calc = Crit[3]
    rhoc_calc = Crit[1]

    Fobj_Tc = abs((Tc_calc-Tc_exp)/Tc_exp)
    Fobj_Pc = abs((Pc_calc-Pc_exp)/Pc_exp)
    Fobj_rhoc    = abs((rhoc_calc-rhoc_exp)/rhoc_exp)
    Fobj      = 10*Fobj_Tc + 10*Fobj_Pc + 10*Fobj_rhoc
    
    print '--------------------------------'
    print 'Parameters:',L__est,phi__est
    print 'Critical Point:',Tc_calc,Pc_calc,rhoc_calc
    print 'Critical deltas:',abs(Tc_calc-Tc_exp),abs(Pc_calc-Pc_exp),abs(rhoc_calc-rhoc_exp)
    print 'Objective Function:',Fobj,Fobj_Tc,Fobj_Pc,Fobj_rhoc
    print '--------------------------------\n'
    
    out = []
    out.append(Fobj)
    #out.append(Tc_calc)
    #out.append(Pc_calc)
    #out.append(rhoc_calc)
    out.append(Fobj_Tc)
    out.append(Fobj_Pc)
    out.append(Fobj_rhoc)
    return out
#=========================================================================================

#Calculate Objective function based on derivative properties------------------------------
def objFunc_dens_deriv(par,argss):

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

    #Particle parameters
    L__est = par[0]
    phi__est = par[1]

    #Modify renormalization parameters
    data.modify_renorm(IDs,par)
    
    #Recover Experimental data for critical point
    T_exp    = data.loadexp3(expfile[0])[0] #Temperatures
    Psat_exp = data.loadexp3(expfile[0])[1] #Saturated Vapor Pressure
    dens_liq_exp = data.loadexp3(expfile[0])[2] #Liquid Saturated Density
    
    Tc_exp = T_exp[len(T_exp)-1]
    Pc_exp = Psat_exp[len(Psat_exp)-1]
    rhoc_exp = dens_liq_exp[len(dens_liq_exp)-1]
    
    #Calculates critical point
    Crit = envelope.PV_findTc2_envelope(EoS,IDs,MR,Tc_exp,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,False,0,0)
    Tc_calc = Crit[0]
    Pc_calc = Crit[3]
    rhoc_calc = Crit[1]
    
    #Calculates speed of sound at critical point
    dp_dat = envelope.PV_deriv_calc_envelope(EoS,IDs,MR,Tc_exp,Tc_exp,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n)

    Fobj_Tc = abs((Tc_calc-Tc_exp)/Tc_exp)
    Fobj_Pc = abs((Pc_calc-Pc_exp)/Pc_exp)
    Fobj_rhoc    = abs((rhoc_calc-rhoc_exp)/rhoc_exp)
    Fobj      = 10*Fobj_Tc + 10*Fobj_Pc + 10*Fobj_rhoc
    
    print '--------------------------------'
    print 'Parameters:',L__est,phi__est
    print 'Critical Point:',Tc_calc,Pc_calc,rhoc_calc
    print 'Critical deltas:',abs(Tc_calc-Tc_exp),abs(Pc_calc-Pc_exp),abs(rhoc_calc-rhoc_exp)
    print 'Objective Function:',Fobj,Fobj_Tc,Fobj_Pc,Fobj_rhoc
    print '--------------------------------\n'
    
    out = []
    out.append(Fobj)
    #out.append(Tc_calc)
    #out.append(Pc_calc)
    #out.append(rhoc_calc)
    out.append(Fobj_Tc)
    out.append(Fobj_Pc)
    out.append(Fobj_rhoc)
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
    nswarm = 10
    nparameter = 2
    ny = 2
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
    bmin[0] = 4E-10
    bmax[0] = 6E-10
    bmin[1] = 0.01
    bmax[1] = 10.0

    #Create Particles
    p = np.empty((nswarm,nparameter))
    for i in range(0,nswarm):
        for j in range(0,nparameter):
            p[i][j] = np.random.uniform(bmin[j],bmax[j])
    #print 'particles'
    #print p
    
    """
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
    
    perc = 0.001
    #4.04904923118e-07;3.07860333899e-05;0.578720592387;0.0232048954443;0.0208354177678
    mini = np.empty((nparameter))
    safe = np.empty((nparameter))
    mini[0] = 4.0521947012852E-07
    mini[1] = 0.0000305941689368634
    mini[2] = 0.351256542099425
    mini[3] = 0.025651688835633
    mini[4] = 0.0117773203079985
    k = 0
    l = 0
    q = 0
    kmax = 1e5
    """
    
    """
    #Calculate Confidence Region
    p[0] = mini
    while k<kmax:
        #for i in range(0,nparameter):
            #p[0][i] = mini[i]*(1.000+perc*np.random.rand()-perc*np.random.rand())
            #p[0][i] = np.random.uniform(bmin[i],bmax[i])
        Fobj = objFunc_dens_Psat_full(p[0],argss)
        if Fobj[0]<=4.090685:
            envelope.report_param([p[0]],'../output/confidence_region.csv')
            l = l+1
            safe = p[0]
            for j in range(0,nparameter):
                p[0][j] = p[0][j]*(1+np.random.rand()*perc-np.random.rand()*perc)
        else:
            for j in range(0,nparameter):
                p[0][j] = mini[j]*(1+np.random.rand()*perc-np.random.rand()*perc)
        k = k+1
        print k,l
    
    
    #Covariance Matrix Calculation
    pp = np.empty((nswarm,nparameter))
    pn = np.empty((nswarm,nparameter))
    p[0] = mini
    pp[0] = p[0]
    pn[0] = p[0]
    T_exp    = data.loadexp2(expfile[0])[0] #Temperatures
    Psat_var = data.loadexp2(expfile[0])[2] #Liquid Saturated Density
    dens_var = data.loadexp2(expfile[0])[4] #Saturated Vapor Pressure
    nexp = len(T_exp)
    r_data = []
    env = []
    
    #envmax = envelope.TV_envelope(T_exp,IDs,EoS,MR,kij,nc,AR,CR,SM,r_data)
    #input('...')
    
    #W matrix calculation
    W = np.zeros((ny*nexp,ny*nexp))
    var = np.hstack((Psat_var,dens_var))
    np.fill_diagonal(W,var)
    W = np.linalg.inv(W)
    
    #J matrix calculation
    J = np.empty((ny*nexp,nparameter))
    Jt = np.empty((nparameter,ny*nexp))
    for i in range(0,nparameter):
        
        #-----POSITIVE PERTURBATION---
        #Modify parameter i
        pp[0][i] = p[0][i]*(1+1e-3)
        
        #Modify parameters in properties data bank
        data.modify_CPA(IDs,pp[0])

        #Calculate new auto-association parameters
        auto = []
        auto = association.CPA_auto(AR,nc,IDs)
        en_auto = auto[0]
        beta_auto = auto[1]
        
        #Calculate answers
        envmax = envelope.TV_envelope(T_exp,IDs,EoS,MR,kij,nc,AR,CR,SM,r_data)
        dens_calcmax = 1/np.array((envmax[3]))
        Psat_calcmax = np.array((envmax[0]))
        
        #-----NEGATIVE PERTURBATION---
        #Modify parameter i
        pn[0][i] = p[0][i]*(1-1e-3)
        
        #Modify parameters in properties data bank
        data.modify_CPA(IDs,pn[0])

        #Calculate new auto-association parameters
        auto = []
        auto = association.CPA_auto(AR,nc,IDs)
        en_auto = auto[0]
        beta_auto = auto[1]
        
        #Calculate answers
        envmin = envelope.TV_envelope(T_exp,IDs,EoS,MR,kij,nc,AR,CR,SM,r_data)
        dens_calcmin = 1/np.array((envmin[3]))
        Psat_calcmin = np.array((envmin[0]))
        
        #Calculate derivatives
        der1 = (Psat_calcmax-Psat_calcmin)/(pp[0][i]-pn[0][i])
        der2 = (dens_calcmax-dens_calcmin)/(pp[0][i]-pn[0][i])
        Jt[i] = np.hstack((der1,der2))
        
        #Back to original, to modify only the next parameter
        pp[0] = p[0]
        pn[0] = p[0]
        print i
    
    #Covariance matrix final
    J = np.transpose(Jt)
    cov1 = np.dot( (np.dot(Jt,W)) , J)
    covinv = np.linalg.inv(cov1)
    covE = objFunc_dens_Psat_full(p[0],argss)[0]/(ny*nexp-nparameter)
    desvE = covE**0.5
    
    COV = covE*covinv
    
    print '----COVARIANCE---'        
    print COV
    print '================='
        
    #Correlation matrix
    CORR = np.zeros((nparameter,nparameter))
    for i in range(0,nparameter):
        for j in range(0,nparameter):
            CORR[i][j] = COV[i][j]/(((COV[i][i])**0.5)*((COV[j][j])**0.5))
    
    print '----CORRELATION---'        
    print CORR
    print '=================='
    
    #Intervalo das respostas estimadas
    #Calculate answers
    env = envelope.TV_envelope(T_exp,IDs,EoS,MR,kij,nc,AR,CR,SM,r_data)
    dens_calc = 1/np.array((env[3]))
    Psat_calc = np.array((env[0]))
    
    t_inf = 2.064
    t_sup = 2.064
    
    print t_inf*(np.diag(covinv)**0.5)*(covE**0.5)
    
    Psat_inf = np.empty((nexp))
    Psat_sup = np.empty((nexp))
    dens_inf = np.empty((nexp))
    dens_sup = np.empty((nexp))
    
    cov_dens = objFunc_dens_Psat_full(p[0],argss)[1]/(ny*nexp-nparameter)
    cov_P    = objFunc_dens_Psat_full(p[0],argss)[2]/(ny*nexp-nparameter)
    
    print cov_dens
    
    for i in range(0,nexp):
        j = i+nexp
        
        Jti = np.transpose(J[i])
        Psat_inf[i] = Psat_calc[i]-t_inf*((np.dot( (np.dot(Jti,covinv)) , J[i]))**0.5)*(covE**0.5)
        Psat_sup[i] = Psat_calc[i]+t_sup*((np.dot( (np.dot(Jti,covinv)) , J[i]))**0.5)*(covE**0.5)
        
        Jtj = np.transpose(J[j])
        dens_inf[i] = dens_calc[i]-t_inf*((np.dot( (np.dot(Jtj,covinv)) , J[j]))**0.5)*(covE**0.5)
        dens_sup[i] = dens_calc[i]+t_sup*((np.dot( (np.dot(Jtj,covinv)) , J[j]))**0.5)*(covE**0.5)

    print Psat_inf
    print Psat_sup
    print dens_inf
    print dens_sup
    print Psat_calc
    print dens_calc
    """
    #Initialize PSO method
    #best = PSO.PSO(nparameter,ndata,nswarm,objFunc_dens_Psat_full,argss,p,bmin,bmax)
    #best_param = best[0]
    #param_list = best[1]
    #param_Fobj = best[2]

    #Modify parameters in properties data bank, using best found solution
    data.modify_CPA(IDs,best[0])
    
    print param_list
    raw_input('...paramlist...')
    print param_Fobj
    raw_input('...Fobj...')

    return best
#====================================================================================== 

#Given initial T, using renormalization method, estimate renormalization parameters----------
def Estimate_Parameters_crit(EoS,IDs,MR,T,Tfinal,stepT,nd,nx,kij,nc,CR,en_auto,beta_auto,SM,n,expfile,estimate_bool,crit_bool,AR):

    #Parameters for PSO
    nswarm = 20
    nparameter = 2
    ny = 2
    ndata = 1

    #Create boundaries
    bmax = np.empty((2))
    bmin = np.empty((2))
    bmin[0] = 5.00E-10
    bmax[0] = 6.00E-10
    bmin[1] = 0.50
    bmax[1] = 1.50
    
    bounds = np.zeros((2,nparameter))
    bounds[0][0] = bmin[0] #min x
    bounds[0][1] = bmin[1] #min y
    bounds[1][0] = bmax[0] #max x
    bounds[1][1] = bmax[1] #max y
    
    #Organize Parameters
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

    #Create Particles
    p = np.empty((nswarm,nparameter))
    for i in range(0,nswarm):
        for j in range(0,nparameter):
            p[i][j] = np.random.uniform(bmin[j],bmax[j])
            #p[i][1] = 1.0
    
    """
    #Create Particles to calculate surface for objective function
    ntot = 50
    p = np.empty((ntot*ntot,nparameter))
    Lrange = np.linspace(bmin[0],bmax[0],ntot)
    phirange = np.linspace(bmin[1],bmax[1],ntot)
    for i in range(0,ntot):
        for j in range(0,ntot):
            k = i*ntot+j
            p[k][0] = Lrange[i]
            p[k][1] = phirange[j]
    nswarm = ntot*ntot
    """
    
    #p[0][0] = 5.45e-10
    #p[0][1] = 1.50
    
    
    #Initialize PSO method - ALL PARAMETERS
    best = PSO.PSO(nparameter,ndata,nswarm,objFunc_dens_Psat_crit,argss,p,bmin,bmax)
    best_param = best[0]
    param_list = best[1]
    param_Fobj = best[2]
    
    #Initialize MDPSO
    #best = PSO.MDPSO(nparameter,ndata,nswarm,objFunc_dens_Psat_crit,argss,p,bounds)
    #best_param = best[0]
    #param_list = best[1]
    #param_Fobj = best[2]
    
    #Initialize PSO method to fit parameter 1 (L)
    #best = PSO.PSO(nparameter,ndata,nswarm,objFunc_dens_Psat_Tc,argss,p,bmin,bmax)
    #best_param = best[0]
    #param_list = best[1]
    #param_Fobj = best[2]

    #Initialize PSO method to fit parameter 2 (phi)
    #Create Particles to calculate surface for objective function
    ntot = 20
    nparameter = 2
    p = np.empty((ntot,nparameter))
    phirange = np.linspace(bmin[1],bmax[1],ntot)
    for i in range(0,ntot):
        p[i][0] = 1e10
        p[i][1] = phirange[i]
    nswarm = ntot
    print 'particles'
    print p
    best = PSO.PSO(nparameter,ndata,nswarm,objFunc_dens_Psat_crit,argss,p,bmin,bmax)
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
