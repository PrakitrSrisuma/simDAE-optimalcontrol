def DAE_MeltingSpeed_Mixed(dSdt_target,Tb0,tspan,nodes,IC):
    from gekko import GEKKO
    import numpy as np
    import matplotlib.pyplot as plt
    import CellThawing_InputParameters as input
    import time


    n1 = input.n1;
    dR = input.dR;
    Ste = input.Ste;
    Tm = input.Tm_d;
    alp = input.alp;
    k1 = input.k1;
    k2 = input.k2;
    k0 = k2/k1;
    U = input.U_d;
    S0 = 1-input.tol;
    IC = np.array(IC);


    # Define Gekko variables
    m = GEKKO()
##    Ts = [m.Var(IC[0:n1]) for i in range(n1)]
##    Tl = [m.Var(IC[n1:2*n1]) for i in range(n1)]
    Ts = m.Array(m.Var,n1);
    Tl = m.Array(m.Var,n1);
    for i in range(0,n1):
        Ts[i].value = IC[i]
        Tl[i].value = IC[i+n1]
        
##    Ts = m.Var(value = IC[0:n1-1]);
##    Tl = m.Var(value = IC[n1:2*n1-1]);
    S = m.Var(IC[2*n1]);
    Tb = m.Var(Tb0);

    # ODEs for solid 
    for i in range(0,n1):
        if i == 0:
            m.Equation(Ts[i].dt()==(4/S**2)*(Ts[i+1]-Ts[i])/dR**2)
        elif i == n1-1:
            m.Equation(Ts[i].dt()==(1/S**2)*(0-2*Ts[i]+Ts[i-1])/dR**2 + (1/(S**2*i*dR))*(0-Ts[i-1])/(2*dR) \
                                   + (i*dR/S)*(0-Ts[i-1])*S.dt()/(2*dR))
        else:
            m.Equation(Ts[i].dt()==(1/S**2)*(Ts[i+1]-2*Ts[i]+Ts[i-1])/dR**2 + (1/(S**2*i*dR))*(Ts[i+1]-Ts[i-1])/(2*dR) \
                                   + (i*dR/S)*(Ts[i+1]-Ts[i-1])*S.dt()/(2*dR))

            

    # ODEs for liquid
    for i in range(0,n1):
        if i == 0:
            m.Equation(Tl[i].dt()==(alp/(1-S)**2)*(Tl[i+1]-2*Tl[i]+0)/(dR**2) + (alp/((1-S)*(S+1*dR*(1-S))))*(Tl[i+1]-0)/(2*dR) \
                                   + ((1-1*dR)/(1-S))*(Tl[i+1]-0)*S.dt()/(2*dR))
        elif i == n1-1:
            m.Equation(Tl[i].dt()==(alp/(1-S)**2)*((U*(1-S)*(2*dR)*(Tb-Tl[i])+Tl[i-1])-2*Tl[i]+Tl[i-1])/(dR**2) + \
                                  (alp/((1-S)*(S+(i+1)*dR*(1-S))))*((U*(1-S)*(2*dR)*(Tb-Tl[i])+Tl[i-1])-Tl[i-1])/(2*dR) \
                                   + ((1-(i+1)*dR)/(1-S))*((U*(1-S)*(2*dR)*(Tb-Tl[i])+Tl[i-1])-Tl[i-1])*S.dt()/(2*dR))
        else:
            m.Equation(Tl[i].dt()==(alp/(1-S)**2)*(Tl[i+1]-2*Tl[i]+Tl[i-1])/(dR**2) + (alp/((1-S)*(S+(i+1)*dR*(1-S))))*(Tl[i+1]-Tl[i-1])/(2*dR) \
                                   + ((1-(i+1)*dR)/(1-S))*(Tl[i+1]-Tl[i-1])*S.dt()/(2*dR))


    # DAEs
    m.Equation(dSdt_target == (Ste/S)*(1/dR)*(-2*Ts[n1-1]+Ts[n1-2]/2) - (k0*Ste/(1-S))*(1/dR)*(2*Tl[0]-Tl[1]/2))
    m.Equation(S.dt()== dSdt_target)

    # DAE solver
##    tspan = np.linspace(0,endtime,20)
    m.time = tspan 
    m.options.IMODE= 7
    m.options.NODES= nodes
    m.options.RTOL = 1e-7
    m.options.OTOL = 1e-7

    # Time measurement
    tic = time.time()
    m.solve(disp=False)
    toc = time.time()
    twall_DAE = toc - tic


    # Data export to MATLAB
##    Tb_opt = np.array(Tb_tmp)
    Tb_opt = np.array(Tb)
##    Tb_opt[Tb_opt>1] = 1
    return [twall_DAE, Tb_opt, tspan]

output_MATLAB = DAE_MeltingSpeed_Mixed(z1,z2,z3,z4,z5)

##plt.plot(m.time, Tb.value, 'b')
##plt.show()
