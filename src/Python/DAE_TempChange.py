def DAE_TempChange(dTdt_target,Tb0):
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

        
    # Define Gekko variables    
    m = GEKKO()
    Ts = m.Array(m.Var,n1);
    Tl = m.Array(m.Var,n1);
    S = m.Var(S0);
    Tb = m.Var(Tb0,fixed_initial=False);
    Tavg = m.Var(0);


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
    m.Equation(S.dt()== (Ste/S)*(1/dR)*(-2*Ts[n1-1]+Ts[n1-2]/2) - (k0*Ste/(1-S))*(1/dR)*(2*Tl[0]-Tl[1]/2))
    m.Equation(Tavg == m.sum(Tl)/n1)
    m.Equation(Tavg.dt() == dTdt_target)

    # DAE solver
    tspan = np.linspace(0,input.endtime,60)
    m.time = tspan;
    m.options.IMODE=7
    m.options.NODES=3
    m.options.RTOL = 1e-7
    m.options.OTOL = 1e-7


    # Time measurement
    tic = time.time()
    m.solve(disp=False)
    toc = time.time()
    twall_DAE = toc - tic


    Tb_opt = np.array(Tb)
    return [twall_DAE, Tb_opt, tspan]

output_MATLAB = DAE_TempChange(x,y)

##plt.plot(m.time, Tb.value, 'b')
##plt.show()
