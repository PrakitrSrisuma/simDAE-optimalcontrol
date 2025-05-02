function outputs = obj_tempchange_pl(Tb, tb, tspan, dTdt_obj, T_ini, option, input)

    Heater_Profile = [Tb, tb];
    [t,T] = ode15s(@(t,T)PDE_2Phases(t,T,Heater_Profile,input), tspan, T_ini, option);
   
    outputs = cal_dTdt_SSE(t,cal_Tavg(T(:,input.n1+1:end-1)),dTdt_obj);

end