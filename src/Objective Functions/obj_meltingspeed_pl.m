function outputs = obj_meltingspeed_pl(Tb, tb, tspan, dSdt_obj, T_ini, option, input)

    Heater_Profile = [Tb, tb];
    [t,T] = ode15s(@(t,T)PDE_2Phases(t,T,Heater_Profile,input), tspan, T_ini, option);
    
    outputs = cal_dSdt_SSE(t,T(:,end),dSdt_obj);

end